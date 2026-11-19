package net.maizegenetics.tassel.io

import kotlinx.io.Buffer
import kotlinx.io.Sink
import kotlinx.io.Source
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.taxa.TaxaList

/**
 * BCF 2.2 genotype files, genotype (`GT`) field only.
 *
 * [read] builds the same [GenotypeTable] as [Vcf.read]: the same taxa, sites, alleles, and calls. Phase is dropped.
 * A [kotlinx.io.Source] is inflated with [decompressing] when the bytes are gzip or BGZF, and read as raw BCF when
 * the magic is at byte 0. A [RandomAccessSource] is BGZF, inflated in parallel by [BgzfReader]. Records may cross
 * block boundaries. The next record starts at the next inflated byte, including when htslib put the header in its
 * own block.
 *
 * [write] emits the records [Vcf.write] prints as BCF 2.2 inside BGZF stored blocks. TASSEL 5 has no BCF writer.
 *
 * [BcfRead.warnings] carries a missing terminal BGZF EOF block, the same text [Vcf.read] records on [VcfRead].
 *
 * @see Vcf.read
 */
object Bcf {
    /**
     * Reads a BCF file from [source].
     *
     * Gzip and BGZF are inflated by [decompressing]. Any other bytes are raw BCF, magic at byte 0. [source] is
     * left open.
     */
    suspend fun read(source: Source, progress: ProgressReporter = ProgressReporter.None): BcfRead {
        val decompressed = source.decompressing()
        val window = InflatedWindow { destination, offset, length ->
            decompressed.source.readAtMostTo(destination, offset, offset + length)
        }
        val table = readRecords(window, progress)
        return BcfRead(table, decompressed.warnings.toList())
    }

    /**
     * Reads a bgzipped BCF from [source].
     *
     * Inflation uses [BgzfReader], so blocks decompress on more than one worker when the host has them.
     * [BcfRead.warnings] carries a missing terminal BGZF EOF block. Raw BCF is not a BGZF source; that path
     * stays on [kotlinx.io.Source].
     */
    suspend fun read(source: RandomAccessSource, progress: ProgressReporter = ProgressReporter.None): BcfRead {
        val reader = BgzfReader(source, progress = progress)
        try {
            val window = InflatedWindow { destination, offset, length ->
                reader.read(destination, offset, length)
            }
            val table = readRecords(window, progress)
            return BcfRead(table, reader.warnings.toList())
        } finally {
            reader.close()
        }
    }

    /**
     * Writes [table] as BCF 2.2 inside BGZF stored blocks.
     *
     * The header is the lines [Vcf.write] emits, plus the PASS filter line, `IDX` on each FORMAT and INFO line,
     * and one `##contig` line per distinct chromosome in site order. Each record uses that site's REF, ALT, and
     * GT indexes. QUAL is missing, FILTER is PASS, INFO is empty, and GT is the only FORMAT field. TASSEL 5 has
     * no BCF writer.
     *
     * Sites are visited through [GenotypeTable.siteBlocks]. [sink] receives the `.bcf` and is left open.
     *
     * @see net.maizegenetics.dna.snp.ExportUtils.writeToVCF
     */
    suspend fun write(table: GenotypeTable, sink: Sink, progress: ProgressReporter = ProgressReporter.None) {
        val header = bcfExportHeader(table)
        val encoder = BcfRecordEncoder(table.numberOfTaxa, header.gtId, header.contigs)
        val genotypes = ByteArray(table.numberOfTaxa)
        val numberOfSites = table.numberOfSites
        BgzfWriter(sink, BgzfLayout.Bcf).use { writer ->
            writer.write(header.preamble)
            writer.endHeader()
            for (block in table.siteBlocks()) {
                for (site in block.sites) {
                    val position = table.positions[site]
                    block.copySiteGenotypes(site, genotypes)
                    val record = vcfSiteRecord(position.knownVariants, block.siteStats(site), genotypes)
                    encoder.encode(position, record)
                    writer.write(encoder.bytes, 0, encoder.size)
                }
                progress.report(block.sites.last + 1, numberOfSites, "Writing BCF")
            }
        }
    }

    private suspend fun readRecords(window: InflatedWindow, progress: ProgressReporter): GenotypeTable {
        val header = window.readHeader()
        val sites = VcfSites(TaxaList.of(header.vcf.taxaNames), progress, "Reading BCF")
        while (true) {
            val position = window.readRecord(header, sites.calls) ?: break
            sites.add(position)
        }
        return sites.finish()
    }
}

/** A BCF file and the warnings found while inflating it. Raw BCF leaves [warnings] empty. */
class BcfRead(val table: GenotypeTable, val warnings: List<String>)

/**
 * Inflated BCF bytes, filled on demand.
 *
 * [readBcfHeader] and [readBcfRecord] pull a whole header or record from a synchronous [Source]. This buffer
 * has that record in hand before the call, so a [BgzfReader] can stay `suspend` and a record can cross a BGZF
 * block without the decoder knowing. The fill is asked for [CHUNK_BYTES] at a time, which is what drives a
 * missing terminal EOF block onto the reader's warning list.
 */
private class InflatedWindow(private val fill: suspend (ByteArray, Int, Int) -> Int) {
    private val pending = Buffer()
    private val chunk = ByteArray(CHUNK_BYTES)
    private var exhausted = false

    suspend fun readHeader(): BcfHeader {
        val preamble = (BCF_MAGIC_LENGTH + UINT32_BYTES).toLong()
        pullUntil(preamble)
        if (pending.size >= preamble) {
            val textLength = peek(preamble.toInt()).uint32Le(BCF_MAGIC_LENGTH)
            if (textLength in 1L..Int.MAX_VALUE.toLong()) pullUntil(preamble + textLength)
        }
        return readBcfHeader(pending)
    }

    suspend fun readRecord(header: BcfHeader, calls: ByteArray): Position? {
        pullUntil(1)
        if (pending.size == 0L) return null
        val prefix = BCF_RECORD_PREFIX_BYTES.toLong()
        pullUntil(prefix)
        if (pending.size < prefix) return readBcfRecord(pending, header, calls)
        val lengths = peek(BCF_RECORD_PREFIX_BYTES)
        val body = bcfBodyLength(lengths.uint32Le(0), lengths.uint32Le(BCF_RECORD_PREFIX_BYTES / 2))
        pullUntil(prefix + body)
        return readBcfRecord(pending, header, calls)
    }

    private suspend fun pullUntil(needed: Long) {
        while (pending.size < needed && !exhausted) {
            val count = fill(chunk, 0, chunk.size)
            if (count <= 0) exhausted = true else pending.write(chunk, 0, count)
        }
    }

    private fun peek(count: Int): ByteArray = pending.peek().use { it.readByteArray(count) }
}

/** Bytes pulled from the inflated stream at a time. A record longer than this is carried across fills. */
private const val CHUNK_BYTES = 64 * 1024
