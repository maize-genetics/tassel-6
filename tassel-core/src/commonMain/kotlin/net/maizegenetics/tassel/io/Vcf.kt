package net.maizegenetics.tassel.io

import kotlinx.io.Source
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList

/**
 * VCF genotype files, genotype (`GT`) field only.
 *
 * [read] takes the variants of a site from its `REF` and `ALT` columns and the calls from the first `FORMAT` field,
 * which TASSEL 5 also requires to be `GT`. Both phased (`0|1`) and unphased (`0/1`) diploid calls are read, and the
 * two are indistinguishable afterwards, because the genotype byte has no phase bit. Every other `FORMAT` field,
 * including `AD`, `DP`, `GQ`, and `PL`, is skipped, so read depth does not survive (tassel-6#14 scope).
 *
 * [read] of a [kotlinx.io.Source] does not decompress: hand it an already-inflated source. [read] of a
 * [RandomAccessSource] inflates BGZF itself, in parallel when the host has more than one worker.
 * A region of an indexed file is [IndexedVcf].
 *
 * @see net.maizegenetics.dna.snp.io.BuilderFromVCF
 */
object Vcf {
    /** The column names TASSEL 5 looks for in the `#CHROM` header line, before the first taxon. */
    val FIXED_COLUMNS: List<String> = listOf(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
    )

    /**
     * Reads a VCF file.
     *
     * Lines starting with `##` do not become sites. `##SAMPLE=` taxa annotations are read and dropped: `Taxon`
     * carries only a name (the maintainer's decision on tassel-6#14, recorded in `docs/porting-ledger.md`).
     * `##contig` ids and lengths are kept on the header from [readHeader], for a later indexed read. The next line
     * is the `#CHROM` header, which names the taxa, and every line after it that does not start with `#` is one site.
     *
     * Sites that are not in [PositionList] order are sorted with their genotypes, as [HapMap.read] does. TASSEL 5
     * instead refuses an unsorted VCF unless the caller asks for `buildAndSortInMemory()`.
     *
     * @throws IllegalArgumentException if the file is empty, the header is missing a column the reader needs, a
     *   position is not an integer, `GT` is not the first `FORMAT` field, a call names an allele the site does not
     *   list, or a line has the wrong number of genotypes.
     */
    suspend fun read(source: Source, progress: ProgressReporter = ProgressReporter.None): GenotypeTable {
        val header = readHeader(source)
        return readRecords(header, progress) { window, offset, length ->
            source.readAtMostTo(window, offset, offset + length)
        }
    }

    /**
     * Reads a bgzipped VCF from [source].
     *
     * Inflation uses [BgzfReader], so blocks decompress on more than one worker when the host has them.
     * [VcfRead.warnings] carries a missing terminal BGZF EOF block. Plain gzip is not a BGZF source;
     * that path stays on [kotlinx.io.Source].
     */
    suspend fun read(source: RandomAccessSource, progress: ProgressReporter = ProgressReporter.None): VcfRead {
        val reader = BgzfReader(source, progress = progress)
        try {
            val header = readVcfHeaderLines { reader.readLine() }
            val table = readRecords(header, progress) { window, offset, length ->
                reader.read(window, offset, length)
            }
            return VcfRead(table, reader.warnings.toList())
        } finally {
            reader.close()
        }
    }

    private suspend fun readRecords(
        header: VcfHeader,
        progress: ProgressReporter,
        fill: suspend (window: ByteArray, offset: Int, length: Int) -> Int,
    ): GenotypeTable {
        val sites = VcfSites(TaxaList.of(header.taxaNames), progress)
        val parser = VcfRecordParser(header)

        var window = ByteArray(RECORD_BUFFER_BYTES)
        var start = 0
        var end = 0
        var eof = false
        while (true) {
            when (val parsed = parser.tryParse(window, start, end, eof, sites.calls)) {
                is VcfParse.End -> break
                is VcfParse.Skip -> start += parsed.consumed
                is VcfParse.Site -> {
                    start += parsed.consumed
                    sites.add(parsed.position)
                }
                VcfParse.Incomplete -> {
                    check(!eof) { "VCF parser asked for more bytes at EOF" }
                    val pending = end - start
                    if (start > 0) {
                        window.copyInto(window, 0, start, end)
                        start = 0
                        end = pending
                    }
                    if (end == window.size) window = window.copyOf(window.size * 2)
                    val read = fill(window, end, window.size - end)
                    if (read <= 0) eof = true else end += read
                }
            }
        }
        return sites.finish()
    }

    /**
     * Reads the header and leaves [source] positioned at the first record.
     *
     * The record parser takes this header and a window of bytes, so it can start at any later record.
     */
    internal fun readHeader(source: Source): VcfHeader = readVcfHeader(source)

    /** Bytes pulled from the source at a time. A record longer than this is carried across fills. */
    private const val RECORD_BUFFER_BYTES = 64 * 1024
}

/**
 * Sites accumulated by [Vcf.read] and [IndexedVcf.read], packed and sorted the same way.
 *
 * [calls] is filled in place by [VcfRecordParser]. [add] copies it before the next record overwrites it.
 */
internal class VcfSites(private val taxa: TaxaList, private val progress: ProgressReporter) {
    val calls = ByteArray(taxa.size)
    private val positions = ArrayList<Position>()
    private var packed = ByteArray(0)
    private var packedUsed = 0

    suspend fun add(position: Position) {
        val sites = positions.size + 1
        require(sites.toLong() * taxa.size <= Int.MAX_VALUE) {
            "VCF file is too large to hold in memory: $sites sites x ${taxa.size} taxa"
        }
        val needed = packedUsed + calls.size
        if (needed > packed.size) {
            var grown = maxOf(packed.size, calls.size).coerceAtLeast(1)
            while (grown < needed) {
                grown = if (grown > Int.MAX_VALUE / 2) needed else grown * 2
            }
            packed = packed.copyOf(grown)
        }
        calls.copyInto(packed, packedUsed)
        packedUsed = needed
        positions.add(position)
        if (positions.size % PROGRESS_SITES == 0) progress.report(positions.size, 0, "Reading VCF")
    }

    suspend fun finish(): GenotypeTable {
        val sorted = (1 until positions.size).all { positions[it - 1] <= positions[it] }
        val exact = if (packedUsed == packed.size) packed else packed.copyOf(packedUsed)
        val combined: ByteArray
        val ordered: List<Position>
        if (sorted) {
            combined = exact
            ordered = positions
        } else {
            val order = positions.indices.sortedBy { positions[it] }
            combined = ByteArray(packedUsed)
            val width = taxa.size
            order.forEachIndexed { site, from ->
                exact.copyInto(combined, site * width, from * width, (from + 1) * width)
            }
            ordered = order.map { positions[it] }
        }
        progress.report(positions.size, positions.size, "Reading VCF")
        return InMemoryGenotypeTable.wrap(taxa, PositionList(ordered), combined)
    }

    private companion object {
        /** Sites between progress reports, and so between suspension points; see `HapMap.PROGRESS_SITES`. */
        const val PROGRESS_SITES = 256
    }
}

/** A bgzipped VCF and the warnings found while inflating it. */
class VcfRead(val table: GenotypeTable, val warnings: List<String>)
