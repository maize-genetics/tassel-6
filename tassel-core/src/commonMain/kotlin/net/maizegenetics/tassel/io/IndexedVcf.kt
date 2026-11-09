package net.maizegenetics.tassel.io

import kotlinx.io.Buffer
import net.maizegenetics.tassel.io.index.GenomicIndex
import net.maizegenetics.tassel.io.index.readGenomicIndex
import net.maizegenetics.tassel.io.index.unsignedLess
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.taxa.TaxaList

/**
 * A bgzipped VCF opened with its CSI or TBI index.
 *
 * [open] reads the header once, from offset 0, and reads the index into memory. The index is told
 * from CSI or TBI by its magic bytes ([readGenomicIndex]). [contigs] are the index names as written,
 * with the length from a `##contig` line when the header id matches that name. [read] returns the
 * sites whose reference span overlaps one region. The table is built the same way as [Vcf.read]:
 * the same site order, allele rules, and errors.
 *
 * [close] closes the BGZF reader only. The caller owns both sources and keeps the data source open
 * for every [read].
 *
 * @see net.maizegenetics.dna.snp.io.BuilderFromVCF
 */
class IndexedVcf internal constructor(
    private val reader: BgzfReader,
    private val index: GenomicIndex,
    private val header: VcfHeader,
    /** Index contig names, in index order, with `##contig` lengths where the header has that id. */
    val contigs: List<Contig>,
    /**
     * `##contig` ids, in header order. Empty when the header lists none, which is the tutorial VCF.
     * An index may name a subset of these: wheat declares `chr6D` and then leaves it out of the CSI.
     */
    val headerContigIds: List<String>,
) : AutoCloseable {
    /**
     * One contig the index names.
     *
     * [length] comes from `##contig` when the header has this [name]. It is null when the header
     * does not, which is the case for the tutorial VCF.
     */
    data class Contig(val name: String, val length: Long?)

    /**
     * The sites in [region].
     *
     * An unknown contig and an inverted region are errors. A known contig with nothing in range, or a
     * region past the index's coordinate limit, is a table with the file's taxa and no sites. A `POS`
     * past [Int.MAX_VALUE] fails when the site becomes a position; the index query and the BGZF seek
     * do not.
     *
     * @throws IllegalArgumentException for the cases above, and for the record errors [Vcf.read] raises.
     */
    suspend fun read(region: VcfRegion, progress: ProgressReporter = ProgressReporter.None): VcfRead {
        if (region.start != null && region.end != null && region.start > region.end) {
            throw IllegalArgumentException(
                "VCF region ${region.contig}:${region.start}-${region.end} is inverted",
            )
        }
        val names = index.contigs.map { it.name }
        val resolved = resolveContigName(region.contig, names)
        val contig = index.contigIndex(resolved)
        val sites = VcfSites(TaxaList.of(header.taxaNames), progress)
        val beg = (region.start ?: 1L) - 1L
        if (beg >= index.maximumCoordinate) return result(sites)
        val endExclusive = region.end ?: index.maximumCoordinate
        val chunks = index.chunks(contig, beg, endExclusive)
        if (chunks.isEmpty()) return result(sites)

        val parser = VcfRecordParser(header)
        for (chunk in chunks) {
            if (!readChunk(chunk.end, chunk.start, resolved, region, parser, sites)) break
        }
        return result(sites)
    }

    private suspend fun readChunk(
        chunkEnd: Long,
        chunkStart: Long,
        resolved: String,
        region: VcfRegion,
        parser: VcfRecordParser,
        sites: VcfSites,
    ): Boolean {
        reader.chunkLimit(chunkEnd)
        reader.seek(chunkStart)
        while (unsignedLess(reader.tell(), chunkEnd)) {
            val line = reader.readLine() ?: return true
            if (line.isEmpty() || line[0] == '#') continue
            val tab = line.indexOf('\t')
            if (tab <= 0) {
                throw IllegalArgumentException("VCF record is missing a contig: ${vcfExcerpt(line)}")
            }
            if (line.substring(0, tab) != resolved) continue
            val span = vcfReferenceSpan(line)
            if (region.end != null && span.pos > region.end) return false
            if (!span.overlaps(region)) continue
            val bytes = line.encodeToByteArray()
            when (val parsed = parser.tryParse(bytes, 0, bytes.size, true, sites.calls)) {
                is VcfParse.Site -> sites.add(parsed.position)
                is VcfParse.Skip -> Unit
                else -> throw IllegalArgumentException("VCF record did not parse: ${vcfExcerpt(line)}")
            }
        }
        return true
    }

    private suspend fun result(sites: VcfSites): VcfRead = VcfRead(sites.finish(), reader.warnings.toList())

    override fun close() {
        reader.close()
    }

    companion object {
        /**
         * Opens [data] with [index].
         *
         * Neither source is closed. [index] is read fully here; [data] stays open for [read].
         */
        suspend fun open(data: RandomAccessSource, index: RandomAccessSource): IndexedVcf {
            val genomic = readGenomicIndex(Buffer().also { it.write(readAll(index)) })
            val reader = BgzfReader(data)
            try {
                reader.exactBlocks(true)
                val header = readVcfHeaderLines { reader.readLine() }
                reader.exactBlocks(false)
                return IndexedVcf(
                    reader,
                    genomic,
                    header,
                    contigsOf(header, genomic),
                    header.contigs.map { it.id },
                )
            } catch (e: Throwable) {
                reader.close()
                throw e
            }
        }
    }
}

private fun contigsOf(header: VcfHeader, index: GenomicIndex): List<IndexedVcf.Contig> {
    val lengths = HashMap<String, Long?>(header.contigs.size)
    for (contig in header.contigs) {
        if (contig.id !in lengths) lengths[contig.id] = contig.length
    }
    return index.contigs.map { IndexedVcf.Contig(it.name, lengths[it.name]) }
}

private suspend fun readAll(source: RandomAccessSource): ByteArray {
    val size = source.size
    if (size > Int.MAX_VALUE) {
        throw IllegalArgumentException("genomic index is larger than ${Int.MAX_VALUE} bytes")
    }
    if (size == 0L) throw IllegalArgumentException("truncated genomic index")
    val bytes = ByteArray(size.toInt())
    var at = 0
    while (at < bytes.size) {
        val count = source.readAt(at.toLong(), bytes, at, bytes.size - at)
        if (count <= 0) throw IllegalArgumentException("truncated genomic index")
        at += count
    }
    return bytes
}
