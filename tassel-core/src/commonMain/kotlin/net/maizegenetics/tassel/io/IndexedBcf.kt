package net.maizegenetics.tassel.io

import kotlinx.io.Buffer
import net.maizegenetics.tassel.io.index.GenomicIndex
import net.maizegenetics.tassel.io.index.readGenomicIndex
import net.maizegenetics.tassel.io.index.unsignedLess
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.taxa.TaxaList

/**
 * A bgzipped BCF opened with its CSI index.
 *
 * [open] reads the header once, from offset 0, and reads the index into memory. Contig names come from the header
 * dictionary, because a BCF CSI stores none (`l_aux` is 0). [contigs] is that list, in header order, with a length
 * only where the `##contig` line has one. [read] returns the sites whose reference span overlaps one region.
 * The table is built the same way as [Bcf.read].
 *
 * [close] closes the BGZF reader only. The caller owns both sources and keeps the data source open for every [read].
 * A TBI, or a CSI that already carries tabix names, is rejected by [readGenomicIndex]. Raw BCF has no index; this
 * reader is BGZF only.
 *
 * @see IndexedVcf
 * @see Bcf.read
 */
internal class IndexedBcf internal constructor(
    private val reader: BgzfReader,
    private val index: GenomicIndex,
    private val header: BcfHeader,
    /** Header contig names, in dictionary order, with `##contig` lengths where the header has them. */
    val contigs: List<IndexedGenotypes.Contig>,
) : AutoCloseable {
    /**
     * The sites in [region].
     *
     * An unknown contig and an inverted region are errors. A known contig with nothing in range, or a region past
     * the index's coordinate limit, is a table with the file's taxa and no sites. A `POS` past [Int.MAX_VALUE]
     * fails when the site becomes a position.
     *
     * @throws IllegalArgumentException for the cases above, and for the record errors [Bcf.read] raises.
     */
    suspend fun read(region: VcfRegion, progress: ProgressReporter = ProgressReporter.None): IndexedGenotypeRead {
        if (region.start != null && region.end != null && region.start > region.end) {
            throw IllegalArgumentException(
                "VCF region ${region.contig}:${region.start}-${region.end} is inverted",
            )
        }
        val names = contigs.map { it.name }
        val resolved = resolveContigName(region.contig, names)
        val contig = index.contigIndex(resolved)
        val sites = VcfSites(TaxaList.of(header.vcf.taxaNames), progress, "Reading BCF")
        val beg = (region.start ?: 1L) - 1L
        if (beg >= index.maximumCoordinate) return result(sites)
        val endExclusive = region.end ?: index.maximumCoordinate
        val chunks = index.chunks(contig, beg, endExclusive)
        if (chunks.isEmpty()) return result(sites)

        for (chunk in chunks) {
            if (!readChunk(chunk.end, chunk.start, resolved, region, sites)) break
        }
        return result(sites)
    }

    private suspend fun readChunk(
        chunkEnd: Long,
        chunkStart: Long,
        resolved: String,
        region: VcfRegion,
        sites: VcfSites,
    ): Boolean {
        reader.chunkLimit(chunkEnd)
        reader.seek(chunkStart)
        while (unsignedLess(reader.tell(), chunkEnd)) {
            val record = readRecord(reader) ?: return true
            if (bcfContigName(record.body, header) != resolved) continue
            val span = bcfReferenceSpan(record.body, record.sharedLength, header)
            if (region.end != null && span.pos > region.end) return false
            if (!span.overlaps(region)) continue
            sites.add(bcfPosition(record.body, record.sharedLength, header, sites.calls))
        }
        return true
    }

    private suspend fun result(sites: VcfSites): IndexedGenotypeRead =
        IndexedGenotypeRead(sites.finish(), reader.warnings.toList())

    override fun close() {
        reader.close()
    }

    companion object {
        /**
         * Opens [data] with [index].
         *
         * Neither source is closed. [index] is read fully here; [data] stays open for [read]. The header is read
         * from offset 0 with one block at a time, so a later region read is the first use of the index chunks.
         */
        suspend fun open(data: RandomAccessSource, index: RandomAccessSource): IndexedBcf {
            val indexBytes = readAll(index)
            val reader = BgzfReader(data)
            try {
                reader.exactBlocks(true)
                val headerBytes = readHeader(reader)
                val header = readBcfHeader(Buffer().also { it.write(headerBytes) })
                reader.exactBlocks(false)
                val contigs = contigsOf(header)
                val genomic = readGenomicIndex(
                    Buffer().also { it.write(indexBytes) },
                    contigs.map { it.name },
                )
                return IndexedBcf(reader, genomic, header, contigs)
            } catch (e: Throwable) {
                reader.close()
                throw e
            }
        }
    }
}

private class RecordBytes(val body: ByteArray, val sharedLength: Int)

private fun contigsOf(header: BcfHeader): List<IndexedGenotypes.Contig> {
    val contigs = ArrayList<IndexedGenotypes.Contig>(header.contigs.size)
    for ((index, contig) in header.contigs.withIndex()) {
        if (contig == null) {
            throw IllegalArgumentException("BCF contig dictionary has an empty slot at $index")
        }
        contigs.add(IndexedGenotypes.Contig(contig.id, contig.length))
    }
    return contigs
}

/** The BCF header at offset 0: magic, `l_text`, and the text, including its trailing NUL. */
private suspend fun readHeader(reader: BgzfReader): ByteArray {
    val preamble = readExact(reader, BCF_MAGIC_LENGTH + UINT32_BYTES, TRUNCATED_BCF_HEADER)
        ?: throw IllegalArgumentException(TRUNCATED_BCF_HEADER)
    val textLength = preamble.uint32Le(BCF_MAGIC_LENGTH)
    if (textLength == 0L || textLength > Int.MAX_VALUE.toLong()) {
        throw IllegalArgumentException(TRUNCATED_BCF_HEADER)
    }
    val text = readExact(reader, textLength.toInt(), TRUNCATED_BCF_HEADER)
        ?: throw IllegalArgumentException(TRUNCATED_BCF_HEADER)
    return preamble + text
}

/**
 * One record body, or null when [reader] is exhausted before the length prefix.
 *
 * The prefix is 8 bytes. The body is `l_shared + l_indiv` and may cross a BGZF block, including the chunk end:
 * the caller stops on [BgzfReader.tell] before the next record, not in the middle of this one.
 */
private suspend fun readRecord(reader: BgzfReader): RecordBytes? {
    val prefix = readExact(reader, BCF_RECORD_PREFIX_BYTES, TRUNCATED_BCF_RECORD) ?: return null
    val shared = prefix.uint32Le(0)
    val indiv = prefix.uint32Le(BCF_RECORD_PREFIX_BYTES / 2)
    val body = readExact(reader, bcfBodyLength(shared, indiv), TRUNCATED_BCF_RECORD)
        ?: throw IllegalArgumentException(TRUNCATED_BCF_RECORD)
    return RecordBytes(body, shared.toInt())
}

/** [count] inflated bytes, or null when nothing remains. A short read is [truncated]. */
private suspend fun readExact(reader: BgzfReader, count: Int, truncated: String): ByteArray? {
    if (count == 0) return ByteArray(0)
    val out = ByteArray(count)
    var at = 0
    while (at < count) {
        val n = reader.read(out, at, count - at)
        if (n <= 0) {
            if (at == 0) return null
            throw IllegalArgumentException(truncated)
        }
        at += n
    }
    return out
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
