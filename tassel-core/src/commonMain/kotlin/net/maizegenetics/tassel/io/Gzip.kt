package net.maizegenetics.tassel.io

import kotlinx.io.Source
import kotlinx.io.buffered
import net.maizegenetics.tassel.io.deflate.GzipRawSource

/**
 * Recorded on [Decompressed.warnings] when a BGZF stream ends on a non-empty block.
 *
 * htslib warns with this sentence and keeps reading. The read succeeds (ADR-0021).
 */
const val MISSING_BGZF_EOF: String = "EOF marker is absent. The input may be truncated"

/**
 * A [Source] for [Vcf.read] or [HapMap.read], plus warnings discovered while inflating it.
 *
 * [warnings] stays empty until [source] is exhausted. For bytes that are not gzip, [source] is the
 * source that was peeked and [warnings] stays empty.
 */
class Decompressed internal constructor(val source: Source, private val gzip: GzipRawSource? = null) {
    val warnings: List<String>
        get() = gzip?.warnings ?: emptyList()
}

/**
 * Returns [this] when the first two bytes are not the gzip magic (`1f 8b`), and a BGZF reader when
 * they are.
 *
 * Detection is the bytes, so a caller never decides compression from a file name. A BGZF member is
 * inflated one block at a time. A plain gzip member, including one longer than a BGZF block, is inflated
 * as bytes arrive, and plain and BGZF members may alternate. CRC32 and ISIZE are checked for every member.
 * A missing terminal BGZF EOF block is recorded on [Decompressed.warnings] once the inflated stream is
 * exhausted, and only when every member was BGZF.
 *
 * The returned source does not close [this]. Keep [this] open while reading, and close it afterwards.
 *
 * @see net.maizegenetics.util.Utils.getBufferedReader
 */
fun Source.decompressing(): Decompressed {
    val gzip = peek().use { preview ->
        preview.request(2) &&
            (preview.readByte().toInt() and 0xFF) == GZIP_ID1 &&
            (preview.readByte().toInt() and 0xFF) == GZIP_ID2
    }
    if (!gzip) return Decompressed(this)
    val raw = GzipRawSource(this)
    return Decompressed(raw.buffered(), raw)
}

private const val GZIP_ID1 = 0x1F
private const val GZIP_ID2 = 0x8B
