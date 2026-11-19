package net.maizegenetics.tassel.io.deflate

import kotlinx.io.Sink

/**
 * One BGZF member (SAMv1 section 4.1) around raw DEFLATE or a stored block.
 *
 * The gzip header is htslib's `g_magic`: MTIME 0, XFL 0, OS 255, XLEN 6, and the `BC` subfield.
 * `BSIZE` is the member length minus one. CRC32 and ISIZE cover [uncompressed], including when
 * [deflate] is a compressed payload. A null payload, or one that would make the member larger than
 * [BGZF_MAX_BLOCK], is written as one final stored block (`BFINAL` 1, `BTYPE` 00) the way
 * `bgzf_compress` does at level 0. Returns the number of bytes written.
 */
internal fun Sink.writeBgzfMember(
    uncompressed: ByteArray,
    from: Int,
    to: Int,
    deflate: ByteArray?,
    crc32: Crc32,
): Int {
    if (from < 0 || to < from || to > uncompressed.size) {
        throw IllegalArgumentException("invalid BGZF member range")
    }
    val count = to - from
    if (count > 0xFFFF) throw IllegalArgumentException("stored block is longer than 65535 bytes")
    val compressed = deflate?.takeIf { it.size <= BGZF_MAX_BLOCK - MEMBER_FRAMING }
    val memberLength = if (compressed == null) {
        count + MEMBER_FRAMING + STORED_HEADER
    } else {
        compressed.size + MEMBER_FRAMING
    }
    if (memberLength > BGZF_MAX_BLOCK) {
        throw IllegalArgumentException("BGZF block is larger than $BGZF_MAX_BLOCK bytes")
    }

    val header = HEADER_TEMPLATE.copyOf()
    val blockSize = memberLength - 1
    header[16] = blockSize.toByte()
    header[17] = (blockSize ushr 8).toByte()
    write(header)

    if (compressed == null) {
        val prefix = ByteArray(STORED_HEADER)
        prefix[0] = STORED_FINAL
        prefix[1] = count.toByte()
        prefix[2] = (count ushr 8).toByte()
        val complement = count.inv() and 0xFFFF
        prefix[3] = complement.toByte()
        prefix[4] = (complement ushr 8).toByte()
        write(prefix)
        if (count > 0) write(uncompressed, from, to)
    } else {
        write(compressed)
    }

    crc32.reset()
    crc32.update(uncompressed, from, to)
    val trailer = ByteArray(GZIP_TRAILER)
    putU32(trailer, 0, crc32.value)
    putU32(trailer, 4, count)
    write(trailer)
    return memberLength
}

private const val MEMBER_FRAMING = BGZF_HEADER_LENGTH + GZIP_TRAILER
private const val STORED_HEADER = 5
private const val STORED_FINAL: Byte = 1

private val HEADER_TEMPLATE: ByteArray = byteArrayOf(
    0x1f,
    0x8b.toByte(),
    0x08,
    0x04,
    0x00,
    0x00,
    0x00,
    0x00,
    0x00,
    0xff.toByte(),
    0x06,
    0x00,
    0x42,
    0x43,
    0x02,
    0x00,
    0x00,
    0x00,
)

private fun putU32(target: ByteArray, at: Int, value: Int) {
    target[at] = value.toByte()
    target[at + 1] = (value ushr 8).toByte()
    target[at + 2] = (value ushr 16).toByte()
    target[at + 3] = (value ushr 24).toByte()
}
