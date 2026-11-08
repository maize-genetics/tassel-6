package net.maizegenetics.tassel.io.deflate

/**
 * BGZF framing shared by the sequential gzip source and the seekable reader (SAMv1 section 4.1).
 *
 * A member with a `BC` extra subfield is one BGZF block of `BSIZE + 1` bytes. The same walk finds that
 * size, the raw DEFLATE slice, and the trailer, so neither caller keeps a second header parser.
 */
internal const val GZIP_FIXED_HEADER = 10
internal const val GZIP_TRAILER = 8
internal const val GZIP_ID1 = 0x1F
internal const val GZIP_ID2 = 0x8B
internal const val GZIP_METHOD_DEFLATE = 8
internal const val GZIP_FHCRC = 2
internal const val GZIP_FEXTRA = 4
internal const val GZIP_FNAME = 8
internal const val GZIP_FCOMMENT = 16
internal const val GZIP_RESERVED = 0xE0
internal const val BGZF_MAX_BLOCK = 65536
internal const val BGZF_MAX_UNCOMPRESSED = 65536

private const val BGZF_ID_B = 'B'.code
private const val BGZF_ID_C = 'C'.code

/**
 * Gzip member header through the extra field.
 *
 * [blockSize] is `BSIZE`, or -1 when the member has no `BC` subfield and is plain gzip.
 * [extraEnd] is the index, relative to the member start, of the first byte after the extra field.
 */
internal class GzipMemberHeader(val flags: Int, val extraEnd: Int, val blockSize: Int)

/**
 * Deflate slice and trailer of one BGZF block.
 *
 * [deflateFrom] and [deflateTo] are indexes into the byte array that was parsed.
 */
internal class BgzfPayload(val deflateFrom: Int, val deflateTo: Int, val expectedCrc: Int, val uncompressedSize: Int)

/**
 * How many bytes the gzip header occupies through the extra field.
 *
 * A negative return is the count still required before that length is known. A positive return that is
 * larger than [BGZF_MAX_BLOCK] is a header the parser will reject.
 */
internal fun gzipHeaderBytesRequired(bytes: ByteArray, offset: Int, available: Int): Int {
    if (available < GZIP_FIXED_HEADER) return -GZIP_FIXED_HEADER
    val flags = u8(bytes, offset + 3)
    if (flags and GZIP_FEXTRA == 0) return GZIP_FIXED_HEADER
    if (available < GZIP_FIXED_HEADER + 2) return -(GZIP_FIXED_HEADER + 2)
    val extraLength = u16(bytes, offset + GZIP_FIXED_HEADER)
    val total = GZIP_FIXED_HEADER + 2 + extraLength
    if (total > BGZF_MAX_BLOCK || available < total) return if (available < total) -total else total
    return total
}

/** Checks the ten-byte gzip header. The extra field, when present, is not required yet. */
internal fun checkGzipFixedHeader(bytes: ByteArray, offset: Int, available: Int) {
    if (available < GZIP_FIXED_HEADER) truncatedBgzf("gzip header")
    if (u8(bytes, offset) != GZIP_ID1 || u8(bytes, offset + 1) != GZIP_ID2) truncatedBgzf("gzip magic")
    if (u8(bytes, offset + 2) != GZIP_METHOD_DEFLATE) {
        throw IllegalArgumentException("gzip compression method is not DEFLATE")
    }
    if (u8(bytes, offset + 3) and GZIP_RESERVED != 0) {
        throw IllegalArgumentException("gzip header has reserved flags")
    }
}

/**
 * Parses the fixed header and the extra field.
 *
 * [available] bytes at [offset] must cover the header through the extra field. Name, comment, and the
 * header CRC are not consumed here; [parseBgzfPayload] skips them once `BSIZE` is known.
 */
internal fun parseGzipMemberHeader(bytes: ByteArray, offset: Int, available: Int): GzipMemberHeader {
    checkGzipFixedHeader(bytes, offset, available)
    val flags = u8(bytes, offset + 3)
    var cursor = GZIP_FIXED_HEADER
    var blockSize = -1
    if (flags and GZIP_FEXTRA != 0) {
        if (available < cursor + 2) truncatedBgzf("gzip extra length")
        val extraLength = u16(bytes, offset + cursor)
        cursor += 2
        if (extraLength > BGZF_MAX_BLOCK - cursor) truncatedBgzf("gzip extra field")
        if (available < cursor + extraLength) truncatedBgzf("gzip extra field")
        val extraEnd = cursor + extraLength
        var sub = cursor
        while (sub + 4 <= extraEnd) {
            val id1 = u8(bytes, offset + sub)
            val id2 = u8(bytes, offset + sub + 1)
            val subLength = u16(bytes, offset + sub + 2)
            val next = sub + 4 + subLength
            if (next > extraEnd) truncatedBgzf("gzip extra subfield")
            if (id1 == BGZF_ID_B && id2 == BGZF_ID_C && subLength == 2) {
                blockSize = u16(bytes, offset + sub + 4)
            }
            sub = next
        }
        cursor = extraEnd
    }
    return GzipMemberHeader(flags, cursor, blockSize)
}

/**
 * Reads the deflate range and the trailer of a BGZF block whose `BSIZE` is already known.
 *
 * The `BSIZE + 1` bytes of the member must sit at [offset].
 */
internal fun parseBgzfPayload(bytes: ByteArray, offset: Int, header: GzipMemberHeader): BgzfPayload {
    val total = header.blockSize + 1
    if (header.blockSize < 0 || total > BGZF_MAX_BLOCK || total < header.extraEnd + GZIP_TRAILER) {
        truncatedBgzf("BGZF block size")
    }
    val trailerAt = total - GZIP_TRAILER
    val deflateAt = findDeflateStart(bytes, offset, header.flags, header.extraEnd, trailerAt)
    if (deflateAt > trailerAt) truncatedBgzf("BGZF block")
    val expectedCrc = u32(bytes, offset + trailerAt)
    val expectedSize = u32(bytes, offset + trailerAt + 4)
    if (expectedSize !in 0..BGZF_MAX_UNCOMPRESSED) truncatedBgzf("BGZF uncompressed size")
    return BgzfPayload(offset + deflateAt, offset + trailerAt, expectedCrc, expectedSize)
}

/**
 * Inflates one BGZF block and checks its CRC32 and ISIZE.
 *
 * [output] receives the inflated bytes at index 0. The same [decoder] and [crc] are reused across blocks.
 */
internal fun inflateBgzf(
    decoder: InflateDecoder,
    crc: Crc32,
    input: ByteArray,
    payload: BgzfPayload,
    output: ByteArray,
): Int {
    val written = try {
        decoder.inflate(input, payload.deflateFrom, payload.deflateTo, output, 0, payload.uncompressedSize)
    } catch (e: IllegalArgumentException) {
        throw IllegalArgumentException("truncated BGZF block (${e.message})")
    }
    crc.reset()
    crc.update(output, 0, written)
    if (crc.value != payload.expectedCrc) {
        throw IllegalArgumentException("gzip CRC32 does not match the inflated bytes")
    }
    if (written != payload.uncompressedSize) {
        throw IllegalArgumentException("gzip ISIZE does not match the inflated bytes")
    }
    return written
}

internal fun truncatedBgzf(where: String): Nothing = throw IllegalArgumentException("truncated BGZF block ($where)")

private fun findDeflateStart(bytes: ByteArray, offset: Int, flags: Int, from: Int, limit: Int): Int {
    var cursor = from
    if (flags and GZIP_FNAME != 0) cursor = skipCString(bytes, offset, cursor, limit)
    if (flags and GZIP_FCOMMENT != 0) cursor = skipCString(bytes, offset, cursor, limit)
    if (flags and GZIP_FHCRC != 0) {
        if (cursor + 2 > limit) truncatedBgzf("gzip header CRC")
        cursor += 2
    }
    return cursor
}

private fun skipCString(bytes: ByteArray, offset: Int, from: Int, limit: Int): Int {
    var cursor = from
    while (cursor < limit) {
        if (bytes[offset + cursor].toInt() == 0) return cursor + 1
        cursor++
    }
    truncatedBgzf("gzip header text")
}

private fun u8(bytes: ByteArray, at: Int): Int = bytes[at].toInt() and 0xFF

private fun u16(bytes: ByteArray, at: Int): Int = u8(bytes, at) or (u8(bytes, at + 1) shl 8)

private fun u32(bytes: ByteArray, at: Int): Int = u16(bytes, at) or (u16(bytes, at + 2) shl 16)
