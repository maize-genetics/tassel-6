package net.maizegenetics.tassel.io.deflate

import kotlinx.io.Buffer
import kotlinx.io.RawSource
import net.maizegenetics.tassel.io.MISSING_BGZF_EOF

/**
 * Sequential BGZF reader (SAMv1 section 4.1) over a gzip member stream (RFC 1952).
 *
 * Each member that carries the `BC` extra subfield is read whole (`BSIZE + 1` bytes) and inflated
 * in one shot. CRC32 and ISIZE are checked for every member. An empty member, including the 28-byte
 * EOF block a `cat` of two bgzipped files leaves in the middle, is verified and then skipped.
 * [readAtMostTo] returns no bytes for that member: a `0` result means EOF to [net.maizegenetics.tassel.io.Vcf.read].
 *
 * A member with no `BC` subfield fails. Plain gzip is tassel-6#85. When the upstream ends on a member
 * boundary and the last member was not empty, [warnings] gains [MISSING_BGZF_EOF]. A short member throws.
 *
 * [close] does not close [upstream]. The caller owns that source and must keep it open while reading.
 *
 * Not safe for concurrent calls.
 */
internal class GzipRawSource(private val upstream: RawSource) : RawSource {
    private val recordedWarnings = ArrayList<String>(1)

    /** Filled when the inflated stream is exhausted. */
    val warnings: List<String> get() = recordedWarnings

    private val decoder = InflateDecoder()
    private val crc = Crc32()
    private val pending = Buffer()
    private val member = ByteArray(MAX_BLOCK)
    private val inflated = ByteArray(MAX_UNCOMPRESSED)

    private var outputPos = 0
    private var outputEnd = 0
    private var sawMember = false
    private var lastWasEmpty = false
    private var finished = false
    private var closed = false

    override fun readAtMostTo(sink: Buffer, byteCount: Long): Long {
        check(!closed) { "Source is closed." }
        require(byteCount >= 0L) { "byteCount: $byteCount" }
        if (byteCount == 0L) return 0L
        if (finished) return -1L

        var copied = 0L
        while (copied < byteCount) {
            if (outputPos == outputEnd && !fill()) {
                return if (copied == 0L) -1L else copied
            }
            val count = minOf(byteCount - copied, (outputEnd - outputPos).toLong()).toInt()
            sink.write(inflated, outputPos, outputPos + count)
            outputPos += count
            copied += count
        }
        return copied
    }

    override fun close() {
        closed = true
    }

    /** Inflates the next non-empty member. Empty members stay inside this loop. */
    private fun fill(): Boolean {
        if (finished) return false
        while (true) {
            when (readMember()) {
                null -> {
                    if (sawMember && !lastWasEmpty) recordedWarnings.add(MISSING_BGZF_EOF)
                    finished = true
                    return false
                }
                0 -> Unit
                else -> return true
            }
        }
    }

    /**
     * Reads one gzip member.
     *
     * @return uncompressed size, `0` for an empty member, or `null` at a clean end of input.
     */
    private fun readMember(): Int? {
        val header = readExact(0, FIXED_HEADER)
        if (header == 0) return null
        if (header < FIXED_HEADER) truncated("gzip header")
        if (u8(0) != GZIP_ID1 || u8(1) != GZIP_ID2) truncated("gzip magic")
        if (u8(2) != METHOD_DEFLATE) {
            throw IllegalArgumentException("gzip compression method is not DEFLATE")
        }
        val flags = u8(3)
        var cursor = FIXED_HEADER
        var blockSize = -1
        if (flags and FEXTRA != 0) {
            if (readExact(cursor, 2) < 2) truncated("gzip extra length")
            val extraLength = u16(cursor)
            cursor += 2
            if (extraLength > MAX_BLOCK - cursor) truncated("gzip extra field")
            if (readExact(cursor, extraLength) < extraLength) truncated("gzip extra field")
            val extraEnd = cursor + extraLength
            var sub = cursor
            while (sub + 4 <= extraEnd) {
                val id1 = u8(sub)
                val id2 = u8(sub + 1)
                val subLength = u16(sub + 2)
                val next = sub + 4 + subLength
                if (next > extraEnd) truncated("gzip extra subfield")
                if (id1 == ID_B && id2 == ID_C && subLength == 2) blockSize = u16(sub + 4)
                sub = next
            }
            cursor = extraEnd
        }
        if (blockSize < 0) plainGzip()
        val total = blockSize + 1
        if (total > MAX_BLOCK || total < cursor + TRAILER) truncated("BGZF block size")
        val rest = total - cursor
        if (readExact(cursor, rest) < rest) truncated("BGZF block")

        val deflateEnd = total - TRAILER
        val deflateAt = findDeflateStart(flags, cursor, deflateEnd)
        if (deflateAt > deflateEnd) truncated("BGZF block")
        val expectedCrc = u32(deflateEnd)
        val expectedSize = u32(deflateEnd + 4)
        if (expectedSize !in 0..MAX_UNCOMPRESSED) truncated("BGZF uncompressed size")

        val written = try {
            decoder.inflate(member, deflateAt, deflateEnd, inflated, 0, expectedSize)
        } catch (e: IllegalArgumentException) {
            throw IllegalArgumentException("truncated BGZF block (${e.message})")
        }
        crc.reset()
        crc.update(inflated, 0, written)
        if (crc.value != expectedCrc) {
            throw IllegalArgumentException("gzip CRC32 does not match the inflated bytes")
        }
        if (written != expectedSize) {
            throw IllegalArgumentException("gzip ISIZE does not match the inflated bytes")
        }
        sawMember = true
        lastWasEmpty = written == 0
        if (written > 0) {
            outputPos = 0
            outputEnd = written
        }
        return written
    }

    private fun findDeflateStart(flags: Int, from: Int, limit: Int): Int {
        var cursor = from
        if (flags and FNAME != 0) cursor = skipCString(cursor, limit)
        if (flags and FCOMMENT != 0) cursor = skipCString(cursor, limit)
        if (flags and FHCRC != 0) {
            if (cursor + 2 > limit) truncated("gzip header CRC")
            cursor += 2
        }
        return cursor
    }

    private fun skipCString(from: Int, limit: Int): Int {
        var cursor = from
        while (cursor < limit) {
            if (member[cursor].toInt() == 0) return cursor + 1
            cursor++
        }
        truncated("gzip header text")
    }

    /** Copies [count] bytes into [member] at [offset]. Returns how many were available. */
    private fun readExact(offset: Int, count: Int): Int {
        var filled = 0
        while (filled < count) {
            if (pending.size == 0L) {
                val read = upstream.readAtMostTo(pending, (count - filled).toLong())
                if (read == -1L) return filled
                if (read == 0L) throw IllegalStateException("upstream read returned no bytes")
            }
            val n = pending.readAtMostTo(member, offset + filled, offset + count)
            if (n <= 0) return filled
            filled += n
        }
        return filled
    }

    private fun u8(at: Int): Int = member[at].toInt() and 0xFF

    private fun u16(at: Int): Int = u8(at) or (u8(at + 1) shl 8)

    private fun u32(at: Int): Int = u16(at) or (u16(at + 2) shl 16)

    private fun truncated(where: String): Nothing = throw IllegalArgumentException("truncated BGZF block ($where)")

    private fun plainGzip(): Nothing = throw IllegalArgumentException(
        "plain gzip is not supported yet (tassel-6#85); this member has no BGZF BC subfield",
    )

    private companion object {
        const val MAX_BLOCK = 65536
        const val MAX_UNCOMPRESSED = 65536
        const val FIXED_HEADER = 10
        const val TRAILER = 8
        const val GZIP_ID1 = 0x1F
        const val GZIP_ID2 = 0x8B
        const val METHOD_DEFLATE = 8
        const val FHCRC = 2
        const val FEXTRA = 4
        const val FNAME = 8
        const val FCOMMENT = 16
        const val ID_B = 'B'.code
        const val ID_C = 'C'.code
    }
}
