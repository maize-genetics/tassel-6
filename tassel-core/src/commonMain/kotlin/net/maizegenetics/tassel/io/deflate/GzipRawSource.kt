package net.maizegenetics.tassel.io.deflate

import kotlinx.io.Buffer
import kotlinx.io.RawSource
import net.maizegenetics.tassel.io.MISSING_BGZF_EOF

/**
 * Sequential gzip reader (RFC 1952) over BGZF members (SAMv1 section 4.1) and plain gzip members.
 *
 * A member that carries the `BC` extra subfield is read whole (`BSIZE + 1` bytes) and inflated in one
 * shot. A member without it is plain gzip: the deflate stream is resumed as input arrives, and only the
 * last 32 KiB of output is kept for back-references. CRC32 and ISIZE are checked for every member. ISIZE
 * is the low 32 bits of the uncompressed length. An empty member, including the 28-byte EOF block a `cat`
 * of two bgzipped files leaves in the middle, is verified and then skipped. [readAtMostTo] returns no
 * bytes for that member: a `0` result means EOF to [net.maizegenetics.tassel.io.Vcf.read].
 *
 * Plain and BGZF members may alternate. [MISSING_BGZF_EOF] is recorded only when every member was BGZF
 * and the last one was not the empty EOF block. A plain gzip file does not warn.
 *
 * When `FHCRC` is set on a plain member, the header CRC16 is checked. BGZF members still skip that field
 * without checking it: a synthetic BGZF block uses a zero CRC16, and the one-shot path is unchanged.
 *
 * [close] does not close [upstream]. The caller owns that source and must keep it open while reading.
 *
 * Not safe for concurrent calls.
 *
 * @see net.maizegenetics.util.Utils.getBufferedReader
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
    private val ahead = ByteArray(AHEAD_FILL + RETURNED_ROOM)
    private val one = ByteArray(1)

    private var aheadPos = 0
    private var aheadEnd = 0
    private var outputPos = 0
    private var outputEnd = 0
    private var sawMember = false
    private var sawBgzf = false
    private var sawPlain = false
    private var lastWasEmpty = false
    private var plainActive = false
    private var plainSize = 0L
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

    /** Inflates the next non-empty member, or the next slice of a plain member. Empty members stay inside this loop. */
    private fun fill(): Boolean {
        if (finished) return false
        while (true) {
            if (plainActive) {
                if (producePlain()) return true
                if (finished) return false
                continue
            }
            when (readMember()) {
                null -> {
                    if (sawBgzf && !sawPlain && !lastWasEmpty) recordedWarnings.add(MISSING_BGZF_EOF)
                    finished = true
                    return false
                }
                0 -> Unit
                else -> return true
            }
        }
    }

    /**
     * Reads one gzip member, or starts a plain member.
     *
     * @return uncompressed size, `0` for an empty member or a plain member that has not produced bytes yet,
     *   or `null` at a clean end of input.
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
        if (flags and RESERVED != 0) {
            throw IllegalArgumentException("gzip header has reserved flags")
        }
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
        if (blockSize < 0) {
            beginPlain(flags, cursor)
            return 0
        }
        sawBgzf = true
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

    /**
     * Header bytes through [headerEnd] are already in [member]. Name, comment, and the optional header CRC
     * follow in the byte stream. The deflate payload is left for [producePlain].
     */
    private fun beginPlain(flags: Int, headerEnd: Int) {
        crc.reset()
        crc.update(member, 0, headerEnd)
        if (flags and FNAME != 0) skipPlainText()
        if (flags and FCOMMENT != 0) skipPlainText()
        if (flags and FHCRC != 0) {
            val stored = readU16Plain()
            if (stored != (crc.value and 0xFFFF)) {
                throw IllegalArgumentException("gzip header CRC does not match")
            }
        }
        crc.reset()
        decoder.resetStream()
        plainSize = 0
        plainActive = true
        sawPlain = true
        sawMember = true
    }

    /** Inflates the next slice of the open plain member. Empty finished members return false so [fill] can continue. */
    private fun producePlain(): Boolean {
        while (true) {
            val noMore = !ensureAhead()
            val step = decoder.inflateSome(ahead, aheadPos, aheadEnd, inflated, 0, inflated.size, noMore)
            aheadPos += step.consumed
            if (step.returned > 0) acceptReturned(step.returned)
            if (step.written > 0) {
                crc.update(inflated, 0, step.written)
                plainSize += step.written
                outputPos = 0
                outputEnd = step.written
            }
            if (step.finished) finishPlainTrailer()
            if (step.written > 0) return true
            if (!plainActive) return false
            if (step.consumed == 0) {
                if (noMore) plainTruncated("deflate")
                throw IllegalStateException("inflate made no progress")
            }
        }
    }

    private fun finishPlainTrailer() {
        if (readExact(0, TRAILER) < TRAILER) plainTruncated("trailer")
        val expectedCrc = u32(0)
        val expectedSize = u32(4)
        if (crc.value != expectedCrc) {
            throw IllegalArgumentException("gzip CRC32 does not match the inflated bytes")
        }
        if (plainSize.toInt() != expectedSize) {
            throw IllegalArgumentException("gzip ISIZE does not match the inflated bytes")
        }
        plainActive = false
        lastWasEmpty = plainSize == 0L
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

    private fun skipPlainText() {
        while (true) {
            val value = readStreamByte()
            if (value < 0) plainTruncated("header text")
            one[0] = value.toByte()
            crc.update(one, 0, 1)
            if (value == 0) return
        }
    }

    private fun readU16Plain(): Int {
        val lo = readStreamByte()
        val hi = readStreamByte()
        if (lo < 0 || hi < 0) plainTruncated("header CRC")
        return lo or (hi shl 8)
    }

    /** Copies [count] bytes into [member] at [offset]. Returns how many were available. */
    private fun readExact(offset: Int, count: Int): Int {
        var filled = 0
        while (filled < count) {
            if (aheadPos < aheadEnd) {
                val n = minOf(count - filled, aheadEnd - aheadPos)
                ahead.copyInto(member, offset + filled, aheadPos, aheadPos + n)
                aheadPos += n
                filled += n
                continue
            }
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

    private fun readStreamByte(): Int {
        if (aheadPos < aheadEnd) return ahead[aheadPos++].toInt() and 0xFF
        if (pending.size == 0L) {
            val read = upstream.readAtMostTo(pending, 1)
            if (read == -1L) return -1
            if (read == 0L) throw IllegalStateException("upstream read returned no bytes")
        }
        val n = pending.readAtMostTo(one, 0, 1)
        if (n <= 0) return -1
        return one[0].toInt() and 0xFF
    }

    private fun ensureAhead(): Boolean {
        if (aheadPos < aheadEnd) return true
        val n = readUpstream(ahead, 0, AHEAD_FILL)
        if (n <= 0) {
            aheadPos = 0
            aheadEnd = 0
            return false
        }
        aheadPos = 0
        aheadEnd = n
        return true
    }

    private fun readUpstream(dest: ByteArray, offset: Int, count: Int): Int {
        if (count == 0) return 0
        if (pending.size == 0L) {
            val read = upstream.readAtMostTo(pending, count.toLong())
            if (read == -1L) return -1
            if (read == 0L) throw IllegalStateException("upstream read returned no bytes")
        }
        val n = pending.readAtMostTo(dest, offset, offset + count)
        return if (n <= 0) -1 else n
    }

    /** [count] bytes from [InflateDecoder.returnedBytes] are older than the unread tail of [ahead]. */
    private fun acceptReturned(count: Int) {
        val remain = aheadEnd - aheadPos
        if (count + remain > ahead.size) throw IllegalStateException("gzip readahead overflow")
        if (aheadPos >= count) {
            aheadPos -= count
            decoder.returnedBytes.copyInto(ahead, aheadPos, 0, count)
            return
        }
        var index = remain - 1
        while (index >= 0) {
            ahead[count + index] = ahead[aheadPos + index]
            index--
        }
        decoder.returnedBytes.copyInto(ahead, 0, 0, count)
        aheadPos = 0
        aheadEnd = count + remain
    }

    private fun u8(at: Int): Int = member[at].toInt() and 0xFF

    private fun u16(at: Int): Int = u8(at) or (u8(at + 1) shl 8)

    private fun u32(at: Int): Int = u16(at) or (u16(at + 2) shl 16)

    private fun truncated(where: String): Nothing = throw IllegalArgumentException("truncated BGZF block ($where)")

    private fun plainTruncated(where: String): Nothing =
        throw IllegalArgumentException("truncated gzip member ($where)")

    private companion object {
        const val MAX_BLOCK = 65536
        const val MAX_UNCOMPRESSED = 65536
        const val AHEAD_FILL = 16 * 1024
        const val RETURNED_ROOM = 16
        const val FIXED_HEADER = 10
        const val TRAILER = 8
        const val GZIP_ID1 = 0x1F
        const val GZIP_ID2 = 0x8B
        const val METHOD_DEFLATE = 8
        const val FHCRC = 2
        const val FEXTRA = 4
        const val FNAME = 8
        const val FCOMMENT = 16
        const val RESERVED = 0xE0
        const val ID_B = 'B'.code
        const val ID_C = 'C'.code
    }
}
