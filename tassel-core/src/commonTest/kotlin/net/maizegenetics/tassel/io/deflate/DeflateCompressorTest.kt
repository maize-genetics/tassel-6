package net.maizegenetics.tassel.io.deflate

import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertNotNull
import kotlin.test.assertTrue

/** Raw DEFLATE from [DeflateCompressor], inflated by [InflateDecoder]. */
class DeflateCompressorTest {
    private val compressor = DeflateCompressor()
    private val decoder = InflateDecoder()

    @Test
    fun emptyBlockIsTheFixedEndOfBlock() {
        val payload = compressor.compress(ByteArray(0), 0, 0)
        assertContentEquals(byteArrayOf(0x03, 0x00), payload)
    }

    @Test
    fun literalsUseAFixedBlock() {
        val plain = "abc".encodeToByteArray()
        val payload = roundTrip(plain)
        assertEquals(FIXED, blockType(payload))
    }

    @Test
    fun aLongRunUsesAFixedBlock() {
        // One literal plus a long match is already a handful of short fixed codes. The dynamic header costs more.
        val plain = ByteArray(400) { 'A'.code.toByte() }
        val payload = roundTrip(plain)
        assertEquals(FIXED, blockType(payload))
        assertTrue(payload.size < 80, "a run of 400 should compress, was ${payload.size} bytes")
    }

    @Test
    fun skewedLiteralsUseADynamicBlock() {
        // No repeated 3-gram, so every byte is a literal. Byte 1 dominates and the fixed codes are 8 bits.
        val plain = ByteArray(254 * 3)
        var cursor = 0
        for (symbol in 2..255) {
            plain[cursor++] = 1
            plain[cursor++] = 1
            plain[cursor++] = symbol.toByte()
        }
        val payload = roundTrip(plain)
        assertEquals(DYNAMIC, blockType(payload))
        assertTrue(payload.size < plain.size, "skewed literals should compress, was ${payload.size} bytes")
    }

    @Test
    fun matchOf258AndASliceRoundTrip() {
        roundTrip(ByteArray(400) { 'Z'.code.toByte() })
        val padded = ByteArray(10) + "abcdef".encodeToByteArray() + ByteArray(4)
        roundTrip(padded, from = 10, to = 16)
    }

    @Test
    fun distance32768RoundTrips() {
        val plain = ByteArray(DeflateCodes.MAX_DISTANCE + 3)
        plain[0] = 1
        plain[1] = 2
        plain[2] = 3
        plain[DeflateCodes.MAX_DISTANCE] = 1
        plain[DeflateCodes.MAX_DISTANCE + 1] = 2
        plain[DeflateCodes.MAX_DISTANCE + 2] = 3
        val payload = roundTrip(plain)
        assertTrue(payload.size < plain.size / 4, "repeated zeros should compress, was ${payload.size} bytes")
    }

    @Test
    fun lengthAndDistanceCodesCoverTheRfcRanges() {
        for (length in 3..DeflateCodes.MAX_MATCH) {
            val index = DeflateCodes.lengthCode[length]
            val extra = length - DeflateCodes.lengthBase[index]
            assertTrue(extra >= 0 && extra < (1 shl DeflateCodes.lengthExtra[index]), "length $length")
        }
        for (distance in 1..DeflateCodes.MAX_DISTANCE) {
            val index = DeflateCodes.distanceCode[distance]
            val extra = distance - DeflateCodes.distanceBase[index]
            assertTrue(extra >= 0 && extra < (1 shl DeflateCodes.distanceExtra[index]), "distance $distance")
        }
        assertEquals(28, DeflateCodes.lengthCode[258])
        assertEquals(29, DeflateCodes.distanceCode[DeflateCodes.MAX_DISTANCE])
    }

    @Test
    fun canonicalCodesMatchTheFixedAlphabet() {
        val lengths = IntArray(288)
        for (symbol in 0..143) lengths[symbol] = 8
        for (symbol in 144..255) lengths[symbol] = 9
        for (symbol in 256..279) lengths[symbol] = 7
        for (symbol in 280..287) lengths[symbol] = 8
        val codes = IntArray(288)
        CodeLengthBuilder().codes(lengths, 288, codes)
        for (symbol in 0 until 288) {
            assertEquals(fixedLiteralCode(symbol).first, codes[symbol], "symbol $symbol")
        }
    }

    @Test
    fun lengthsStayWithinFifteenBits() {
        val freq = IntArray(32)
        freq[0] = 1
        freq[1] = 1
        var previous = 1
        var current = 1
        for (symbol in 2 until 17) {
            val next = previous + current
            freq[symbol] = next
            previous = current
            current = next
        }
        val builder = CodeLengthBuilder()
        val unlimited = IntArray(32)
        builder.lengths(freq, 17, 32, unlimited)
        assertTrue(unlimited.max() > 15, "fibonacci frequencies should need a length past 15")
        val limited = IntArray(32)
        builder.lengths(freq, 17, 15, limited)
        assertTrue(limited.max() <= 15)
        assertComplete(limited, 17)
    }

    private fun roundTrip(plain: ByteArray, from: Int = 0, to: Int = plain.size): ByteArray {
        val payload = compressor.compress(plain, from, to)
        assertNotNull(payload)
        val size = to - from
        val out = ByteArray(size)
        val written = decoder.inflate(payload, 0, payload.size, out, 0, size)
        assertEquals(size, written)
        assertContentEquals(plain.copyOfRange(from, to), out)
        return payload
    }

    private fun blockType(payload: ByteArray): Int {
        val bits = payload[0].toInt() and 0xFF
        return (bits ushr 1) and 3
    }

    private fun assertComplete(lengths: IntArray, symbolCount: Int) {
        val counts = IntArray(16)
        var used = 0
        for (symbol in 0 until symbolCount) {
            val length = lengths[symbol]
            if (length != 0) {
                counts[length]++
                used++
            }
        }
        if (used <= 1) return
        var left = 1
        for (length in 1..15) {
            left = (left shl 1) - counts[length]
            assertTrue(left >= 0, "over-subscribed")
        }
        assertEquals(0, left, "incomplete")
    }

    private companion object {
        const val FIXED = 1
        const val DYNAMIC = 2
    }
}
