package net.maizegenetics.tassel.io.deflate

import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class InflateDecoderTest {
    @Test
    fun emptyFixedStreamIsTheRfcBytes() {
        val writer = BitWriter()
        writer.writeLsb(1, 1)
        writer.writeLsb(1, 2)
        writer.writeMsb(0, 7)
        assertContentEquals(byteArrayOf(0x03, 0x00), writer.toByteArray())
        assertContentEquals(byteArrayOf(), inflate(byteArrayOf(0x03, 0x00)))
    }

    @Test
    fun emptyStoredBlock() {
        assertContentEquals(byteArrayOf(), inflate(storedBlock(byteArrayOf(), last = true)))
    }

    @Test
    fun storedBlocksIncludingTheMaximumLength() {
        assertContentEquals("hi".encodeToByteArray(), inflate(storedBlock("hi".encodeToByteArray(), last = true)))
        val big = ByteArray(65_535) { (it and 0xFF).toByte() }
        assertContentEquals(big, inflate(storedBlock(big, last = true), big.size))
    }

    @Test
    fun fixedLiteralsAndAnOverlappingMatch() {
        val hello = "Hi".encodeToByteArray()
        assertContentEquals(hello, inflate(fixedBlock(hello, last = true)))
        val overlapped = byteArrayOf('a'.code.toByte()) + ByteArray(3) { 'a'.code.toByte() }
        val stream = fixedBlock(byteArrayOf('a'.code.toByte()), last = true, match = Match(3, 1))
        assertContentEquals(overlapped, inflate(stream))
    }

    @Test
    fun length258AndDistance32768() {
        val lead = ByteArray(32_768) { ((it * 31) + 7).toByte() }
        val writer = BitWriter()
        writeFixed(writer, lead, last = true, match = Match(258, 32_768))
        val expected = ByteArray(lead.size + 258)
        lead.copyInto(expected)
        lead.copyInto(expected, lead.size, 0, 258)
        assertContentEquals(expected, inflate(writer.toByteArray(), expected.size))
    }

    @Test
    fun fastLoopCopiesAnOverlappingMatch() {
        val lead = ByteArray(10) { 'q'.code.toByte() }
        val tail = ByteArray(260) { 'z'.code.toByte() }
        val ordered = BitWriter()
        ordered.writeLsb(1, 1)
        ordered.writeLsb(1, 2)
        for (byte in lead) {
            val (code, length) = fixedLiteralCode(byte.toInt() and 0xFF)
            ordered.writeMsb(code, length)
        }
        writeFixedMatch(ordered, 20, 1)
        for (byte in tail) {
            val (code, length) = fixedLiteralCode(byte.toInt() and 0xFF)
            ordered.writeMsb(code, length)
        }
        val (endCode, endLength) = fixedLiteralCode(256)
        ordered.writeMsb(endCode, endLength)
        val expected = ByteArray(10 + 20 + tail.size)
        lead.copyInto(expected)
        for (index in 0 until 20) expected[10 + index] = expected[10 + index - 1]
        tail.copyInto(expected, 30)
        assertContentEquals(expected, inflate(ordered.toByteArray(), expected.size))
    }

    @Test
    fun dynamicBlockWithASingleDistanceCode() {
        // Plaintext "A". Literal codes: 'A' and end-of-block, one bit each. One unused distance code of length 1.
        val writer = BitWriter()
        writer.writeLsb(1, 1)
        writer.writeLsb(2, 2)
        writer.writeLsb(0, 5)
        writer.writeLsb(0, 5)
        writer.writeLsb(14, 4)
        val codeLengthLengths = IntArray(18)
        codeLengthLengths[3] = 1
        codeLengthLengths[17] = 1
        for (length in codeLengthLengths) writer.writeLsb(length, 3)
        repeat(65) { writer.writeLsb(0, 1) }
        writer.writeLsb(1, 1)
        repeat(190) { writer.writeLsb(0, 1) }
        writer.writeLsb(1, 1)
        writer.writeLsb(1, 1)
        writer.writeLsb(0, 1)
        writer.writeLsb(1, 1)
        assertContentEquals("A".encodeToByteArray(), inflate(writer.toByteArray()))
    }

    @Test
    fun twoBlocksInOneStreamAndASliceOfALargerBuffer() {
        val writer = BitWriter()
        writeStored(writer, "ab".encodeToByteArray(), last = false)
        writeFixed(writer, "c".encodeToByteArray(), last = true)
        assertContentEquals("abc".encodeToByteArray(), inflate(writer.toByteArray()))

        val stream = byteArrayOf(0x03, 0x00)
        val padded = byteArrayOf(0x7F) + stream + byteArrayOf(0x7F)
        val out = ByteArray(4) { 0x5A }
        val written = InflateDecoder().inflate(padded, 1, 3, out, 1, 2)
        assertEquals(0, written)
        assertContentEquals(byteArrayOf(0x5A, 0x5A, 0x5A, 0x5A), out)
    }

    @Test
    fun oneDecoderReusesItsTablesAcrossBlockTypes() {
        val decoder = InflateDecoder()
        val stored = storedBlock("xy".encodeToByteArray(), last = true)
        val fixed = fixedBlock("z".encodeToByteArray(), last = true)
        val dynamic = dynamicLiteralA()
        assertContentEquals("xy".encodeToByteArray(), inflate(stored, decoder = decoder))
        assertContentEquals("z".encodeToByteArray(), inflate(fixed, decoder = decoder))
        assertContentEquals("A".encodeToByteArray(), inflate(dynamic, decoder = decoder))
        val again = ByteArray(2)
        assertEquals(2, decoder.inflate(stored, 0, stored.size, again, 0, 2))
        assertContentEquals("xy".encodeToByteArray(), again)
    }

    @Test
    fun aFailedDecodeDoesNotStopTheNextCall() {
        val decoder = InflateDecoder()
        assertFailsWith<IllegalArgumentException> { inflate(byteArrayOf(0x07), decoder = decoder) }
        assertContentEquals(byteArrayOf(), inflate(byteArrayOf(0x03, 0x00), decoder = decoder))
    }

    @Test
    fun rejectsEachStrictError() {
        assertRejects(byteArrayOf(0x07), 8, "invalid block type")
        val mismatch = BitWriter()
        mismatch.writeLsb(1, 1)
        mismatch.writeLsb(0, 2)
        mismatch.align()
        mismatch.writeShortLe(1)
        mismatch.writeShortLe(0)
        mismatch.writeByte('a'.code)
        assertRejects(mismatch.toByteArray(), 8, "stored block length")
        assertRejects(overSubscribed(), 8, "over-subscribed")
        assertRejects(incomplete(), 8, "incomplete code lengths")
        assertRejects(fixedSymbol(286), 8, "invalid length symbol")
        assertRejects(fixedDistance(30), 8, "invalid distance symbol")
        assertRejects(fixedBlock(byteArrayOf(), last = true, match = Match(3, 1)), 8, "distance past the start")
        assertRejects(fixedBlock("a".encodeToByteArray(), last = true), 0, "output overflow")
        assertRejects(byteArrayOf(0x03), 8, "truncated input")
        assertRejects(byteArrayOf(0x03, 0x00, 0x00), 8, "trailing input")
        assertRejects(codeLengthRepeatWithoutAPrevious(), 8, "invalid length symbol")
        assertRejects(codeLengthRepeatPastTheEnd(), 8, "invalid length symbol")
    }

    @Test
    fun rejectsARangeOutsideEitherArray() {
        val failure = assertFailsWith<IllegalArgumentException> {
            InflateDecoder().inflate(byteArrayOf(0x03, 0x00), 0, 3, ByteArray(1), 0, 1)
        }
        assertEquals(true, failure.message?.contains("invalid inflate range"), failure.message)
    }

    private fun inflate(
        input: ByteArray,
        capacity: Int = input.size.coerceAtLeast(1) * 8 + 64,
        decoder: InflateDecoder = InflateDecoder(),
    ): ByteArray {
        val output = ByteArray(capacity)
        val written = decoder.inflate(input, 0, input.size, output, 0, capacity)
        return output.copyOf(written)
    }

    private fun assertRejects(input: ByteArray, capacity: Int, snippet: String) {
        val output = ByteArray(capacity)
        val failure = assertFailsWith<IllegalArgumentException> {
            InflateDecoder().inflate(input, 0, input.size, output, 0, capacity)
        }
        assertTrue(failure.message?.contains(snippet) == true, failure.message)
    }

    private fun fixedSymbol(symbol: Int): ByteArray {
        val writer = BitWriter()
        writer.writeLsb(1, 1)
        writer.writeLsb(1, 2)
        val (code, length) = fixedLiteralCode(symbol)
        writer.writeMsb(code, length)
        return writer.toByteArray()
    }

    private fun fixedDistance(symbol: Int): ByteArray {
        val writer = BitWriter()
        writer.writeLsb(1, 1)
        writer.writeLsb(1, 2)
        val (code, length) = fixedLiteralCode('a'.code)
        writer.writeMsb(code, length)
        val (lengthCode, lengthBits) = fixedLiteralCode(257)
        writer.writeMsb(lengthCode, lengthBits)
        writer.writeMsb(symbol, 5)
        return writer.toByteArray()
    }

    /** Every literal length is 1, so the literal tree is over-subscribed. The precode is the single symbol `1`. */
    private fun overSubscribed(): ByteArray {
        val writer = BitWriter()
        writer.writeLsb(1, 1)
        writer.writeLsb(2, 2)
        writer.writeLsb(0, 5)
        writer.writeLsb(0, 5)
        writer.writeLsb(14, 4)
        repeat(17) { writer.writeLsb(0, 3) }
        writer.writeLsb(1, 3)
        repeat(258) { writer.writeLsb(0, 1) }
        return writer.toByteArray()
    }

    /** Two literal codes of length 2, which does not fill the tree. */
    private fun incomplete(): ByteArray {
        val writer = BitWriter()
        writer.writeLsb(1, 1)
        writer.writeLsb(2, 2)
        writer.writeLsb(0, 5)
        writer.writeLsb(0, 5)
        writer.writeLsb(13, 4)
        val lengths = IntArray(17)
        lengths[3] = 1
        lengths[16] = 1
        for (length in lengths) writer.writeLsb(length, 3)
        writer.writeLsb(1, 1)
        writer.writeLsb(1, 1)
        repeat(256) { writer.writeLsb(0, 1) }
        return writer.toByteArray()
    }

    private fun codeLengthRepeatWithoutAPrevious(): ByteArray {
        val writer = BitWriter()
        writer.writeLsb(1, 1)
        writer.writeLsb(2, 2)
        writer.writeLsb(0, 5)
        writer.writeLsb(0, 5)
        writer.writeLsb(0, 4)
        writer.writeLsb(1, 3)
        writer.writeLsb(0, 3)
        writer.writeLsb(0, 3)
        writer.writeLsb(0, 3)
        writer.writeLsb(0, 1)
        return writer.toByteArray()
    }

    private fun codeLengthRepeatPastTheEnd(): ByteArray {
        val writer = BitWriter()
        writer.writeLsb(1, 1)
        writer.writeLsb(2, 2)
        writer.writeLsb(0, 5)
        writer.writeLsb(0, 5)
        writer.writeLsb(0, 4)
        writer.writeLsb(0, 3)
        writer.writeLsb(0, 3)
        writer.writeLsb(1, 3)
        writer.writeLsb(0, 3)
        writer.writeLsb(0, 1)
        writer.writeLsb(127, 7)
        writer.writeLsb(0, 1)
        writer.writeLsb(127, 7)
        return writer.toByteArray()
    }

    private fun dynamicLiteralA(): ByteArray {
        val writer = BitWriter()
        writer.writeLsb(1, 1)
        writer.writeLsb(2, 2)
        writer.writeLsb(0, 5)
        writer.writeLsb(0, 5)
        writer.writeLsb(14, 4)
        val codeLengthLengths = IntArray(18)
        codeLengthLengths[3] = 1
        codeLengthLengths[17] = 1
        for (length in codeLengthLengths) writer.writeLsb(length, 3)
        repeat(65) { writer.writeLsb(0, 1) }
        writer.writeLsb(1, 1)
        repeat(190) { writer.writeLsb(0, 1) }
        writer.writeLsb(1, 1)
        writer.writeLsb(1, 1)
        writer.writeLsb(0, 1)
        writer.writeLsb(1, 1)
        return writer.toByteArray()
    }
}
