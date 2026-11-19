package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.io.deflate.BGZF_BLOCK_SIZE
import net.maizegenetics.tassel.io.deflate.BGZF_EOF_BLOCK
import net.maizegenetics.tassel.io.deflate.bgzfBlocks
import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNotEquals
import kotlin.test.assertTrue

/**
 * Stored BGZF against htslib 1.25's block splits (ADR-0022).
 *
 * The oracle files are compressed, so the comparison is the uncompressed size of each member,
 * including the empty EOF member, and the inflated bytes. Stored members are larger on disk.
 */
class BgzfWriterTest {
    @Test
    fun emptyInputIsTheEofMarker() {
        val written = writeBlocks(ByteArray(0))
        assertContentEquals(BGZF_EOF_BLOCK, written)
        assertEquals(28, written.size)
    }

    @Test
    fun smallVcfIsAHeaderBlockAndOneBodyBlock() = runTest {
        val plain = TestData.readBytes(SOURCE_VCF)
        val written = writeBlocks(plain)
        assertEquals(listOf(165, 70), isizes(TestData.readBytes(NO_EOF)))
        assertEquals(listOf(165, 70, 0), isizes(written))
        assertContentEquals(plain, roundTrip(written))
        assertFraming(written, firstUncompressed = 165)
    }

    @Test
    fun longLineAndLongHeaderMatchTheLayoutFixtures() = runTest {
        assertLayout(LONG_LINE, listOf(113, 65_280, 258, 0), splitInsideLine = true)
        assertLayout(LONG_HEADER, listOf(65_251, 440, 29, 0), splitInsideLine = false)
    }

    @Test
    fun tutorialMatchesTheHtslibBlockList() = runTest {
        val plain = TestData.readBytes(PLAIN_VCF)
        val written = writeBlocks(plain)
        val expected = isizes(TestData.readBytes(TUTORIAL_GZ))
        assertEquals(288, expected.size)
        assertEquals(19_000, expected.first())
        assertEquals(0, expected.last())
        assertEquals(expected, isizes(written))
        assertContentEquals(plain, roundTrip(written))
        assertContentEquals(BGZF_EOF_BLOCK, written.copyOfRange(written.size - 28, written.size))
        val hash = sha256(written)
        assertEquals(TUTORIAL_SHA256, hash, hash)
    }

    @Test
    fun tutorialBytesDoNotDependOnTheWriteSize() {
        val plain = TestData.readBytes(PLAIN_VCF)
        val whole = writeBlocks(plain)
        assertContentEquals(whole, writeBlocks(plain, chunk = 4096))
        assertContentEquals(whole, writeBlocks(plain, chunk = 1))
    }

    @Test
    fun flushDoesNotEndABlock() {
        val plain = inflate(TestData.readBytes(LONG_LINE))
        val sink = Buffer()
        val writer = BgzfWriter(sink, BgzfLayout.Vcf)
        writer.write(plain, 0, plain.size - 10)
        val told = writer.tell()
        val framed = sink.size
        assertTrue(framed > 0L, "a complete block should already be framed")
        writer.flush()
        assertEquals(told, writer.tell())
        assertEquals(framed, sink.size)
        writer.write(plain, plain.size - 10, plain.size)
        writer.close()
        assertContentEquals(writeBlocks(plain), sink.readByteArray())
    }

    @Test
    fun closeLeavesTheSinkOpen() {
        val sink = Buffer()
        val writer = BgzfWriter(sink, BgzfLayout.Vcf)
        writer.write("##fileformat=VCFv4.2\nchr1\t1\n".encodeToByteArray())
        writer.close()
        val framed = sink.size
        sink.writeByte(9)
        assertEquals(framed + 1, sink.size)
        writer.close()
        assertEquals(framed + 1, sink.size)
        assertFailsWith<IllegalStateException> { writer.write(byteArrayOf(1)) }
        assertFailsWith<IllegalStateException> { writer.tell() }
        assertFailsWith<IllegalStateException> { writer.flush() }
        val bytes = sink.readByteArray()
        assertEquals(9, bytes.last().toInt())
        assertContentEquals(BGZF_EOF_BLOCK, bytes.copyOfRange(bytes.size - 29, bytes.size - 1))
    }

    @Test
    fun endHeaderEndsTheVcfHeaderOnce() {
        val plain = TestData.readBytes(SOURCE_VCF)
        val header = headerLength(plain)
        assertEquals(165, header)
        val sink = Buffer()
        val writer = BgzfWriter(sink, BgzfLayout.Vcf)
        writer.write(plain, 0, header)
        writer.endHeader()
        val afterHeader = bgzfVirtualOffset(165 + STORED_OVERHEAD.toLong(), 0)
        assertEquals(afterHeader, writer.tell())
        writer.endHeader()
        assertEquals(afterHeader, writer.tell())
        writer.write(plain, header, header + 1)
        assertEquals(bgzfVirtualOffset(165 + STORED_OVERHEAD.toLong(), 1), writer.tell())
        writer.endHeader()
        assertEquals(bgzfVirtualOffset(165 + STORED_OVERHEAD.toLong(), 1), writer.tell())
        writer.write(plain, header + 1, plain.size)
        writer.close()
        assertContentEquals(writeBlocks(plain), sink.readByteArray())
    }

    @Test
    fun bcfBlocksHold65280BytesAndMaySplitARecord() {
        val packed = writeBlocks(ByteArray(BGZF_BLOCK_SIZE + 100), BgzfLayout.Bcf)
        assertEquals(listOf(BGZF_BLOCK_SIZE, 100, 0), isizes(packed))

        val header = ByteArray(100) { 'H'.code.toByte() }
        val body = ByteArray(BGZF_BLOCK_SIZE + 50)
        var index = 0
        while (index < body.size) {
            val line = minOf(7, body.size - index)
            for (offset in 0 until line) body[index + offset] = if (offset == 6) NEWLINE else 'a'.code.toByte()
            index += line
        }
        assertNotEquals(NEWLINE, body[BGZF_BLOCK_SIZE - 1])
        val sink = Buffer()
        val writer = BgzfWriter(sink, BgzfLayout.Bcf)
        writer.write(header)
        writer.endHeader()
        writer.write(body)
        writer.close()
        assertEquals(listOf(100, BGZF_BLOCK_SIZE, 50, 0), isizes(sink.readByteArray()))

        val longHeader = Buffer()
        BgzfWriter(longHeader, BgzfLayout.Bcf).use { longWriter ->
            longWriter.write(ByteArray(BGZF_BLOCK_SIZE + 40))
            longWriter.endHeader()
            longWriter.write(ByteArray(10))
        }
        assertEquals(listOf(BGZF_BLOCK_SIZE, 40, 10, 0), isizes(longHeader.readByteArray()))
    }

    @Test
    fun tellTracksTheOpenBcfBlock() {
        val sink = Buffer()
        val writer = BgzfWriter(sink, BgzfLayout.Bcf)
        assertEquals(0L, writer.tell())
        writer.write(ByteArray(10))
        assertEquals(bgzfVirtualOffset(0, 10), writer.tell())
        writer.write(ByteArray(BGZF_BLOCK_SIZE - 10))
        val member = (BGZF_BLOCK_SIZE + STORED_OVERHEAD).toLong()
        assertEquals(member, sink.size)
        assertEquals(bgzfVirtualOffset(member, 0), writer.tell())
        writer.write(ByteArray(5))
        val told = writer.tell()
        assertEquals(bgzfVirtualOffset(member, 5), told)
        writer.flush()
        assertEquals(told, writer.tell())
        assertEquals(member, sink.size)
        writer.close()
    }

    @Test
    fun compressorFallbackStoresTheBlock() {
        val plain = TestData.readBytes(SOURCE_VCF)
        val stored = writeBlocks(plain)
        var calls = 0
        val counting = BgzfCompressor { _, from, to ->
            calls++
            assertTrue(to > from)
            null
        }
        assertContentEquals(stored, writeBlocks(plain, compressor = counting))
        assertEquals(2, calls, "the EOF member is not compressed")
        val oversized = BgzfCompressor { _, _, _ -> ByteArray(70_000) }
        assertContentEquals(stored, writeBlocks(plain, compressor = oversized))
    }

    @Test
    fun sha256OfAbcMatchesTheVector() {
        assertEquals(
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
            sha256("abc".encodeToByteArray()),
        )
        assertEquals(
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
            sha256(ByteArray(0)),
        )
    }

    private suspend fun assertLayout(path: String, sizes: List<Int>, splitInsideLine: Boolean) {
        val oracle = TestData.readBytes(path)
        val plain = inflate(oracle)
        val written = writeBlocks(plain)
        assertEquals(sizes, isizes(oracle))
        assertEquals(sizes, isizes(written))
        assertContentEquals(plain, roundTrip(written))
        if (splitInsideLine) {
            val middle = plain.copyOfRange(sizes[0], sizes[0] + sizes[1])
            assertEquals(BGZF_BLOCK_SIZE, middle.size)
            assertNotEquals(NEWLINE, middle.last())
        }
    }

    private suspend fun roundTrip(stored: ByteArray): ByteArray {
        val fromReader = ByteArrayRandomAccessSource(stored).use { source ->
            BgzfReader(source).use { it.readAll() }
        }
        val fromGzip = inflate(stored)
        assertContentEquals(fromReader, fromGzip)
        return fromReader
    }

    private fun inflate(bytes: ByteArray): ByteArray =
        Buffer().also { it.write(bytes) }.decompressing().source.readByteArray()

    private fun isizes(file: ByteArray): List<Int> = bgzfBlocks(file).map { it.uncompressedSize }

    private fun assertFraming(written: ByteArray, firstUncompressed: Int) {
        assertEquals(0x1f, written[0].toInt() and 0xFF)
        assertEquals(0x8b, written[1].toInt() and 0xFF)
        assertEquals(8, written[2].toInt() and 0xFF)
        assertEquals(4, written[3].toInt() and 0xFF)
        assertEquals(0, written[4].toInt() or written[5].toInt() or written[6].toInt() or written[7].toInt())
        assertEquals(0, written[8].toInt())
        assertEquals(0xFF, written[9].toInt() and 0xFF)
        assertEquals(6, written[10].toInt() and 0xFF)
        assertEquals(0, written[11].toInt())
        assertEquals('B'.code, written[12].toInt() and 0xFF)
        assertEquals('C'.code, written[13].toInt() and 0xFF)
        assertEquals(2, written[14].toInt() and 0xFF)
        val blockSize = u16(written, 16)
        assertEquals(firstUncompressed + STORED_OVERHEAD - 1, blockSize)
        assertEquals(1, written[18].toInt() and 0xFF)
        assertEquals(firstUncompressed, u16(written, 19))
        assertEquals(firstUncompressed.inv() and 0xFFFF, u16(written, 21))
        assertContentEquals(BGZF_EOF_BLOCK, written.copyOfRange(written.size - 28, written.size))
    }

    private fun headerLength(plain: ByteArray): Int {
        var index = 0
        while (index < plain.size && (plain[index] == HASH || plain[index] == AT)) {
            while (index < plain.size && plain[index] != NEWLINE) index++
            if (index < plain.size) index++
        }
        return index
    }

    private suspend fun BgzfReader.readAll(): ByteArray {
        val buffer = Buffer()
        val chunk = ByteArray(64 * 1024)
        while (true) {
            val count = read(chunk, 0, chunk.size)
            if (count == 0) break
            buffer.write(chunk, 0, count)
        }
        return buffer.readByteArray()
    }

    private companion object {
        const val PLAIN_VCF = "tutorial/mdp_genotype_w_depth.vcf"
        const val TUTORIAL_GZ = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
        const val LONG_LINE = "htslib/bgzf-layout/long-line.vcf.gz"
        const val LONG_HEADER = "htslib/bgzf-layout/long-header.vcf.gz"
        const val SOURCE_VCF = "htslib/variants/source.vcf"
        const val NO_EOF = "htslib/variants/no-eof.vcf.gz"
        const val STORED_OVERHEAD = 31

        /** Stored BGZF of `tutorial/mdp_genotype_w_depth.vcf`. The same bytes on every target. */
        const val TUTORIAL_SHA256 = "898755e364dcfd7b5acd8c32a3a981a5d4206d721fae9f4c5cde16b5ea8256a0"

        val NEWLINE = '\n'.code.toByte()
        val HASH = '#'.code.toByte()
        val AT = '@'.code.toByte()

        fun writeBlocks(
            plain: ByteArray,
            layout: BgzfLayout = BgzfLayout.Vcf,
            chunk: Int = plain.size.coerceAtLeast(1),
            compressor: BgzfCompressor = BgzfCompressor.Stored,
        ): ByteArray {
            val sink = Buffer()
            BgzfWriter(sink, layout, compressor).use { writer ->
                var offset = 0
                while (offset < plain.size) {
                    val end = minOf(plain.size, offset + chunk)
                    writer.write(plain, offset, end)
                    offset = end
                }
            }
            return sink.readByteArray()
        }

        fun u16(bytes: ByteArray, at: Int): Int =
            (bytes[at].toInt() and 0xFF) or ((bytes[at + 1].toInt() and 0xFF) shl 8)
    }
}

private fun sha256(message: ByteArray): String {
    val bitLength = message.size.toLong() shl 3
    val zeroPad = (56 - (message.size + 1) % 64 + 64) % 64
    val data = ByteArray(message.size + 1 + zeroPad + 8)
    message.copyInto(data)
    data[message.size] = 0x80.toByte()
    for (shift in 0 until 8) {
        data[data.size - 1 - shift] = ((bitLength ushr (8 * shift)) and 0xFF).toByte()
    }
    val words = IntArray(64)
    var h0 = 0x6a09e667
    var h1 = 0xbb67ae85.toInt()
    var h2 = 0x3c6ef372
    var h3 = 0xa54ff53a.toInt()
    var h4 = 0x510e527f
    var h5 = 0x9b05688c.toInt()
    var h6 = 0x1f83d9ab.toInt()
    var h7 = 0x5be0cd19.toInt()
    var offset = 0
    while (offset < data.size) {
        for (index in 0 until 16) {
            val at = offset + index * 4
            words[index] = ((data[at].toInt() and 0xFF) shl 24) or
                ((data[at + 1].toInt() and 0xFF) shl 16) or
                ((data[at + 2].toInt() and 0xFF) shl 8) or
                (data[at + 3].toInt() and 0xFF)
        }
        for (index in 16 until 64) {
            val earlier = words[index - 15]
            val recent = words[index - 2]
            val small = earlier.rotateRight(7) xor earlier.rotateRight(18) xor (earlier ushr 3)
            val large = recent.rotateRight(17) xor recent.rotateRight(19) xor (recent ushr 10)
            words[index] = words[index - 16] + small + words[index - 7] + large
        }
        var a = h0
        var b = h1
        var c = h2
        var d = h3
        var e = h4
        var f = h5
        var g = h6
        var h = h7
        for (index in 0 until 64) {
            val sum1 = e.rotateRight(6) xor e.rotateRight(11) xor e.rotateRight(25)
            val choose = (e and f) xor (e.inv() and g)
            val first = h + sum1 + choose + SHA256_K[index] + words[index]
            val sum0 = a.rotateRight(2) xor a.rotateRight(13) xor a.rotateRight(22)
            val majority = (a and b) xor (a and c) xor (b and c)
            h = g
            g = f
            f = e
            e = d + first
            d = c
            c = b
            b = a
            a = first + sum0 + majority
        }
        h0 += a
        h1 += b
        h2 += c
        h3 += d
        h4 += e
        h5 += f
        h6 += g
        h7 += h
        offset += 64
    }
    return hex(h0) + hex(h1) + hex(h2) + hex(h3) + hex(h4) + hex(h5) + hex(h6) + hex(h7)
}

private fun hex(value: Int): String {
    val digits = "0123456789abcdef"
    val chars = CharArray(8)
    for (index in 0 until 8) {
        chars[7 - index] = digits[(value ushr (4 * index)) and 0xF]
    }
    return chars.concatToString()
}

private val SHA256_K: IntArray = intArrayOf(
    0x428a2f98, 0x71374491, 0xb5c0fbcf.toInt(), 0xe9b5dba5.toInt(),
    0x3956c25b, 0x59f111f1, 0x923f82a4.toInt(), 0xab1c5ed5.toInt(),
    0xd807aa98.toInt(), 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe.toInt(), 0x9bdc06a7.toInt(), 0xc19bf174.toInt(),
    0xe49b69c1.toInt(), 0xefbe4786.toInt(), 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152.toInt(), 0xa831c66d.toInt(), 0xb00327c8.toInt(), 0xbf597fc7.toInt(),
    0xc6e00bf3.toInt(), 0xd5a79147.toInt(), 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e.toInt(), 0x92722c85.toInt(),
    0xa2bfe8a1.toInt(), 0xa81a664b.toInt(), 0xc24b8b70.toInt(), 0xc76c51a3.toInt(),
    0xd192e819.toInt(), 0xd6990624.toInt(), 0xf40e3585.toInt(), 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814.toInt(), 0x8cc70208.toInt(),
    0x90befffa.toInt(), 0xa4506ceb.toInt(), 0xbef9a3f7.toInt(), 0xc67178f2.toInt(),
)
