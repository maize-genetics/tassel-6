package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.io.deflate.BGZF_BLOCK_SIZE
import net.maizegenetics.tassel.io.deflate.BGZF_EOF_BLOCK
import net.maizegenetics.tassel.io.deflate.bgzfBlocks
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.sha256
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
