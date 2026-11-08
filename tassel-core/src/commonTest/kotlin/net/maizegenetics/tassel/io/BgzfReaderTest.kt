package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.io.deflate.GZIP_FIXED_HEADER
import net.maizegenetics.tassel.io.deflate.parseBgzfPayload
import net.maizegenetics.tassel.io.deflate.parseGzipMemberHeader
import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Virtual-offset seeks against the htslib 1.25 fixtures (ADR-0021).
 *
 * The `.csi` walker is test-local. The index reader is #88.
 */
class BgzfReaderTest {
    @Test
    fun tutorialBytesMatchThePlainVcf() = runTest {
        val expected = TestData.readBytes(PLAIN_VCF)
        assertContentEquals(expected, readAll(TUTORIAL_GZ))
        assertEquals(emptyList(), warningsAfter(TUTORIAL_GZ))
        ByteArrayRandomAccessSource(TestData.readBytes(TUTORIAL_GZ)).use { source ->
            BgzfReader(source).use { reader ->
                assertContentEquals(expected, reader.readAll())
                assertEquals(emptyList(), reader.warnings)
            }
        }
    }

    @Test
    fun wheatAndLargePositionMatchSequentialGzip() = runTest {
        assertContentEquals(sequential(WHEAT_GZ), readAll(WHEAT_GZ))
        assertContentEquals(sequential(LARGE_POS_GZ), readAll(LARGE_POS_GZ))
    }

    @Test
    fun seekingToEachTutorialLineReadsThatLine() = runTest {
        assertEachLineSeek(TUTORIAL_GZ, TestData.readBytes(PLAIN_VCF))
    }

    @Test
    fun seekingToEachWheatAndLargePositionLineReadsThatLine() = runTest {
        assertEachLineSeek(WHEAT_GZ, sequential(WHEAT_GZ))
        assertEachLineSeek(LARGE_POS_GZ, sequential(LARGE_POS_GZ))
    }

    @Test
    fun csiOffsetsLandOnALineOrTheEnd() = runTest {
        assertCsiOffsets(TUTORIAL_GZ, minShift = 14, depth = 8)
        assertCsiOffsets(WHEAT_GZ, minShift = 14, depth = 6)
        assertCsiOffsets(LARGE_POS_GZ, minShift = 14, depth = 6)
    }

    @Test
    fun uoffsetAtBlockLengthContinuesAtTheNextBlock() = runTest {
        val block = firstShortBlock(TestData.readBytes(TUTORIAL_GZ))
        withReader(TUTORIAL_GZ) { reader ->
            val endOfBlock = bgzfVirtualOffset(block.start.toLong(), block.length)
            reader.seek(endOfBlock)
            val next = bgzfVirtualOffset(block.next.toLong(), 0)
            assertEquals(next, reader.tell())
            val line = reader.readLine()
            reader.seek(next)
            assertEquals(line, reader.readLine())

            reader.seek(bgzfVirtualOffset(block.start.toLong(), 0))
            val whole = ByteArray(block.length)
            assertEquals(block.length, reader.read(whole, 0, block.length))
            assertEquals(next, reader.tell())
        }
    }

    @Test
    fun seekRejectsACoffsetThatIsNotABlockStart() = runTest {
        val virtual = bgzfVirtualOffset(1, 0)
        val failure = assertFailsWith<IllegalArgumentException> {
            withReader(TUTORIAL_GZ) { it.seek(virtual) }
        }
        assertContains(failure.message ?: "", "not a block start")
        assertContains(failure.message ?: "", virtual.toString())
    }

    @Test
    fun seekRejectsAUoffsetPastTheBlock() = runTest {
        val block = firstShortBlock(TestData.readBytes(TUTORIAL_GZ))
        val virtual = bgzfVirtualOffset(block.start.toLong(), block.length + 1)
        val failure = assertFailsWith<IllegalArgumentException> {
            withReader(TUTORIAL_GZ) { it.seek(virtual) }
        }
        assertContains(failure.message ?: "", "past the block length")
        assertContains(failure.message ?: "", virtual.toString())
    }

    @Test
    fun skipsTheEmptyBlockInTheMiddleOfJoinedBgzip() = runTest {
        assertContentEquals(TestData.readBytes(SOURCE_VCF), readAll(MID_EOF))
        assertEquals(emptyList(), warningsAfter(MID_EOF))
    }

    @Test
    fun warnsWhenTheEofBlockIsMissing() = runTest {
        assertContentEquals(TestData.readBytes(SOURCE_VCF), readAll(NO_EOF))
        assertEquals(listOf(MISSING_BGZF_EOF), warningsAfter(NO_EOF))
    }

    @Test
    fun rejectsABadCrc() = runTest {
        val failure = assertFailsWith<IllegalArgumentException> { readAll(BAD_CRC) }
        assertContains(failure.message ?: "", "CRC")
    }

    @Test
    fun reportsProgressBeforeTheFileIsExhausted() = runTest {
        val compressed = TestData.readBytes(TUTORIAL_GZ)
        var reports = 0
        var first = -1L
        ByteArrayRandomAccessSource(compressed).use { source ->
            BgzfReader(source) { progress ->
                reports++
                if (first < 0) first = progress.completed
            }.use { reader ->
                reader.readAll()
            }
        }
        assertTrue(reports > 1, "reports $reports")
        assertTrue(first in 1 until compressed.size.toLong(), "first report at $first of ${compressed.size}")
    }

    private suspend fun assertEachLineSeek(path: String, plain: ByteArray) {
        val lines = ArrayList<Pair<Long, String>>()
        withReader(path) { reader ->
            while (true) {
                val at = reader.tell()
                val line = reader.readLine() ?: break
                lines.add(at to line)
            }
            val text = lines.joinToString("\n") { it.second }
            assertEquals(plain.decodeToString().removeSuffix("\n"), text, path)
            for ((at, line) in lines) {
                reader.seek(at)
                assertEquals(at, reader.tell(), "$path at $at")
                assertEquals(line, reader.readLine(), "$path at $at")
            }
        }
    }

    private suspend fun assertCsiOffsets(path: String, minShift: Int, depth: Int) {
        val starts = HashSet<Long>()
        withReader(path) { reader ->
            while (true) {
                val at = reader.tell()
                if (reader.readLine() == null) break
                starts.add(at)
            }
        }
        val index = csiOffsets(sequential("$path.csi"))
        assertEquals(minShift, index.minShift, path)
        assertEquals(depth, index.depth, path)
        assertTrue(index.offsets.isNotEmpty(), path)
        val one = ByteArray(1)
        withReader(path) { reader ->
            for (offset in index.offsets) {
                reader.seek(offset)
                if (reader.tell() in starts) continue
                assertEquals(0, reader.read(one, 0, 1), "$path offset $offset")
            }
        }
    }

    private suspend fun readAll(path: String): ByteArray = withReader(path) { it.readAll() }

    private suspend fun warningsAfter(path: String): List<String> = withReader(path) { reader ->
        reader.readAll()
        reader.warnings
    }

    private suspend fun <T> withReader(path: String, body: suspend (BgzfReader) -> T): T =
        openRandomAccessFile("${TestData.root}/$path").use { source ->
            BgzfReader(source).use { body(it) }
        }

    private companion object {
        const val PLAIN_VCF = "tutorial/mdp_genotype_w_depth.vcf"
        const val TUTORIAL_GZ = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
        const val WHEAT_GZ = "htslib/wheat/wheat.vcf.gz"
        const val LARGE_POS_GZ = "htslib/large-pos/large_pos.vcf.gz"
        const val SOURCE_VCF = "htslib/variants/source.vcf"
        const val MID_EOF = "htslib/variants/mid-eof.vcf.gz"
        const val NO_EOF = "htslib/variants/no-eof.vcf.gz"
        const val BAD_CRC = "htslib/variants/bad-crc.vcf.gz"

        fun sequential(path: String): ByteArray =
            Buffer().also { it.write(TestData.readBytes(path)) }.decompressing().source.readByteArray()

        suspend fun BgzfReader.readAll(): ByteArray {
            val buffer = Buffer()
            val chunk = ByteArray(64 * 1024)
            while (true) {
                val count = read(chunk, 0, chunk.size)
                if (count == 0) break
                buffer.write(chunk, 0, count)
            }
            return buffer.readByteArray()
        }

        fun firstShortBlock(file: ByteArray): ShortBlock {
            var start = 0
            while (start + GZIP_FIXED_HEADER <= file.size) {
                val header = parseGzipMemberHeader(file, start, file.size - start)
                check(header.blockSize >= 0) { "plain gzip at $start" }
                val total = header.blockSize + 1
                val payload = parseBgzfPayload(file, start, header)
                val next = start + total
                if (payload.uncompressedSize in 1..65534 && next + GZIP_FIXED_HEADER <= file.size) {
                    val nextHeader = parseGzipMemberHeader(file, next, file.size - next)
                    val nextPayload = parseBgzfPayload(file, next, nextHeader)
                    if (nextPayload.uncompressedSize > 0) {
                        return ShortBlock(start, payload.uncompressedSize, next)
                    }
                }
                start = next
            }
            error("no short BGZF block followed by another data block")
        }

        fun csiOffsets(bytes: ByteArray): CsiOffsets {
            require(bytes.size >= 16 && bytes[0] == 'C'.code.toByte() && bytes[1] == 'S'.code.toByte()) {
                "CSI magic"
            }
            require(bytes[2] == 'I'.code.toByte() && bytes[3] == 1.toByte()) { "CSI magic" }
            var cursor = 4
            val minShift = i32(bytes, cursor)
            cursor += 4
            val depth = i32(bytes, cursor)
            cursor += 4
            val aux = i32(bytes, cursor)
            cursor += 4
            require(aux >= 0) { "CSI aux length $aux" }
            cursor = skip(bytes, cursor, aux)
            val references = i32(bytes, cursor)
            cursor += 4
            require(references in 0..1_000_000) { "CSI references $references" }
            // htslib's meta bin is ((1 << (3 * depth + 3)) - 1) / 7 + 1. Its second chunk is
            // (n_mapped, n_unmapped), not a virtual offset. The first chunk is the contig's span.
            val meta = ((1L shl (3 * depth + 3)) - 1) / 7 + 1
            val offsets = ArrayList<Long>()
            repeat(references) {
                val bins = i32(bytes, cursor)
                cursor += 4
                require(bins in 0..1_000_000) { "CSI bins $bins" }
                repeat(bins) {
                    val bin = u32(bytes, cursor)
                    cursor += 4
                    offsets.add(u64(bytes, cursor))
                    cursor += 8
                    val chunks = i32(bytes, cursor)
                    cursor += 4
                    require(chunks in 0..1_000_000) { "CSI chunks $chunks" }
                    repeat(chunks) { chunk ->
                        val beg = u64(bytes, cursor)
                        cursor += 8
                        val end = u64(bytes, cursor)
                        cursor += 8
                        if (!(bin == meta && chunk == 1)) {
                            offsets.add(beg)
                            offsets.add(end)
                        }
                    }
                }
            }
            return CsiOffsets(minShift, depth, offsets)
        }

        fun skip(bytes: ByteArray, cursor: Int, count: Int): Int {
            require(cursor >= 0 && count >= 0 && cursor <= bytes.size - count) { "CSI truncated at $cursor" }
            return cursor + count
        }

        fun u32(bytes: ByteArray, at: Int): Long {
            skip(bytes, at, 4)
            return (bytes[at].toLong() and 0xFF) or
                ((bytes[at + 1].toLong() and 0xFF) shl 8) or
                ((bytes[at + 2].toLong() and 0xFF) shl 16) or
                ((bytes[at + 3].toLong() and 0xFF) shl 24)
        }

        fun i32(bytes: ByteArray, at: Int): Int {
            skip(bytes, at, 4)
            return (bytes[at].toInt() and 0xFF) or
                ((bytes[at + 1].toInt() and 0xFF) shl 8) or
                ((bytes[at + 2].toInt() and 0xFF) shl 16) or
                ((bytes[at + 3].toInt() and 0xFF) shl 24)
        }

        fun u64(bytes: ByteArray, at: Int): Long {
            skip(bytes, at, 8)
            var value = 0L
            for (index in 0 until 8) {
                value = value or ((bytes[at + index].toLong() and 0xFF) shl (8 * index))
            }
            return value
        }
    }
}

private class ShortBlock(val start: Int, val length: Int, val next: Int)

private class CsiOffsets(val minShift: Int, val depth: Int, val offsets: List<Long>)
