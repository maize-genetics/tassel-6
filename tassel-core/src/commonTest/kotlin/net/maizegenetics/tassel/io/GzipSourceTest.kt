package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.RawSource
import kotlinx.io.buffered
import kotlinx.io.readByteArray
import kotlinx.io.readString
import kotlinx.io.writeString
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.io.deflate.Crc32
import net.maizegenetics.tassel.io.deflate.bgzfBlocks
import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Sequential BGZF reads against the htslib 1.25 fixtures (ADR-0021) and the TASSEL 5 HapMap export.
 */
class GzipSourceTest {
    @Test
    fun writesTheTassel5HapMapExportOfTheTutorialVcfGz() = runTest {
        val opened = open(TUTORIAL_GZ)
        val written = Buffer().also { HapMap.write(Vcf.read(opened.source), it) }.readString().lines()
        assertEquals(emptyList(), opened.warnings, "the tutorial bgzip ends on an EOF block")

        val goldenLines = TestData.readText(GOLDEN).replace("\r\n", "\n").lines()
        val skipped = goldenLines.takeWhile { it.startsWith("##") }
        assertEquals(281, skipped.size, "skipped ##SAMPLE= lines")
        assertTrue(skipped.all { it.startsWith("##SAMPLE=<ID=") }, "every skipped line is a ##SAMPLE= annotation")
        val golden = goldenLines.drop(skipped.size)
        assertTrue(golden.first().startsWith("rs#\talleles\tchrom\tpos\t"), "first compared line")
        assertEquals(golden.size, written.size, "number of lines")
        written.forEachIndexed { index, line ->
            if (line != golden[index]) reportLineDifference(skipped.size + index, golden[index], line)
        }
    }

    @Test
    fun readsBgzippedHalvesJoinedWithCat() {
        val (bytes, warnings) = inflated(MID_EOF)
        assertContentEquals(TestData.readBytes(SOURCE), bytes)
        assertEquals(emptyList(), warnings)
    }

    @Test
    fun warnsWhenTheEofBlockIsMissing() {
        val opened = open(NO_EOF)
        assertEquals(emptyList(), opened.warnings, "the warning is recorded when the stream ends")
        assertContentEquals(TestData.readBytes(SOURCE), opened.source.readByteArray())
        assertEquals(listOf(MISSING_BGZF_EOF), opened.warnings)
    }

    @Test
    fun rejectsATruncatedBlock() {
        assertFailure(TRUNCATED, "truncated")
    }

    @Test
    fun rejectsABadCrc() {
        assertFailure(BAD_CRC, "CRC")
    }

    @Test
    fun rejectsAnIsizeThatDoesNotMatch() {
        val file = TestData.readBytes(NO_EOF).copyOf()
        val last = bgzfBlocks(file).last { it.uncompressedSize > 0 }
        val at = last.to + 4
        val isize = u32(file, at)
        putU32(file, at, isize + 1)
        val failure = assertFailsWith<IllegalArgumentException> {
            Buffer().also { it.write(file) }.decompressing().source.readByteArray()
        }
        assertContains(failure.message ?: "", "ISIZE")
    }

    @Test
    fun rejectsPlainGzip() {
        assertFailure(PLAIN, "tassel-6#85")
        assertFailure(CONCATENATED, "tassel-6#85")
    }

    @Test
    fun reportsProgressBeforeTheCompressedInputIsExhausted() = runTest {
        val compressed = TestData.readBytes(TUTORIAL_GZ)
        val counting = CountingRawSource(compressed)
        val opened = counting.buffered().decompressing()
        var positionAtFirstReport = -1
        Vcf.read(opened.source) {
            if (positionAtFirstReport < 0) positionAtFirstReport = counting.position
        }
        assertTrue(
            positionAtFirstReport in 1 until compressed.size,
            "first progress report at compressed byte $positionAtFirstReport of ${compressed.size}",
        )
    }

    @Test
    fun readsABgzippedHapMap() = runTest {
        val text = hapMapText()
        val gz = bgzfMember(text.encodeToByteArray()) + bgzfMember(ByteArray(0))
        val opened = Buffer().also { it.write(gz) }.decompressing()
        val table = HapMap.read(opened.source)
        assertEquals(listOf("T1"), table.taxa.names)
        assertEquals(1, table.numberOfSites)
        assertEquals(Nucleotide.genotype("AA"), table.genotype(0, 0))
        assertEquals(emptyList(), opened.warnings)
    }

    @Test
    fun skipsAnEmptyBlockBetweenHalves() {
        val text = "abcdefghij".encodeToByteArray()
        val gz = bgzfMember(text.copyOfRange(0, 4)) +
            bgzfMember(ByteArray(0)) +
            bgzfMember(text.copyOfRange(4, text.size)) +
            bgzfMember(ByteArray(0))
        val (bytes, warnings) = inflated(gz)
        assertContentEquals(text, bytes)
        assertEquals(emptyList(), warnings)
    }

    @Test
    fun warnsWhenASyntheticMemberHasNoEofBlock() {
        val text = "abcdef".encodeToByteArray()
        val (bytes, warnings) = inflated(bgzfMember(text))
        assertContentEquals(text, bytes)
        assertEquals(listOf(MISSING_BGZF_EOF), warnings)
    }

    @Test
    fun skipsANameCommentAndHeaderCrc() {
        val text = "hello".encodeToByteArray()
        val gz = bgzfMember(text, name = "sample", comment = "note", headerCrc = true) + bgzfMember(ByteArray(0))
        val (bytes, warnings) = inflated(gz)
        assertContentEquals(text, bytes)
        assertEquals(emptyList(), warnings)
    }

    @Test
    fun leavesAnUncompressedSourceUnchanged() {
        val text = "not gzip"
        val opened = Buffer().also { it.writeString(text) }.decompressing()
        assertEquals(text, opened.source.readString())
        assertEquals(emptyList(), opened.warnings)
    }

    private fun open(path: String): Decompressed = Buffer().also { it.write(TestData.readBytes(path)) }.decompressing()

    private fun inflated(path: String): Pair<ByteArray, List<String>> {
        val opened = open(path)
        return opened.source.readByteArray() to opened.warnings
    }

    private fun inflated(bytes: ByteArray): Pair<ByteArray, List<String>> {
        val opened = Buffer().also { it.write(bytes) }.decompressing()
        return opened.source.readByteArray() to opened.warnings
    }

    private fun assertFailure(path: String, snippet: String) {
        val failure = assertFailsWith<IllegalArgumentException> { inflated(path) }
        assertContains(failure.message ?: "", snippet)
    }

    private fun hapMapText(): String {
        val header = (HapMap.FIXED_COLUMNS + "T1").joinToString("\t")
        val site = listOf("s1", "A/C", "1", "10", "+", "NA", "NA", "NA", "NA", "NA", "NA", "A").joinToString("\t")
        return "$header\n$site\n"
    }

    private fun reportLineDifference(line: Int, expected: String, actual: String): Nothing {
        val column = expected.commonPrefixWith(actual).length
        throw AssertionError(
            "Line ${line + 1} differs at column ${column + 1}:\n" +
                "  expected: ...${expected.drop(column).take(40)}\n" +
                "  actual:   ...${actual.drop(column).take(40)}",
        )
    }

    /**
     * One BGZF member whose deflate payload is a stored block. Enough to build a HapMap the production
     * code has no compressor for.
     */
    private fun bgzfMember(
        data: ByteArray,
        name: String? = null,
        comment: String? = null,
        headerCrc: Boolean = false,
    ): ByteArray {
        val stored = ByteArray(5 + data.size)
        stored[0] = 1
        stored[1] = (data.size and 0xFF).toByte()
        stored[2] = ((data.size ushr 8) and 0xFF).toByte()
        val complement = data.size xor 0xFFFF
        stored[3] = (complement and 0xFF).toByte()
        stored[4] = ((complement ushr 8) and 0xFF).toByte()
        data.copyInto(stored, 5)

        val nameBytes = name?.encodeToByteArray()?.plus(0.toByte()) ?: ByteArray(0)
        val commentBytes = comment?.encodeToByteArray()?.plus(0.toByte()) ?: ByteArray(0)
        val headerCrcBytes = if (headerCrc) byteArrayOf(0, 0) else ByteArray(0)
        val total = 18 + nameBytes.size + commentBytes.size + headerCrcBytes.size + stored.size + 8
        val block = ByteArray(total)
        block[0] = 0x1F
        block[1] = 0x8B.toByte()
        block[2] = 8
        var flags = 4
        if (name != null) flags = flags or 8
        if (comment != null) flags = flags or 16
        if (headerCrc) flags = flags or 2
        block[3] = flags.toByte()
        block[9] = 0xFF.toByte()
        block[10] = 6
        block[12] = 'B'.code.toByte()
        block[13] = 'C'.code.toByte()
        block[14] = 2
        val bsize = total - 1
        block[16] = (bsize and 0xFF).toByte()
        block[17] = ((bsize ushr 8) and 0xFF).toByte()
        var cursor = 18
        nameBytes.copyInto(block, cursor)
        cursor += nameBytes.size
        commentBytes.copyInto(block, cursor)
        cursor += commentBytes.size
        headerCrcBytes.copyInto(block, cursor)
        cursor += headerCrcBytes.size
        stored.copyInto(block, cursor)
        cursor += stored.size
        val checksum = Crc32().also { it.update(data, 0, data.size) }.value
        putU32(block, cursor, checksum)
        putU32(block, cursor + 4, data.size)
        return block
    }

    private fun u32(file: ByteArray, at: Int): Int {
        fun u16(index: Int) = (file[index].toInt() and 0xFF) or ((file[index + 1].toInt() and 0xFF) shl 8)
        return u16(at) or (u16(at + 2) shl 16)
    }

    private fun putU32(dest: ByteArray, at: Int, value: Int) {
        dest[at] = (value and 0xFF).toByte()
        dest[at + 1] = ((value ushr 8) and 0xFF).toByte()
        dest[at + 2] = ((value ushr 16) and 0xFF).toByte()
        dest[at + 3] = ((value ushr 24) and 0xFF).toByte()
    }

    private class CountingRawSource(private val data: ByteArray) : RawSource {
        var position: Int = 0
            private set

        override fun readAtMostTo(sink: Buffer, byteCount: Long): Long {
            if (byteCount == 0L) return 0L
            if (position >= data.size) return -1L
            val count = minOf(byteCount, (data.size - position).toLong()).toInt()
            sink.write(data, position, position + count)
            position += count
            return count.toLong()
        }

        override fun close() = Unit
    }

    private companion object {
        const val TUTORIAL_GZ = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
        const val GOLDEN = "golden/vcf/mdp_genotype_w_depth.hmp.txt"
        const val SOURCE = "htslib/variants/source.vcf"
        const val MID_EOF = "htslib/variants/mid-eof.vcf.gz"
        const val NO_EOF = "htslib/variants/no-eof.vcf.gz"
        const val TRUNCATED = "htslib/variants/truncated.vcf.gz"
        const val BAD_CRC = "htslib/variants/bad-crc.vcf.gz"
        const val PLAIN = "htslib/variants/plain.vcf.gz"
        const val CONCATENATED = "htslib/variants/concatenated.vcf.gz"
    }
}
