package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.RawSource
import kotlinx.io.buffered
import kotlinx.io.readByteArray
import kotlinx.io.readString
import kotlinx.io.writeString
import net.maizegenetics.tassel.genotype.GenotypeTable
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
    fun readsPlainAndConcatenatedGzip() = runTest {
        val plainVcf = TestData.readBytes(SOURCE)
        val expected = Vcf.read(Buffer().also { it.write(plainVcf) })
        for (path in listOf(PLAIN, CONCATENATED)) {
            val opened = open(path)
            assertContentEquals(plainVcf, opened.source.readByteArray())
            assertEquals(emptyList(), opened.warnings, path)
            assertSameTable(expected, Vcf.read(open(path).source))
        }
    }

    @Test
    fun readsPlainGzipOneByteAtATime() {
        val compressed = TestData.readBytes(PLAIN)
        val opened = OneByteRawSource(compressed).buffered().decompressing()
        assertContentEquals(TestData.readBytes(SOURCE), opened.source.readByteArray())
        assertEquals(emptyList(), opened.warnings)
    }

    @Test
    fun readsPlainHeaderFieldsAndChecksTheHeaderCrc() {
        val text = "hello".encodeToByteArray()
        val gz = plainMember(text, name = "sample", comment = "note", extra = byteArrayOf(1, 2), headerCrc = true)
        val (bytes, warnings) = inflated(gz)
        assertContentEquals(text, bytes)
        assertEquals(emptyList(), warnings)
    }

    @Test
    fun readsPlainMembersOnEitherSideOfBgzf() {
        val text = "abcdefghij".encodeToByteArray()
        val forward = bgzfMember(text.copyOfRange(0, 3)) +
            plainMember(text.copyOfRange(3, 7)) +
            bgzfMember(text.copyOfRange(7, text.size)) +
            bgzfMember(ByteArray(0))
        val backward = plainMember(text.copyOfRange(0, 4)) +
            bgzfMember(text.copyOfRange(4, 8)) +
            plainMember(text.copyOfRange(8, text.size))
        val (forwardBytes, forwardWarnings) = inflated(forward)
        val (backwardBytes, backwardWarnings) = inflated(backward)
        assertContentEquals(text, forwardBytes)
        assertContentEquals(text, backwardBytes)
        assertEquals(emptyList(), forwardWarnings)
        assertEquals(emptyList(), backwardWarnings)
    }

    @Test
    fun rejectsAPlainHeaderCrcTrailerCrcAndAShortMember() {
        val text = "hello".encodeToByteArray()
        assertPlainFailure(plainMember(text, headerCrc = true, corruptHeaderCrc = true), "header CRC")
        assertPlainFailure(plainMember(text, trailerCrc = 1), "CRC32")
        assertPlainFailure(plainMember(text, trailerSize = text.size + 1), "ISIZE")
        assertPlainFailure(plainMember(text).copyOf(12), "truncated")
        val reserved = plainMember(text).also { it[3] = (it[3].toInt() or 0x20).toByte() }
        assertPlainFailure(reserved, "reserved")
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

    private fun assertSameTable(expected: GenotypeTable, actual: GenotypeTable) {
        assertEquals(expected.taxa.names, actual.taxa.names)
        assertEquals(expected.numberOfTaxa, actual.numberOfTaxa)
        assertEquals(expected.numberOfSites, actual.numberOfSites)
        for (site in 0 until expected.numberOfSites) {
            assertContentEquals(expected.siteGenotypes(site), actual.siteGenotypes(site), "site $site")
        }
    }

    private fun assertPlainFailure(bytes: ByteArray, snippet: String) {
        val failure = assertFailsWith<IllegalArgumentException> { inflated(bytes) }
        assertContains(failure.message ?: "", snippet)
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

    private fun plainMember(
        data: ByteArray,
        name: String? = null,
        comment: String? = null,
        extra: ByteArray? = null,
        headerCrc: Boolean = false,
        corruptHeaderCrc: Boolean = false,
        trailerCrc: Int? = null,
        trailerSize: Int? = null,
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
        val extraBytes = extra ?: ByteArray(0)
        val headerBody = ByteArray(
            10 + (if (extra != null) 2 + extraBytes.size else 0) + nameBytes.size + commentBytes.size,
        )
        headerBody[0] = 0x1F
        headerBody[1] = 0x8B.toByte()
        headerBody[2] = 8
        var flags = 0
        if (extra != null) flags = flags or 4
        if (name != null) flags = flags or 8
        if (comment != null) flags = flags or 16
        if (headerCrc) flags = flags or 2
        headerBody[3] = flags.toByte()
        var cursor = 10
        if (extra != null) {
            headerBody[cursor] = (extraBytes.size and 0xFF).toByte()
            headerBody[cursor + 1] = ((extraBytes.size ushr 8) and 0xFF).toByte()
            extraBytes.copyInto(headerBody, cursor + 2)
            cursor += 2 + extraBytes.size
        }
        nameBytes.copyInto(headerBody, cursor)
        cursor += nameBytes.size
        commentBytes.copyInto(headerBody, cursor)

        val crc16 = Crc32().also { it.update(headerBody, 0, headerBody.size) }.value and 0xFFFF
        val storedCrc = if (corruptHeaderCrc) crc16 xor 0xFF else crc16
        val headerCrcBytes = if (headerCrc) {
            byteArrayOf((storedCrc and 0xFF).toByte(), ((storedCrc ushr 8) and 0xFF).toByte())
        } else {
            ByteArray(0)
        }
        val checksum = trailerCrc ?: Crc32().also { it.update(data, 0, data.size) }.value
        val size = trailerSize ?: data.size
        val block = ByteArray(headerBody.size + headerCrcBytes.size + stored.size + 8)
        headerBody.copyInto(block)
        headerCrcBytes.copyInto(block, headerBody.size)
        stored.copyInto(block, headerBody.size + headerCrcBytes.size)
        val trailerAt = block.size - 8
        putU32(block, trailerAt, checksum)
        putU32(block, trailerAt + 4, size)
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

    private class OneByteRawSource(private val data: ByteArray) : RawSource {
        private var position = 0

        override fun readAtMostTo(sink: Buffer, byteCount: Long): Long {
            if (byteCount == 0L) return 0L
            if (position >= data.size) return -1L
            sink.write(byteArrayOf(data[position]))
            position++
            return 1L
        }

        override fun close() = Unit
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
