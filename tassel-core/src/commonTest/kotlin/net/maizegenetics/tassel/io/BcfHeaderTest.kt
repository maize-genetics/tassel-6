package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * BCF 2.2 header dictionaries against the htslib 1.25 fixtures (ADR-0022).
 *
 * `bcftools view` prints the stored header with `,IDX=` removed. The expected dictionary is the first appearance
 * of each id in that text. The explicit-`IDX` cases are built here, because every #118 file numbers without gaps.
 */
class BcfHeaderTest {
    @Test
    fun everyHtslibHeaderMatchesTheViewText() = runTest {
        for ((bcf, view) in VIEWED) {
            val bytes = bcfHeaderBytes(bcf)
            val header = parse(bytes)
            val stripped = withoutIdx(headerLines(bytes))
            assertEquals(viewHeader(view), stripped, bcf)
            assertMatchesAppearance(header, stripped, bcf)
            assertKnown(bcf, header)
        }
        for (bcf in HEADER_ONLY) {
            val bytes = bcfHeaderBytes(bcf)
            val header = parse(bytes)
            assertMatchesAppearance(header, withoutIdx(headerLines(bytes)), bcf)
            assertKnown(bcf, header)
        }
    }

    @Test
    fun rawStoredAndCompressedSmallHeadersMatch() = runTest {
        val compressed = parse(bcfHeaderBytes(SMALL))
        assertSame(compressed, parse(bcfHeaderBytes(SMALL_RAW)))
        assertSame(compressed, parse(bcfHeaderBytes(SMALL_STORED)))
    }

    @Test
    fun badMagicAndBcf21UseTheHtslibMessages() = runTest {
        val bad = bcfHeaderBytes(BAD_MAGIC)
        assertEquals(0, bad[0].toInt() and 0xff)
        assertEquals('C'.code, bad[1].toInt() and 0xff)
        assertEquals('F'.code, bad[2].toInt() and 0xff)
        val badMagic = assertFailsWith<IllegalArgumentException> { parse(bad) }
        assertEquals("Invalid BCF2 magic string", badMagic.message)

        val old = bcfHeaderBytes(BCF21)
        assertEquals("BCF", old.decodeToString(0, 3))
        assertEquals(2, old[3].toInt() and 0xff)
        assertEquals(1, old[4].toInt() and 0xff)
        val version = assertFailsWith<IllegalArgumentException> { parse(old) }
        assertEquals("Invalid BCF2 magic string: only BCFv2.2 is supported", version.message)
    }

    @Test
    fun truncatedRecordStillParsesTheHeader() = runTest {
        assertSame(parse(bcfHeaderBytes(SMALL)), parse(bcfHeaderBytes(TRUNCATED_RECORD)))
    }

    @Test
    fun explicitIdxLeavesGapsAndKeepsPassAtZero() {
        val header = header(
            """##FILTER=<ID=PASS,Description="All filters passed",IDX=0>""",
            """##INFO=<ID=DP,Number=1,Type=Integer,Description="depth",IDX=2>""",
            """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",IDX=4>""",
            """##contig=<ID=chr1,length=10,IDX=0>""",
            """##contig=<ID=chr3,length=30,IDX=2>""",
        )
        assertEquals(listOf("PASS", null, "DP", null, "GT"), header.strings)
        assertEquals(4, header.gtId)
        assertEquals(listOf("chr1", null, "chr3"), header.contigs.map { it?.id })
        assertEquals(listOf(10L, null, 30L), header.contigs.map { it?.length })
        assertEquals(listOf("chr1", "chr3"), header.vcf.contigs.map { it.id })
    }

    @Test
    fun missingIdxNumbersInOrderOfAppearance() {
        val header = header(
            """##FILTER=<ID=PASS,Description="All filters passed">""",
            """##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">""",
            """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">""",
            """##contig=<ID=chr1,length=10>""",
            """##contig=<ID=chr2,length=20>""",
        )
        assertEquals(listOf("PASS", "DP", "GT"), header.strings)
        assertEquals(2, header.gtId)
        assertEquals(listOf("chr1", "chr2"), header.contigs.map { it?.id })
        assertEquals(listOf(10L, 20L), header.contigs.map { it?.length })
    }

    @Test
    fun passIsZeroWhenItIsNotTheFirstLine() {
        val header = header(
            """##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">""",
            """##FILTER=<ID=PASS,Description="All filters passed">""",
            """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">""",
            """##contig=<ID=chr1,length=10>""",
        )
        assertEquals(listOf("PASS", "DP", "GT"), header.strings)
        assertEquals(2, header.gtId)
    }

    @Test
    fun implicitIdDoesNotFillAGap() {
        val header = header(
            """##FILTER=<ID=PASS,Description="All filters passed",IDX=0>""",
            """##INFO=<ID=DP,Number=1,Type=Integer,Description="depth",IDX=2>""",
            """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">""",
            """##contig=<ID=chr1,length=10,IDX=0>""",
        )
        assertEquals(listOf("PASS", null, "DP", "GT"), header.strings)
        assertEquals(3, header.gtId)
    }

    @Test
    fun aHeaderWithoutPassDoesNotInventOne() {
        val header = header(
            """##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">""",
            """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">""",
            """##contig=<ID=chr1,length=10>""",
        )
        assertEquals(listOf("DP", "GT"), header.strings)
        assertEquals(1, header.gtId)
    }

    @Test
    fun conflictingIdxIsRejected() {
        val failure = assertFailsWith<IllegalArgumentException> {
            header(
                """##FILTER=<ID=PASS,Description="All filters passed",IDX=0>""",
                """##INFO=<ID=DP,Number=1,Type=Integer,Description="depth",IDX=2>""",
                """##INFO=<ID=AF,Number=1,Type=Float,Description="frequency",IDX=2>""",
                """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",IDX=1>""",
                """##contig=<ID=chr1,length=10,IDX=0>""",
            )
        }
        assertContains(failure.message ?: "", "Conflicting IDX=2")
    }

    @Test
    fun vcfHeaderStillRequiresATaxon() {
        val text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        val failure = assertFailsWith<IllegalArgumentException> {
            readVcfHeader(Buffer().also { it.write(text.encodeToByteArray()) })
        }
        assertContains(failure.message ?: "", "at least one taxon")
    }

    private fun assertMatchesAppearance(header: BcfHeader, lines: List<String>, path: String) {
        val (strings, contigs) = appearance(lines)
        assertEquals(strings, header.strings.map { checkNotNull(it) }, path)
        assertEquals(contigs, header.contigs.map { checkNotNull(it).id }, path)
        assertEquals(header.vcf.contigs.map { it.id }, contigs, path)
        val gt = strings.indexOf("GT")
        assertEquals(if (gt < 0) null else gt, header.gtId, path)
    }

    private fun assertKnown(path: String, header: BcfHeader) {
        when (path) {
            SMALL -> {
                assertEquals(listOf("PASS", "GT"), header.strings)
                assertEquals(1, header.gtId)
                assertEquals("chr1", checkNotNull(header.contigs.single()).id)
                assertEquals(1000L, header.contigs.single()?.length)
                assertEquals(listOf("s1"), header.vcf.taxaNames)
            }
            BOUNDARY -> {
                assertEquals(listOf("PASS", "GT"), header.strings)
                assertEquals(2147483647L, header.contigs.single()?.length)
            }
            EDGE -> {
                assertEquals(141, header.strings.size)
                assertEquals("PASS", header.strings[0])
                assertEquals("q10", header.strings[1])
                assertEquals(136, header.gtId)
                assertEquals("GT", header.strings[136])
            }
            NO_SAMPLE -> {
                assertEquals(listOf("PASS", "NS"), header.strings)
                assertNull(header.gtId)
                assertEquals(emptyList<String>(), header.vcf.taxaNames)
                assertEquals("chr1", checkNotNull(header.contigs.single()).id)
            }
            TUTORIAL -> {
                assertEquals(8, header.strings.size)
                assertEquals("DP", header.strings[3])
                assertEquals("PL", header.strings[5])
                assertEquals("AF", header.strings[7])
                assertEquals(1, header.gtId)
                assertEquals((1..10).map { it.toString() }, header.contigs.map { checkNotNull(it).id })
                assertEquals(281, header.vcf.taxaNames.size)
                assertEquals("33-16", header.vcf.taxaNames.first())
            }
            WHEAT -> {
                assertEquals(listOf("PASS", "END", "SVLEN", "SVTYPE", "GT"), header.strings)
                assertEquals(4, header.gtId)
                assertEquals(21, header.contigs.size)
                assertEquals("chr7D", header.contigs.last()?.id)
                assertTrue(header.contigs.any { it?.id == "chr6D" })
            }
        }
    }

    private fun assertSame(expected: BcfHeader, actual: BcfHeader) {
        assertEquals(expected.strings, actual.strings)
        assertEquals(expected.gtId, actual.gtId)
        assertEquals(expected.contigs, actual.contigs)
        assertEquals(expected.vcf.taxaNames, actual.vcf.taxaNames)
    }

    private fun header(vararg meta: String): BcfHeader {
        val text = buildString {
            append("##fileformat=VCFv4.3\n")
            for (line in meta) append(line).append('\n')
            append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n")
        }
        return parse(rawBcf(text))
    }

    private fun parse(bytes: ByteArray): BcfHeader = readBcfHeader(Buffer().also { it.write(bytes) })

    private suspend fun bcfHeaderBytes(path: String): ByteArray {
        val file = TestData.readBytes(path)
        if (file.size < 2 || file[0] != GZIP_ID0 || file[1] != GZIP_ID1) return file
        return ByteArrayRandomAccessSource(file).use { source ->
            BgzfReader(source).use { reader ->
                val prefix = readExact(reader, 9)
                val textLength = uint32Le(prefix, 5).toInt()
                prefix + readExact(reader, textLength)
            }
        }
    }

    private suspend fun readExact(reader: BgzfReader, count: Int): ByteArray {
        val out = ByteArray(count)
        var at = 0
        while (at < count) {
            val n = reader.read(out, at, count - at)
            check(n > 0) { "BCF header ended after $at of $count bytes" }
            at += n
        }
        return out
    }

    private companion object {
        const val SMALL = "htslib/bcf/small/small.bcf"
        const val SMALL_RAW = "htslib/bcf/small/small.raw.bcf"
        const val SMALL_STORED = "htslib/bcf/small/small.stored.bcf"
        const val BOUNDARY = "htslib/bcf/boundary/pos2147483647.bcf"
        const val EDGE = "htslib/bcf/edge/edge.bcf"
        const val NO_SAMPLE = "htslib/bcf/edge/no-sample.bcf"
        const val TUTORIAL = "htslib/bcf/tutorial/mdp_genotype_w_depth.bcf"
        const val WHEAT = "htslib/bcf/wheat/wheat.bcf"
        const val BAD_MAGIC = "htslib/bcf/corrupt/bad-magic.bcf"
        const val BCF21 = "htslib/bcf/corrupt/bcf21.bcf"
        const val TRUNCATED_RECORD = "htslib/bcf/corrupt/truncated-record.bcf"
        val GZIP_ID0 = 0x1f.toByte()
        val GZIP_ID1 = 0x8b.toByte()
        val IDX = Regex(",IDX=\\d+")
        val VIEWED = listOf(
            SMALL to "htslib/bcf/small/small.view.vcf",
            BOUNDARY to "htslib/bcf/boundary/pos2147483647.view.vcf",
            EDGE to "htslib/bcf/edge/edge.view.vcf",
            NO_SAMPLE to "htslib/bcf/edge/no-sample.view.vcf",
        )
        val HEADER_ONLY = listOf(TUTORIAL, WHEAT)

        fun rawBcf(text: String): ByteArray {
            val body = text.encodeToByteArray()
            val withNul = ByteArray(body.size + 1)
            body.copyInto(withNul)
            val out = ByteArray(9 + withNul.size)
            out[0] = 'B'.code.toByte()
            out[1] = 'C'.code.toByte()
            out[2] = 'F'.code.toByte()
            out[3] = 2
            out[4] = 2
            val length = withNul.size
            out[5] = (length and 0xff).toByte()
            out[6] = ((length ushr 8) and 0xff).toByte()
            out[7] = ((length ushr 16) and 0xff).toByte()
            out[8] = ((length ushr 24) and 0xff).toByte()
            withNul.copyInto(out, 9)
            return out
        }

        fun headerLines(bytes: ByteArray): List<String> {
            val length = uint32Le(bytes, 5).toInt()
            return linesOf(bytes.decodeToString(9, 9 + length - 1))
        }

        fun linesOf(text: String): List<String> {
            val lines = text.split('\n').map { it.trimEnd('\r') }
            return if (lines.isNotEmpty() && lines.last().isEmpty()) lines.dropLast(1) else lines
        }

        fun withoutIdx(lines: List<String>): List<String> = lines.map { line ->
            if (line.startsWith("##") && line.endsWith('>')) line.replace(IDX, "") else line
        }

        fun viewHeader(path: String): List<String> {
            val lines = ArrayList<String>()
            for (line in linesOf(TestData.readText(path))) {
                if (!line.startsWith("#")) break
                lines.add(line)
            }
            return lines
        }

        fun appearance(lines: List<String>): Pair<List<String>, List<String>> {
            val strings = ArrayList<String>()
            val seen = HashSet<String>()
            val contigs = ArrayList<String>()
            val seenContigs = HashSet<String>()
            for (line in lines) {
                val id = idOf(line) ?: continue
                when {
                    line.startsWith("##contig=<") -> if (seenContigs.add(id)) contigs.add(id)
                    line.startsWith("##FILTER=<") ||
                        line.startsWith("##INFO=<") ||
                        line.startsWith("##FORMAT=<") -> if (seen.add(id)) strings.add(id)
                }
            }
            return strings to contigs
        }

        fun idOf(line: String): String? {
            val at = line.indexOf("<ID=")
            if (at < 0) return null
            val start = at + 4
            val comma = line.indexOf(',', start)
            val end = if (comma < 0) line.indexOf('>', start) else comma
            if (end < 0) return null
            return line.substring(start, end).ifEmpty { null }
        }

        fun uint32Le(bytes: ByteArray, at: Int): Long = (bytes[at].toLong() and 0xff) or
            ((bytes[at + 1].toLong() and 0xff) shl 8) or
            ((bytes[at + 2].toLong() and 0xff) shl 16) or
            ((bytes[at + 3].toLong() and 0xff) shl 24)
    }
}
