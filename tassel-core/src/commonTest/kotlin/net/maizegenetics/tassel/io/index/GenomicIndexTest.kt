package net.maizegenetics.tassel.io.index

import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.io.decompressing
import net.maizegenetics.tassel.io.deflate.bgzfBlocks
import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * CSI and TBI queries against the htslib 1.25 fixtures (ADR-0021).
 *
 * Record virtual offsets come from a scan of the BGZF blocks. The seekable reader is #87.
 */
class GenomicIndexTest {
    @Test
    fun parsesEveryFixtureIndex() {
        val tutorialCsi = open("$TUTORIAL.csi")
        val tutorialTbi = open("$TUTORIAL.tbi")
        val wheat = open("$WHEAT.csi")
        val large = open("$LARGE.csi")

        assertEquals(14, tutorialCsi.minShift)
        assertEquals(8, tutorialCsi.depth)
        assertEquals(1L shl 38, tutorialCsi.maximumCoordinate)
        assertEquals(14, tutorialTbi.minShift)
        assertEquals(5, tutorialTbi.depth)
        assertEquals(1L shl 29, tutorialTbi.maximumCoordinate)
        assertEquals(14, wheat.minShift)
        assertEquals(6, wheat.depth)
        assertEquals(1L shl 32, wheat.maximumCoordinate)
        assertEquals(14, large.minShift)
        assertEquals(6, large.depth)

        assertEquals(TestData.readText(TUTORIAL_NAMES).trimEnd().lines(), tutorialCsi.contigs.map { it.name })
        assertEquals(tutorialCsi.contigs.map { it.name }, tutorialTbi.contigs.map { it.name })
        assertCounts(TUTORIAL_COUNTS, tutorialCsi)
        assertCounts(TUTORIAL_COUNTS, tutorialTbi)
        assertCounts(WHEAT_COUNTS, wheat)
        assertCounts(LARGE_COUNTS, large)
        assertEquals(wheat.contigs.map { it.name }, counts(WHEAT_COUNTS).map { it.first })
    }

    @Test
    fun reg2binAgreesWithTheCsiSpecification() {
        assertEquals(4681L, reg2bin(0, 1, 14, 5))
        assertEquals(4681L, reg2bin(16383, 16384, 14, 5))
        assertEquals(4682L, reg2bin(16384, 16385, 14, 5))
        assertEquals(37448L, reg2bin((1L shl 29) - 1, 1L shl 29, 14, 5))
        assertEquals(listOf(0L, 1L, 9L, 73L, 585L, 4681L), reg2bins(0, 1, 14, 5))
        assertEquals(listOf(0L, 8L, 72L, 584L, 4680L, 37448L), reg2bins((1L shl 29) - 1, 1L shl 29, 14, 5))
        val full = reg2bins(0, 1L shl 29, 14, 5)
        assertEquals(37449, full.size)
        assertEquals(0L, full.first())
        assertEquals(37448L, full.last())

        assertEquals(37449L, reg2bin(0, 1, 14, 6))
        assertEquals(70217L, reg2bin(1L shl 29, (1L shl 29) + 1, 14, 6))
        assertEquals(74070L, reg2bin(600_000_000, 600_000_001, 14, 6))
        assertEquals(listOf(0L, 2L, 17L, 144L, 1157L, 9258L, 74070L), reg2bins(600_000_000, 600_000_001, 14, 6))

        assertEquals(2396745L, reg2bin(0, 1, 14, 8))
        assertEquals(19173960L, reg2bin((1L shl 38) - 2, (1L shl 38) - 1, 14, 8))
        assertEquals(
            listOf(0L, 1L, 9L, 73L, 585L, 4681L, 37449L, 299593L, 2396745L),
            reg2bins(0, 1, 14, 8),
        )
        assertEquals(emptyList(), reg2bins(10, 10, 14, 5))
    }

    @Test
    fun queriesCoverEveryOracleRecord() {
        val tutorialLines = lineOffsets(TestData.readBytes(TUTORIAL))
        val tutorialCsi = open("$TUTORIAL.csi")
        val tutorialTbi = open("$TUTORIAL.tbi")
        for (region in regions("htslib/tutorial/regions")) {
            val coveredByCsi = coveredLines(region, tutorialCsi, tutorialLines)
            val coveredByTbi = coveredLines(region, tutorialTbi, tutorialLines)
            assertEquals(coveredByCsi, coveredByTbi, region.stem)
        }
        cover("htslib/tutorial/regions", tutorialCsi, tutorialLines)
        cover("htslib/wheat/regions", open("$WHEAT.csi"), lineOffsets(TestData.readBytes(WHEAT)))
        cover("htslib/large-pos/regions", open("$LARGE.csi"), lineOffsets(TestData.readBytes(LARGE)))
    }

    @Test
    fun smallRegionsReadASmallPartOfTheFile() {
        val tutorial = TestData.readBytes(TUTORIAL)
        val wheat = TestData.readBytes(WHEAT)
        val tutorialIndex = open("$TUTORIAL.csi")
        val wheatIndex = open("$WHEAT.csi")
        val reports = listOf(
            report("tutorial 1:157104", tutorialIndex, "1", 157103, 157104, tutorial.size),
            report("tutorial 1:157104-10068726", tutorialIndex, "1", 157103, 10_068_726, tutorial.size),
            report("wheat chr3B:1000", wheatIndex, "chr3B", 999, 1000, wheat.size),
            report("wheat chr3B:180100-180200", wheatIndex, "chr3B", 180099, 180200, wheat.size),
            report("wheat chr3B:600010000-600020000", wheatIndex, "chr3B", 600_009_999, 600_020_000, wheat.size),
        )
        reports.forEach { println(it) }
        assertTrue(blocksTouched(tutorialIndex, "1", 157103, 157104) <= 2)
        assertTrue(blocksTouched(wheatIndex, "chr3B", 999, 1000) <= 2)
        assertTrue(blocksTouched(wheatIndex, "chr3B", 180099, 180200) <= 2)
        assertTrue(blocksTouched(wheatIndex, "chr3B", 600_009_999, 600_020_000) <= 4)
        assertTrue(spanOf(tutorialIndex, "1", 157103, 10_068_726) < tutorial.size / 5)
    }

    @Test
    fun aRegionPastTheLimitOrAnEmptyContigHasNoChunks() {
        val tbi = open("$TUTORIAL.tbi")
        assertEquals(emptyList(), tbi.chunks(0, 1L shl 29, (1L shl 29) + 1))
        val inverted = assertFailsWith<IllegalArgumentException> { tbi.chunks(0, 100, 50) }
        assertContains(inverted.message ?: "", "before start")
    }

    @Test
    fun offsetsAboveTwoToTheThirtyOneRoundTrip() {
        val index = open("$TUTORIAL.csi")
        val chunks = index.chunks(index.contigIndex("1"), 299_170_077, 299_999_999)
        assertTrue(chunks.isNotEmpty())
        assertTrue(chunks.any { it.start > Int.MAX_VALUE.toLong() || it.end > Int.MAX_VALUE.toLong() })

        val synthetic = syntheticCsi(bin = 1L shl 31, start = 3_000_000_000L, end = 3_000_001_000L)
        val parsed = readGenomicIndex(Buffer().also { it.write(synthetic) })
        assertEquals(1L, parsed.contigs.single().mapped)
        val hit = parsed.chunks(0, 0, 1)
        assertEquals(1, hit.size)
        assertEquals(3_000_000_000L, hit[0].start)
        assertTrue(hit[0].start > Int.MAX_VALUE.toLong())
    }

    @Test
    fun aBadIndexFailsWithAMessage() {
        val magic = assertFailsWith<IllegalArgumentException> {
            readGenomicIndex(Buffer().also { it.write(byteArrayOf(1, 2, 3, 4)) })
        }
        assertContains(magic.message ?: "", "magic")
        val truncated = assertFailsWith<IllegalArgumentException> {
            readGenomicIndex(Buffer().also { it.write(TestData.readBytes("$TUTORIAL.csi").copyOf(12)) })
        }
        assertContains(truncated.message ?: "", "truncated")
        val format = assertFailsWith<IllegalArgumentException> {
            readGenomicIndex(Buffer().also { it.write(syntheticCsi(format = 99)) })
        }
        assertContains(format.message ?: "", "unknown tabix format")
        val nameless = assertFailsWith<IllegalArgumentException> {
            readGenomicIndex(Buffer().also { it.write(syntheticCsi(names = false)) })
        }
        assertContains(nameless.message ?: "", "no contig names")
    }

    private fun cover(directory: String, index: GenomicIndex, lines: Map<String, List<Long>>) {
        for (region in regions(directory)) {
            if (region.inverted || index.contigIndex(region.contig) < 0) {
                assertEquals(emptyList(), region.records, region.stem)
                continue
            }
            val offsets = coveredLines(region, index, lines)
            for (record in region.records) {
                assertTrue(siteKey(record) in offsets, "${region.stem} misses ${siteKey(record)}")
            }
        }
    }

    private fun coveredLines(region: Region, index: GenomicIndex, lines: Map<String, List<Long>>): Set<String> {
        val contig = index.contigIndex(region.contig)
        if (region.inverted) {
            assertTrue(contig >= 0, region.stem)
            assertFailsWith<IllegalArgumentException> { index.chunks(contig, region.beg, region.end) }
            return emptySet()
        }
        if (contig < 0) return emptySet()
        val end = if (region.whole) index.maximumCoordinate else region.end
        val chunks = index.chunks(contig, region.beg, end)
        val covered = HashSet<String>()
        for ((text, offsets) in lines) {
            if (offsets.any { offset -> chunks.any { chunk -> inside(offset, chunk) } }) covered.add(siteKey(text))
        }
        return covered
    }

    private fun siteKey(line: String): String {
        val first = line.indexOf('\t')
        val second = line.indexOf('\t', first + 1)
        return line.substring(0, second)
    }

    private fun inside(offset: Long, chunk: BgzfChunk): Boolean =
        !unsignedLess(offset, chunk.start) && unsignedLess(offset, chunk.end)

    private fun spanOf(index: GenomicIndex, contig: String, beg: Long, end: Long): Long {
        val chunks = index.chunks(index.contigIndex(contig), beg, end)
        if (chunks.isEmpty()) return 0
        val from = chunks.minOf { it.start ushr 16 }
        val to = chunks.maxOf { it.end ushr 16 }
        return (to - from).coerceAtLeast(1)
    }

    private fun blocksTouched(index: GenomicIndex, contig: String, beg: Long, end: Long): Long {
        val chunks = index.chunks(index.contigIndex(contig), beg, end)
        return chunks.sumOf { chunk ->
            val from = chunk.start ushr 16
            val to = chunk.end ushr 16
            if (to == from) 1 else to - from + if ((chunk.end and 0xFFFFL) == 0L) 0 else 1
        }
    }

    private fun report(label: String, index: GenomicIndex, contig: String, beg: Long, end: Long, file: Int): String {
        val chunks = index.chunks(index.contigIndex(contig), beg, end)
        val from = chunks.minOf { it.start ushr 16 }
        val to = chunks.maxOf { it.end ushr 16 }
        val where = if (from == to) "one BGZF block" else "${to - from} compressed bytes"
        return "selectivity $label chunks=${chunks.size} $where of $file"
    }

    private fun open(path: String): GenomicIndex =
        readGenomicIndex(Buffer().also { it.write(TestData.readBytes(path)) })

    private fun assertCounts(path: String, index: GenomicIndex) {
        for ((name, count) in counts(path)) {
            val contig = index.contigs[index.contigIndex(name)]
            assertEquals(count, contig.mapped, name)
        }
    }

    private fun counts(path: String): List<Pair<String, Long>> =
        TestData.readText(path).trimEnd().lines().filter { it.isNotEmpty() }.map { line ->
            val columns = line.split('\t')
            columns[0] to columns[2].toLong()
        }

    private fun regions(directory: String): List<Region> {
        val names = when (directory) {
            "htslib/tutorial/regions" -> TUTORIAL_REGIONS
            "htslib/wheat/regions" -> WHEAT_REGIONS
            else -> LARGE_REGIONS
        }
        return names.map { stem ->
            val text = TestData.readText("$directory/$stem.record.txt").trimEnd()
            val records = if (text.isEmpty()) emptyList() else text.lines()
            val parts = stem.split('_')
            val beg = parts.getOrNull(parts.size - 2)?.toLongOrNull()
            val end = parts.last().toLongOrNull()
            if (parts.size >= 3 && beg != null && end != null) {
                val contig = parts.dropLast(2).joinToString("_")
                Region(stem, contig, beg - 1, end, whole = false, inverted = end < beg, records = records)
            } else {
                Region(stem, stem, 0, 0, whole = true, inverted = false, records = records)
            }
        }
    }

    private fun lineOffsets(gz: ByteArray): Map<String, List<Long>> {
        val plain = Buffer().also { it.write(gz) }.decompressing().source.readByteArray()
        val blocks = bgzfBlocks(gz)
        var blockStart = 0
        val spans = blocks.map { block ->
            val span = blockStart to block.uncompressedSize
            blockStart = block.to + 8
            span
        }
        val lines = HashMap<String, MutableList<Long>>()
        var position = 0
        var block = 0
        var origin = 0
        while (position < plain.size) {
            var end = position
            while (end < plain.size && plain[end] != LINE_FEED) end++
            while (block < spans.size && position >= origin + spans[block].second) {
                origin += spans[block].second
                block++
            }
            if (block < spans.size && spans[block].second > 0 && plain[position] != HASH) {
                val virtual = (spans[block].first.toLong() shl 16) or (position - origin).toLong()
                val text = plain.decodeToString(position, end)
                lines.getOrPut(text) { ArrayList() }.add(virtual)
            }
            position = if (end < plain.size) end + 1 else end
        }
        return lines
    }

    /**
     * A one-reference CSI. [bin] `1L shl 31` checks that a bin above 2^31 survives as an unsigned value,
     * and the chunk offset is above 2^31 as well. The queried bin is the depth-5 leaf that holds position 0.
     */
    private fun syntheticCsi(
        bin: Long = 4681L,
        start: Long = 3_000_000_000L,
        end: Long = 3_000_001_000L,
        format: Int = 2,
        names: Boolean = true,
    ): ByteArray {
        val out = Writer()
        out.bytes(byteArrayOf('C'.code.toByte(), 'S'.code.toByte(), 'I'.code.toByte(), 1))
        out.i32(14)
        out.i32(5)
        if (!names) {
            out.i32(0)
            out.i32(1)
            out.i32(0)
            return out.toArray()
        }
        val name = "c".encodeToByteArray() + byteArrayOf(0)
        out.i32(28 + name.size)
        out.i32(format)
        out.i32(1)
        out.i32(2)
        out.i32(0)
        out.i32('#'.code)
        out.i32(0)
        out.i32(name.size)
        out.bytes(name)
        out.i32(1)
        out.i32(if (bin == 4681L) 2 else 3)
        out.u32(4681L)
        out.u64(0)
        out.i32(1)
        out.u64(start)
        out.u64(end)
        if (bin != 4681L) {
            out.u32(bin)
            out.u64(0)
            out.i32(1)
            out.u64(start + 10)
            out.u64(end + 10)
        }
        out.u32(37450L)
        out.u64(0)
        out.i32(2)
        out.u64(start)
        out.u64(end)
        out.u64(1)
        out.u64(0)
        out.u64(0)
        return out.toArray()
    }

    private class Region(
        val stem: String,
        val contig: String,
        val beg: Long,
        val end: Long,
        val whole: Boolean,
        val inverted: Boolean,
        val records: List<String>,
    )

    private class Writer {
        private val bytes = ArrayList<Byte>()

        fun i32(value: Int) {
            bytes.add((value and 0xFF).toByte())
            bytes.add(((value ushr 8) and 0xFF).toByte())
            bytes.add(((value ushr 16) and 0xFF).toByte())
            bytes.add(((value ushr 24) and 0xFF).toByte())
        }

        fun u32(value: Long) = i32(value.toInt())

        fun u64(value: Long) {
            var remaining = value
            repeat(8) {
                bytes.add((remaining and 0xFF).toByte())
                remaining = remaining ushr 8
            }
        }

        fun bytes(value: ByteArray) {
            value.forEach { bytes.add(it) }
        }

        fun toArray(): ByteArray = bytes.toByteArray()
    }

    private companion object {
        const val TUTORIAL = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
        const val WHEAT = "htslib/wheat/wheat.vcf.gz"
        const val LARGE = "htslib/large-pos/large_pos.vcf.gz"
        const val TUTORIAL_NAMES = "htslib/tutorial/contigs.txt"
        const val TUTORIAL_COUNTS = "htslib/tutorial/index-s.txt"
        const val WHEAT_COUNTS = "htslib/wheat/index-s.txt"
        const val LARGE_COUNTS = "htslib/large-pos/index-s.txt"
        const val LINE_FEED = '\n'.code.toByte()
        const val HASH = '#'.code.toByte()

        val TUTORIAL_REGIONS = listOf(
            "1_157104_157104",
            "1_157104_10068726",
            "1_299170078_299999999",
            "99_1_1000",
            "1_10068726_157104",
        )
        val WHEAT_REGIONS = listOf(
            "chr3B_1000_1000",
            "chr3B_1000_5000000",
            "chr3B",
            "chr3B_599000000_601000000",
            "chr3B_700000000_700000000",
            "chr3B_536870002_536870010",
            "chr3B_180100_180200",
            "chr3B_600010000_600020000",
            "chr3B_800000000_851934019",
            "chr6D_1_495380293",
            "chrUn_1_1000",
            "chr3B_5000_1000",
        )
        val LARGE_REGIONS = listOf("chr1_3000000000_3000000000")
    }
}
