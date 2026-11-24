package net.maizegenetics.tassel.io.index

import kotlinx.io.Buffer
import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * BCF CSI files from #118 carry no contig names (`l_aux` is 0). The names passed here are the header
 * order recorded by that issue: tutorial chromosomes 1 through 10, the wheat `##contig` order in
 * `scripts/generate-htslib-fixtures.sh`, and `chr1` for the boundary record. The test does not read a BCF.
 *
 * Region notes are the `.record.txt` files #118 checked as byte-identical to `bcftools view -H -r` on the
 * BCF. A non-empty note must produce chunks. An empty note on a contig that has other records need not:
 * a CSI bin can overlap a range whose records do not.
 */
class BcfCsiIndexTest {
    @Test
    fun tutorialContigsMatchTheBcfIndexStats() {
        val index = open(TUTORIAL_CSI, TUTORIAL_NAMES)
        assertEquals(14, index.minShift)
        assertEquals(6, index.depth)
        assertEquals(1L shl 32, index.maximumCoordinate)
        assertEquals(TUTORIAL_NAMES, index.contigs.map { it.name })
        assertMapped(TUTORIAL_COUNTS, index)
    }

    @Test
    fun wheatContigsMatchTheHeaderIncludingEmptyOnes() {
        val index = open(WHEAT_CSI, WHEAT_HEADER)
        assertEquals(14, index.minShift)
        assertEquals(6, index.depth)
        assertEquals(1L shl 32, index.maximumCoordinate)
        assertEquals(WHEAT_HEADER, index.contigs.map { it.name })
        val stats = counts(WHEAT_COUNTS)
        var previous = -1
        for ((name, mapped) in stats) {
            val slot = index.contigIndex(name)
            assertTrue(slot > previous, name)
            previous = slot
            assertEquals(mapped, index.contigs[slot].mapped, name)
        }
        val named = stats.map { it.first }.toSet()
        val absent = WHEAT_HEADER.filter { it !in named }
        assertEquals(16, absent.size)
        assertTrue("chr6D" in absent)
        for (name in absent) {
            assertEquals(0L, index.contigs[index.contigIndex(name)].mapped, name)
        }
    }

    @Test
    fun boundaryRecordSitsInItsOwnBin() {
        val index = open(BOUNDARY_CSI, listOf("chr1"))
        assertEquals(14, index.minShift)
        assertEquals(6, index.depth)
        assertEquals("chr1", index.contigs.single().name)
        assertEquals(1L, index.contigs.single().mapped)
        val slot = index.contigIndex("chr1")
        assertTrue(index.chunks(slot, 2_147_483_646, 2_147_483_647).isNotEmpty())
        assertEquals(emptyList(), index.chunks(slot, 0, 1))
    }

    @Test
    fun tutorialRegionsFollowTheHtslibNotes() {
        val index = open(TUTORIAL_CSI, TUTORIAL_NAMES)
        for (stem in TUTORIAL_REGIONS) assertRegion(index, "htslib/tutorial/regions", stem)
    }

    @Test
    fun wheatRegionsFollowTheHtslibNotes() {
        val index = open(WHEAT_CSI, WHEAT_HEADER)
        for (stem in WHEAT_REGIONS) assertRegion(index, "htslib/wheat/regions", stem)
    }

    @Test
    fun aNamelessCsiStillNeedsNames() {
        for (path in listOf(TUTORIAL_CSI, WHEAT_CSI, BOUNDARY_CSI)) {
            val failure = assertFailsWith<IllegalArgumentException> {
                readGenomicIndex(Buffer().also { it.write(TestData.readBytes(path)) })
            }
            assertContains(failure.message ?: "", "no contig names")
        }
    }

    @Test
    fun anIndexThatAlreadyHasNamesRejectsTheSuppliedList() {
        val csi = assertFailsWith<IllegalArgumentException> {
            readGenomicIndex(Buffer().also { it.write(TestData.readBytes("$VCF_TUTORIAL.csi")) }, listOf("1"))
        }
        assertEquals("CSI index already contains contig names", csi.message)
        val tbi = assertFailsWith<IllegalArgumentException> {
            readGenomicIndex(Buffer().also { it.write(TestData.readBytes("$VCF_TUTORIAL.tbi")) }, listOf("1"))
        }
        assertEquals("TBI index already contains contig names", tbi.message)
    }

    @Test
    fun aCountMismatchNamesTheShortSide() {
        val shortNames = assertFailsWith<IllegalArgumentException> {
            open(WHEAT_CSI, WHEAT_HEADER.dropLast(1))
        }
        assertEquals("contig name list is short: 20 names for 21 references", shortNames.message)
        val shortIndex = assertFailsWith<IllegalArgumentException> {
            open(WHEAT_CSI, WHEAT_HEADER + "chrUn")
        }
        assertEquals("genomic index is short: 21 references for 22 contig names", shortIndex.message)
        val duplicate = assertFailsWith<IllegalArgumentException> {
            open(WHEAT_CSI, WHEAT_HEADER.dropLast(1) + "chr7B")
        }
        assertContains(duplicate.message ?: "", "repeats contig")
    }

    private fun assertRegion(index: GenomicIndex, directory: String, stem: String) {
        val records = TestData.readBytes("$directory/$stem.record.txt")
        val parts = stem.split('_')
        val beg = parts.getOrNull(parts.size - 2)?.toLongOrNull()
        val end = parts.last().toLongOrNull()
        if (beg != null && end != null && parts.size >= 3) {
            val contig = parts.dropLast(2).joinToString("_")
            val slot = index.contigIndex(contig)
            if (contig == "99" || contig == "chrUn") {
                assertEquals(-1, slot, stem)
                assertEquals(0, records.size, stem)
                return
            }
            assertTrue(slot >= 0, stem)
            if (end < beg) {
                assertEquals(0, records.size, stem)
                val failure = assertFailsWith<IllegalArgumentException> { index.chunks(slot, beg - 1, end) }
                assertContains(failure.message ?: "", "before start")
                return
            }
            val chunks = index.chunks(slot, beg - 1, end)
            if (contig == "chr6D") {
                assertEquals(emptyList(), chunks, stem)
                assertEquals(0, records.size, stem)
                return
            }
            if (records.isNotEmpty()) assertTrue(chunks.isNotEmpty(), stem)
            return
        }
        val slot = index.contigIndex(stem)
        assertTrue(slot >= 0, stem)
        if (records.isNotEmpty()) {
            assertTrue(index.chunks(slot, 0, index.maximumCoordinate).isNotEmpty(), stem)
        }
    }

    private fun assertMapped(path: String, index: GenomicIndex) {
        for ((name, mapped) in counts(path)) {
            assertEquals(mapped, index.contigs[index.contigIndex(name)].mapped, name)
        }
    }

    private fun counts(path: String): List<Pair<String, Long>> =
        TestData.readText(path).trimEnd().lines().filter { it.isNotEmpty() }.map { line ->
            val columns = line.split('\t')
            columns[0] to columns[2].toLong()
        }

    private fun open(path: String, names: List<String>): GenomicIndex =
        readGenomicIndex(Buffer().also { it.write(TestData.readBytes(path)) }, names)

    private companion object {
        const val TUTORIAL_CSI = "htslib/bcf/tutorial/mdp_genotype_w_depth.bcf.csi"
        const val WHEAT_CSI = "htslib/bcf/wheat/wheat.bcf.csi"
        const val BOUNDARY_CSI = "htslib/bcf/boundary/pos2147483647.bcf.csi"
        const val TUTORIAL_COUNTS = "htslib/bcf/tutorial/index-s.txt"
        const val WHEAT_COUNTS = "htslib/bcf/wheat/index-s.txt"
        const val VCF_TUTORIAL = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"

        val TUTORIAL_NAMES = (1..10).map { it.toString() }
        val WHEAT_HEADER = listOf(
            "chr1A", "chr1B", "chr1D",
            "chr2A", "chr2B", "chr2D",
            "chr3A", "chr3B", "chr3D",
            "chr4A", "chr4B", "chr4D",
            "chr5A", "chr5B", "chr5D",
            "chr6A", "chr6B", "chr6D",
            "chr7A", "chr7B", "chr7D",
        )
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
    }
}
