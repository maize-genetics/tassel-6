package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readString
import kotlinx.io.writeString
import net.maizegenetics.tassel.filter.TaxaFilter
import net.maizegenetics.tassel.filter.filterTaxa
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import net.maizegenetics.tassel.testing.GoldenManifest
import net.maizegenetics.tassel.testing.SuspendLazy
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.Tolerance
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals

/**
 * TASSEL 5 VCF export, against `testdata/golden/vcf-export/` (ADR-0009). The inputs are HapMap, so the reference
 * allele is unknown and REF is the major allele.
 */
class VcfWriteTest {
    @Test
    fun writesTheTutorialExport() = runTest {
        assertEquals(Tolerance.Exact(), GoldenManifest.load().tolerance(TUTORIAL_GOLDEN))
        assertMatchesGolden(TUTORIAL_GOLDEN, write(tutorial()))
    }

    @Test
    fun writesTheIncludedTaxaExport() = runTest {
        assertEquals(Tolerance.Exact(), GoldenManifest.load().tolerance(INCLUDE_GOLDEN))
        val included = filterTaxa(
            tutorial(),
            TaxaFilter(taxaNames = listOf("W22", "B73", "CML247", "33-16"), includeTaxa = true),
        )
        assertMatchesGolden(INCLUDE_GOLDEN, write(included))
    }

    @Test
    fun writesTheEdgeExport() = runTest {
        assertEquals(Tolerance.Exact(), GoldenManifest.load().tolerance(EDGES_GOLDEN))
        val table = HapMap.read(Buffer().also { it.write(TestData.readBytes(EDGES)) })
        assertMatchesGolden(EDGES_GOLDEN, write(table))
    }

    @Test
    fun readingBackTheTutorialGivesTheSameTable() = runTest {
        val table = tutorial()
        val reread = Vcf.read(Buffer().also { it.writeString(write(table)) })
        assertEquals(table.taxa.names, reread.taxa.names)
        assertEquals(table.numberOfSites, reread.numberOfSites)
        for (site in 0 until table.numberOfSites) {
            assertEquals(table.positions[site], reread.positions[site], "site $site")
            assertContentEquals(table.siteGenotypes(site), reread.siteGenotypes(site), "site $site")
        }
    }

    @Test
    fun reportsProgressAtTheEndOfTheOnlyBlock() = runTest {
        // 3,093 sites fit in one block, so a block-boundary report fires once, after the write.
        val reports = mutableListOf<Progress>()
        Vcf.write(tutorial(), Buffer()) { reports += it }
        assertEquals(listOf(Progress(3093, 3093, "Writing VCF")), reports)
    }

    @Test
    fun writesNPrefixWhenTheSiteListsNoVariants() = runTest {
        // No alleles column, so the no-known-variant branch runs. A and the gap tie, and A wins by allele code.
        val table = InMemoryGenotypeTable.of(
            TaxaList.of("t1", "t2"),
            PositionList(listOf(Position(Chromosome.of("1"), 10, snpId = "gap"))),
            byteArrayOf(Nucleotide.genotype("AA"), Nucleotide.genotype("--")),
        )
        val record = write(table).lineSequence().single { it.startsWith("1\t") }
        assertEquals("1\t10\tgap\tNA\tN\t.\tPASS\t.\tGT\t0/0\t1/1", record)
    }

    private suspend fun write(table: GenotypeTable): String = Buffer().also { Vcf.write(table, it) }.readString()

    private fun assertMatchesGolden(goldenPath: String, written: String) {
        val golden = normalizeLineEndings(TestData.readText(goldenPath)).lines()
        val actual = normalizeLineEndings(written).lines()
        assertEquals(golden.size, actual.size, "number of lines in $goldenPath")
        actual.forEachIndexed { index, line ->
            if (line != golden[index]) reportLineDifference(goldenPath, index, golden[index], line)
        }
    }

    private fun reportLineDifference(source: String, line: Int, expected: String, actual: String): Nothing {
        val column = expected.commonPrefixWith(actual).length
        throw AssertionError(
            "$source line ${line + 1} differs at column ${column + 1}:\n" +
                "  expected: ...${expected.drop(column).take(40)}\n" +
                "  actual:   ...${actual.drop(column).take(40)}",
        )
    }

    private companion object {
        const val TUTORIAL = "tutorial/mdp_genotype.hmp.txt"
        const val TUTORIAL_GOLDEN = "golden/vcf-export/mdp_genotype.vcf"
        const val INCLUDE_GOLDEN = "golden/vcf-export/include-list.vcf"
        const val EDGES = "tutorial/vcf_export_edges.hmp.txt"
        const val EDGES_GOLDEN = "golden/vcf-export/vcf_export_edges.vcf"

        val tutorial = SuspendLazy {
            HapMap.read(Buffer().also { it.write(TestData.readBytes(TUTORIAL)) })
        }
    }
}
