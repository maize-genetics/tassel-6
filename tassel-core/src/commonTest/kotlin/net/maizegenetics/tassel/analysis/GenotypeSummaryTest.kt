package net.maizegenetics.tassel.analysis

import kotlinx.io.Buffer
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.table.ResultTable
import net.maizegenetics.tassel.taxa.TaxaList
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.Tolerance
import net.maizegenetics.tassel.testing.assertMatchesGolden
import net.maizegenetics.tassel.testing.goldenTolerance
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * The genotype summary of the tutorial HapMap file against the TASSEL 5 tables in `testdata/golden/summary/`
 * (ADR-0009). Each table is checked twice: against the rounded `.txt` export, which pins the layout, the placeholder
 * text, and TASSEL's five-significant-digit formatting, and against the `.full.txt` export, which pins the values.
 *
 * @see net.maizegenetics.analysis.data.GenotypeSummaryPlugin
 */
class GenotypeSummaryTest {
    @Test
    fun overviewMatchesGolden() = assertMatchesGoldens("overview", tutorialSummary.overview)

    @Test
    fun allelesMatchGolden() = assertMatchesGoldens("alleles", tutorialSummary.alleles)

    @Test
    fun sitesMatchGolden() = assertMatchesGoldens("sites", tutorialSummary.sites)

    @Test
    fun taxaMatchGolden() = assertMatchesGoldens("taxa", tutorialSummary.taxa)

    @Test
    fun tablesAreShapedAsTassel5() {
        assertEquals(
            listOf("Overall Summary", "Allele Summary", "Site Summary", "Taxa Summary"),
            tutorialSummary.tables.map { it.title },
        )
        assertEquals(15, tutorialSummary.overview!!.rowCount)
        assertEquals(3093, tutorialSummary.sites!!.rowCount)
        assertEquals(281, tutorialSummary.taxa!!.rowCount)
        // Fifteen fixed columns, four each for alleles 3 to 6, then six more.
        assertEquals(37, tutorialSummary.sites!!.columnCount)
    }

    @Test
    fun countsAreDoublesInTheOverviewAndLongsInTheAlleleTable() {
        // TASSEL 5 holds every overview value as a Double, which is why its full-precision export says "281.0".
        assertTrue(tutorialSummary.overview!!.rows.all { it[1] is Double }, "overview values")
        // The allele table holds its counts as Long, so they print without a decimal point.
        assertTrue(tutorialSummary.alleles!!.rows.all { it[1] is Long }, "allele counts")
    }

    @Test
    fun majorMinorPairRowsHaveNoFrequency() {
        val pairs = tutorialSummary.alleles!!.rows.filter { (it[0] as String).contains(':') }
        assertEquals(16, pairs.size)
        pairs.forEach { assertNull(it[3], "frequency of ${it[0]}") }
        // Every site is counted once, by its major and minor allele.
        assertEquals(3093L, pairs.sumOf { it[1] as Long })
    }

    @Test
    fun eachReportCanBeTurnedOff() {
        val onlySites = genotypeSummary(tutorial, overview = false, taxaSummary = false)
        assertNull(onlySites.overview)
        assertNull(onlySites.alleles)
        assertNull(onlySites.taxa)
        assertEquals(listOf("Site Summary"), onlySites.tables.map { it.title })
        val withoutSites = genotypeSummary(tutorial, siteSummary = false, taxaSummary = false)
        assertEquals(listOf("Overall Summary", "Allele Summary"), withoutSites.tables.map { it.title })
        // TASSEL 5 recomputes the overview totals in its own pass when the site summary is off. The tutorial file has
        // no rare or undefined allele codes, so that pass gives the same numbers as the site pass.
        assertEquals(tutorialSummary.overview!!.rows, withoutSites.overview!!.rows)
        val nothing = genotypeSummary(tutorial, overview = false, siteSummary = false, taxaSummary = false)
        assertTrue(nothing.tables.isEmpty())
    }

    @Test
    fun countsOneSidedMissingGenotypesAsHeterozygous() {
        // TASSEL 5's GenotypeTable.heterozygousCount(site) asks whether the two allele codes differ, so A/N counts as
        // heterozygous even though AlleleStats.heterozygousCount does not count it.
        val summary = genotypeSummary(oneSite("AA", "AN", "AC", "NN"))
        val sites = summary.sites!!
        assertEquals(2, sites[0, sites.columnIndex("Number Heterozygous")])
        assertEquals(2.0 / 3.0, sites[0, sites.columnIndex("Proportion Heterozygous")])
        // A/N calls one gamete and N/N calls neither, so three of the eight gametes are missing.
        assertEquals(3, sites[0, sites.columnIndex("Gametes Missing")])
        val taxa = summary.taxa!!
        val heterozygous = taxa.columnIndex("Number Heterozygous")
        assertEquals(listOf(0, 1, 1, 0), (0 until 4).map { taxa[it, heterozygous] })
        val overview = summary.overview!!
        assertEquals(listOf("Number Heterozygous", 2.0), overview.rows[12])
        assertEquals(listOf("Number Missing", 1.0), overview.rows[5])
    }

    @Test
    fun reportsAbsentAllelesAsTassel5Placeholders() {
        val summary = genotypeSummary(oneSite("AA", "AA"))
        val sites = summary.sites!!
        assertEquals("NA", sites[0, sites.columnIndex("Minor Allele")])
        assertEquals(0, sites[0, sites.columnIndex("Minor Allele Gametes")])
        assertEquals(0.0, sites[0, sites.columnIndex("Minor Allele Proportion")])
        assertEquals(0.0, sites[0, sites.columnIndex("Minor Allele Frequency")])
        assertEquals("TBD", sites[0, sites.columnIndex("Inbreeding Coefficient")])
        assertEquals("TBD", sites[0, sites.columnIndex("Inbreeding Coefficient Scaled by Missing")])
        val taxa = summary.taxa!!
        assertEquals("Inbreeding Coefficient", taxa[0, taxa.columnIndex("Inbreeding Coefficient")])
        assertEquals("ICSBM", taxa[0, taxa.columnIndex("Inbreeding Coefficient Scaled by Missing")])
        // A site with one allele is counted under "A:A", as TASSEL 5's majorMinorCounts counts it.
        assertEquals(listOf("A:A", 1L, 1.0, null), summary.alleles!!.rows[1])
    }

    /**
     * Compares one table with both of its goldens: the rounded export as written, and the full-precision export
     * against a rendering that prints each `Double` with all its digits. The manifest tolerances (`exact` and
     * `abs 1e-12`) come from `testdata/golden/MANIFEST`.
     */
    private fun assertMatchesGoldens(name: String, table: ResultTable?) {
        assertEquals(Tolerance.Exact, goldenTolerance("golden/summary/$name.txt"))
        assertEquals(Tolerance.Absolute(1e-12), goldenTolerance("golden/summary/$name.full.txt"))
        assertMatchesGolden("golden/summary/$name.txt", table!!.toDelimitedString())
        assertMatchesGolden("golden/summary/$name.full.txt", table.toFullPrecisionString())
    }

    /** The table with every `Double` written in full, as `scripts/golden/DoubleFormat.java` makes TASSEL 5 write it. */
    private fun ResultTable.toFullPrecisionString(): String = buildString {
        columnNames.joinTo(this, "\t")
        append('\n')
        rows.forEach { cells ->
            cells.joinTo(this, "\t") { it?.toString() ?: "" }
            append('\n')
        }
    }

    /** One site over [calls], so the hand-checked cases need no file. */
    private fun oneSite(vararg calls: String): GenotypeTable = InMemoryGenotypeTable.of(
        TaxaList.of(List(calls.size) { "T$it" }),
        PositionList(listOf(Position(Chromosome.of("1"), 100, "s1"))),
        ByteArray(calls.size) { Nucleotide.genotype(calls[it]) },
    )

    private companion object {
        val tutorial: GenotypeTable by lazy {
            HapMap.read(Buffer().also { it.write(TestData.readBytes("tutorial/mdp_genotype.hmp.txt")) })
        }

        /** Summarizing 281 taxa x 3,093 sites is slow on Node and Native, so every test shares one run. */
        val tutorialSummary: GenotypeSummary by lazy { genotypeSummary(tutorial) }
    }
}
