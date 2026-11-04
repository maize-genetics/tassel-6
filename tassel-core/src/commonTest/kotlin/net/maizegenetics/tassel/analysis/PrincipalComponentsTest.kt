package net.maizegenetics.tassel.analysis

import kotlinx.io.Buffer
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.table.PhenotypeTable
import net.maizegenetics.tassel.table.ResultTable
import net.maizegenetics.tassel.table.doubleToString
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.assertDelimitedMatches
import net.maizegenetics.tassel.testing.assertMatchesGolden
import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.goldenTolerance
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * PCA on the tutorial HapMap, against TASSEL 5's `PrincipalComponentsPlugin` output in
 * `testdata/golden/pca/` (ADR-0009).
 *
 * The eigenvector and principal-component goldens are compared whole. The two eigenvalue goldens
 * carry the MANIFEST's `rows-unordered` keyword, because TASSEL 5 returns its eigenvalues in an
 * order no sorted decomposition can reproduce ([tassel5ReturnsEigenvaluesOutOfOrderAfterTheSixtyThirdRow]),
 * so [eigenvaluesMatchGoldenInAnyRowOrder] compares them at that tolerance, which for
 * `eigenvalues.txt` means cell-exact text with rows paired by value rather than by line. The two
 * columns a row's position determines are checked where TASSEL 5's order is already sorted
 * ([eigenvalueIndexAndCumulativeProportionMatchGoldenWhileTassel5IsSorted]), and the one eigenvalue
 * that is mathematically zero is left to [theZeroEigenvalueIsLostInTheGramMatrix].
 *
 * @see net.maizegenetics.analysis.data.PrincipalComponentsPlugin
 * @see net.maizegenetics.stats.PCA.PrinComp
 */
class PrincipalComponentsTest {
    @Test
    fun eigenvectorsMatchGoldenLayout() {
        assertMatchesGolden(EIGENVECTORS, eigenvectors.toDelimitedString())
    }

    @Test
    fun eigenvectorsMatchGoldenValuesUpToSign() {
        assertMatchesGolden(EIGENVECTORS_FULL, fullPrecision(eigenvectors))
    }

    @Test
    fun principalComponentsMatchGoldenLayout() {
        assertPhenotypeMatchesGolden(PCS, pcs.toDelimitedString())
    }

    @Test
    fun principalComponentsMatchGoldenValuesUpToSign() {
        assertPhenotypeMatchesGolden(PCS_FULL, pcs.toDelimitedString())
    }

    /**
     * The `eigenvalue` and `proportion of total` columns do not depend on where a row sits, so
     * those two columns are compared over 280 of the 281 rows under the `rows-unordered` tolerance
     * the MANIFEST records, which pairs rows by value rather than by line. For `eigenvalues.txt`
     * that tolerance is `exact`, so each cell compared must match TASSEL 5's rounded text exactly;
     * for the full-precision twin it is `rel 1e-9`.
     *
     * The row left out is the one [theZeroEigenvalueIsLostInTheGramMatrix] checks instead.
     */
    @Test
    fun eigenvaluesMatchGoldenInAnyRowOrder() {
        eachEigenvalueGolden { golden, actual ->
            assertDelimitedMatches(
                expected = columns(withoutTheZeroEigenvalue(goldenText(golden)), VALUE, PROPORTION),
                actual = columns(withoutTheZeroEigenvalue(actual), VALUE, PROPORTION),
                tolerance = goldenTolerance(golden),
                source = golden,
            )
        }
    }

    /**
     * The centered matrix has rank 280, so its smallest singular value is zero and both sides are
     * reporting round-off. The Gram matrix squares the condition number, so where TASSEL 5 squares
     * a singular value of 3.1e-15 to get 9.3e-30, `G`'s own noise floor is around `eps * trace`,
     * and no decomposition of `G` can tell 3.5e-15 from zero. Both are zero to every digit that
     * means anything; neither is zero in the sense `rel 1e-9` asks for, which is why this one value
     * is checked here rather than in [eigenvaluesMatchGoldenInAnyRowOrder]. Dropping it is the only
     * thing that test does beyond the MANIFEST tolerance, and this test is why.
     */
    @Test
    fun theZeroEigenvalueIsLostInTheGramMatrix() {
        val golden = dataRows(goldenText(EIGENVALUES_FULL)).map { it[VALUE].toDouble() }.min()
        val actual = eigenvalues[eigenvalues.rowCount - 1, VALUE] as Double
        assertTrue(golden < NOISE_FLOOR, "TASSEL 5's smallest eigenvalue is $golden")
        assertTrue(actual in 0.0..NOISE_FLOOR, "the smallest eigenvalue is $actual, outside 0 .. $NOISE_FLOOR")
    }

    /**
     * The `PC` index and the running `cumulative proportion` are exactly what a row's position
     * means, so [eigenvaluesMatchGoldenInAnyRowOrder] cannot compare them. They are an oracle over
     * the leading rows where TASSEL 5's own order is already descending, and nowhere else. The `PC`
     * index leads this projection, so the canonical order `rows-unordered` sorts into is the order
     * TASSEL 5 wrote.
     */
    @Test
    fun eigenvalueIndexAndCumulativeProportionMatchGoldenWhileTassel5IsSorted() {
        eachEigenvalueGolden { golden, actual ->
            val rows = sortedPrefix(goldenText(golden))
            assertDelimitedMatches(
                expected = firstRows(columns(goldenText(golden), PC, CUMULATIVE), rows),
                actual = firstRows(columns(actual, PC, CUMULATIVE), rows),
                tolerance = goldenTolerance(golden),
                source = "$golden, first $rows rows",
            )
        }
    }

    /**
     * TASSEL 5's eigenvalues come from EJML's `SvdImplicitQrDecompose`, which does not sort its
     * singular values, so from row 63 on the golden interleaves a descending run of large values
     * with an ascending run of small ones. The Gram-matrix route this port takes returns them
     * sorted, which is also what the plugin's own help text promises ("sorted high to low").
     *
     * This test is the record of why the two eigenvalue goldens carry `rows-unordered` in the
     * MANIFEST. If a regenerated golden ever comes out sorted, it fails, and the keyword and
     * [eigenvalueIndexAndCumulativeProportionMatchGoldenWhileTassel5IsSorted] can both go.
     */
    @Test
    fun tassel5ReturnsEigenvaluesOutOfOrderAfterTheSixtyThirdRow() {
        val golden = dataRows(goldenText(EIGENVALUES_FULL)).map { it[VALUE].toDouble() }
        assertEquals(281, golden.size)
        assertEquals(281, golden.distinct().size, "ties would make the sorted comparison ambiguous")
        assertEquals(63, sortedPrefix(goldenText(EIGENVALUES_FULL)), "leading rows that are the largest, in order")
        assertTrue(golden[63] < golden[64], "row 63 starts an ascending run of the small eigenvalues")
        assertTrue(golden[63] < golden[280 - 1], "row 63 is smaller than rows that follow it")
    }

    @Test
    fun eigenvaluesCoverEveryTaxonAndSumToOne() {
        assertEquals(tutorial.numberOfTaxa, eigenvalues.rowCount)
        val columns = listOf("PC", "eigenvalue", "proportion of total", "cumulative proportion")
        assertEquals(columns, eigenvalues.columnNames)
        assertEquals("0", eigenvalues[0, 0])
        assertEquals("280", eigenvalues[280, 0])
        assertTrue(abs(eigenvalues[280, CUMULATIVE] as Double - 1.0) < 1e-12)
        val proportions = (0 until eigenvalues.rowCount).sumOf { eigenvalues[it, PROPORTION] as Double }
        assertTrue(abs(proportions - 1.0) < 1e-12, "the proportions summed to $proportions")
    }

    @Test
    fun eigenvectorsHaveOneRowPerSiteNamedByItsPosition() {
        assertEquals(tutorial.numberOfSites, eigenvectors.rowCount)
        assertEquals(listOf("Trait") + (1..5).map { "Eigenvector$it" }, eigenvectors.columnNames)
        val traits = (0 until eigenvectors.rowCount).map { eigenvectors[it, 0] }
        assertEquals(tutorial.positions.map { it.name }, traits)
        assertTrue(eigenvectors[0, 1] is Double, "a ResultTable cell must be a Double, never a Float")
    }

    /** Each eigenvector column is a unit vector over the sites, as the columns of `V` are. */
    @Test
    fun eigenvectorColumnsAreUnitLength() {
        for (component in 1..5) {
            val length = (0 until eigenvectors.rowCount).sumOf {
                val value = eigenvectors[it, component] as Double
                value * value
            }
            assertTrue(abs(length - 1.0) < 1e-12, "column $component has squared length $length")
        }
    }

    @Test
    fun principalComponentsAreFloatCovariatesOverEveryTaxon() {
        assertEquals(6, pcs.numberOfAttributes)
        assertEquals(tutorial.numberOfTaxa, pcs.numberOfObservations)
        assertEquals(PhenotypeTable.AttributeType.taxa, pcs.attributes[0].type)
        assertEquals(tutorial.taxa.names, pcs.attributes[0].values)
        assertEquals((1..5).map { "PC$it" }, pcs.attributes.drop(1).map { it.name })
        assertTrue(pcs.attributes.drop(1).all { it.type == PhenotypeTable.AttributeType.covariate })
        assertTrue(pcs.attributes[1].values.all { it is Float }, "TASSEL 5 stores PCs as floats")
    }

    /**
     * `nTaxa - 1` times the squared length of a principal component is its eigenvalue, which ties
     * the `Float` scores back to the `Double` table they came from.
     */
    @Test
    fun principalComponentLengthsRecoverTheEigenvalues() {
        for (component in 1..5) {
            val scores = pcs.attributes[component].values.map { (it as Float).toDouble() }
            val variance = scores.sumOf { it * it } / (tutorial.numberOfTaxa - 1)
            val expected = eigenvalues[component - 1, VALUE] as Double
            assertTrue(abs(variance - expected) < 1e-5, "PC$component variance $variance, eigenvalue $expected")
        }
    }

    @Test
    fun reportsProgressOverEverySiteOfEveryPass() {
        val passes = reports.also { result }.groupBy { it.message }
        assertEquals(3, passes.size, "expected one pass each for probabilities, the Gram matrix, and V")
        passes.forEach { (message, pass) ->
            assertTrue(pass.all { it.total == tutorial.numberOfSites.toLong() }, message)
            assertEquals(tutorial.numberOfSites.toLong(), pass.last().completed, message)
            assertEquals(pass.map { it.completed }.sorted(), pass.map { it.completed }, message)
        }
    }

    @Test
    fun componentsAreCappedAtTheNumberOfSingularValues() {
        val wide = principalComponents(fourSites, numberOfComponents = 10)
        assertEquals(fourSites.numberOfTaxa + 1, wide.principalComponents.numberOfAttributes)
        assertEquals(fourSites.numberOfTaxa, wide.eigenvalues?.rowCount)
    }

    @Test
    fun eachTableCanBeTurnedOff() {
        val only = principalComponents(fourSites, reportEigenvalues = false, reportEigenvectors = false)
        assertNull(only.eigenvalues)
        assertNull(only.eigenvectors)
        assertEquals(emptyList(), only.tables)
        assertEquals(principalComponents(fourSites).principalComponents.toDelimitedString(), only.toText())
        assertEquals(listOf(eigenvalues, eigenvectors), result.tables)
    }

    /** Correlation PCA scales every column to unit variance, so the eigenvalues sum to the sites. */
    @Test
    fun correlationScalesEveryColumnToUnitVariance() {
        val correlation = principalComponents(fourSites, useCovariance = false, numberOfComponents = 2)
        val table = correlation.eigenvalues!!
        val total = (0 until table.rowCount).sumOf { table[it, VALUE] as Double }
        assertTrue(abs(total - fourSites.numberOfSites) < 1e-9, "eigenvalues summed to $total")
        assertEquals(3, correlation.principalComponents.numberOfAttributes)
    }

    @Test
    fun rejectsInputTooSmallToDecompose() {
        assertFailsWith<IllegalArgumentException> { principalComponents(hapMap("AACC")) }
        assertFailsWith<IllegalArgumentException> { principalComponents(fourSites, numberOfComponents = 0) }
    }

    private companion object {
        const val TUTORIAL = "tutorial/mdp_genotype.hmp.txt"
        const val EIGENVALUES = "golden/pca/eigenvalues.txt"
        const val EIGENVALUES_FULL = "golden/pca/eigenvalues.full.txt"
        const val EIGENVECTORS = "golden/pca/eigenvectors.txt"
        const val EIGENVECTORS_FULL = "golden/pca/eigenvectors.full.txt"
        const val PCS = "golden/pca/pcs.txt"
        const val PCS_FULL = "golden/pca/pcs.full.txt"

        /** Column indices of the eigenvalue table. */
        const val PC = 0
        const val VALUE = 1
        const val PROPORTION = 2
        const val CUMULATIVE = 3

        /** Lines of a `Phenotype` export ahead of the column names: `<Phenotype>` and the types. */
        const val PREAMBLE = 2

        /** Roughly `eps * trace(G)`, below which an eigenvalue of the Gram matrix is only noise. */
        const val NOISE_FLOOR = 1e-12

        val tutorial: GenotypeTable by lazy {
            HapMap.read(Buffer().also { it.write(TestData.readBytes(TUTORIAL)) })
        }

        /** Progress from the one run of [result], which is the only pass over the tutorial data. */
        val reports = mutableListOf<Progress>()

        val result: PrincipalComponentsResult by lazy {
            principalComponents(tutorial, progress = ProgressReporter { reports += it })
        }
        val eigenvalues: ResultTable get() = result.eigenvalues!!
        val eigenvectors: ResultTable get() = result.eigenvectors!!
        val pcs: PhenotypeTable get() = result.principalComponents

        /** Four taxa over four sites, enough to exercise the correlation path. */
        val fourSites: GenotypeTable by lazy { hapMap("AACC", "ACCC", "CCAC", "ACAA") }

        fun PrincipalComponentsResult.toText(): String = principalComponents.toDelimitedString()

        /** Runs [check] for the rounded eigenvalue golden and for its full-precision twin. */
        fun eachEigenvalueGolden(check: (golden: String, actual: String) -> Unit) {
            check(EIGENVALUES, eigenvalues.toDelimitedString())
            check(EIGENVALUES_FULL, fullPrecision(eigenvalues))
        }

        /**
         * A `Phenotype` export is not a table: its first two lines hold one and `n` cells. They are
         * compared as text, and the rest against the golden under the MANIFEST tolerance.
         */
        fun assertPhenotypeMatchesGolden(golden: String, actual: String) {
            val expected = normalizeLineEndings(goldenText(golden))
            assertEquals(expected.lines().take(PREAMBLE), actual.lines().take(PREAMBLE), "$golden preamble")
            assertDelimitedMatches(
                expected = expected.lines().drop(PREAMBLE).joinToString("\n"),
                actual = actual.lines().drop(PREAMBLE).joinToString("\n"),
                tolerance = goldenTolerance(golden),
                source = golden,
            )
        }

        /**
         * [table] with every `Double` spelled in full, which is what the `tassel5-full` goldens
         * hold. `ResultTable.writeDelimited` deliberately rounds through `DoubleFormat`, so it
         * cannot write the `.full.txt` twins.
         */
        fun fullPrecision(table: ResultTable): String = buildString {
            append(table.columnNames.joinToString("\t")).append('\n')
            table.rows.forEach { cells ->
                append(cells.joinToString("\t") { if (it is Double) doubleToString(it) else it.toString() })
                append('\n')
            }
        }

        fun header(text: String): String = normalizeLineEndings(text).lineSequence().first()

        fun dataRows(text: String): List<List<String>> =
            normalizeLineEndings(text).lines().dropLastWhile { it.isEmpty() }.drop(1).map { it.split("\t") }

        fun rejoin(header: String, rows: List<List<String>>): String =
            (listOf(header) + rows.map { it.joinToString("\t") }).joinToString("\n", postfix = "\n")

        /** [text] with only [keep] of its columns, in that order. */
        fun columns(text: String, vararg keep: Int): String {
            val head = header(text).split("\t")
            return rejoin(keep.joinToString("\t") { head[it] }, dataRows(text).map { row -> keep.map { row[it] } })
        }

        /**
         * [text] without the row holding the smallest eigenvalue, the one row the two sides cannot
         * agree on; see [theZeroEigenvalueIsLostInTheGramMatrix]. The row is found by value and not
         * by position, because neither file's row order is meaningful.
         */
        fun withoutTheZeroEigenvalue(text: String): String {
            val rows = dataRows(text)
            val zero = rows.indices.minBy { rows[it][VALUE].toDouble() }
            return rejoin(header(text), rows.filterIndexed { index, _ -> index != zero })
        }

        /** The first [rows] data rows of [text], with its header. */
        fun firstRows(text: String, rows: Int): String = rejoin(header(text), dataRows(text).take(rows))

        /** How many leading data rows of [text] already hold the largest eigenvalues, in order. */
        fun sortedPrefix(text: String): Int {
            val values = dataRows(text).map { it[VALUE].toDouble() }
            val sorted = values.sortedDescending()
            var rows = 0
            while (rows < values.size && values[rows] == sorted[rows]) rows++
            return rows
        }

        /** One HapMap site per character of each [taxa] string, with taxa named `t0`, `t1`, ... */
        fun hapMap(vararg taxa: String): GenotypeTable {
            val names = taxa.indices.map { "t$it" }
            val header = (HapMap.FIXED_COLUMNS + names).joinToString("\t")
            val lines = taxa[0].indices.map { site ->
                val fixed = listOf("s$site", "A/C", "1", "${site + 1}", "+", "NA", "NA", "NA", "NA", "NA", "NA")
                (fixed + taxa.map { "${it[site]}${it[site]}" }).joinToString("\t")
            }
            val text = (listOf(header) + lines).joinToString("\n", postfix = "\n")
            return HapMap.read(Buffer().also { it.write(text.encodeToByteArray()) })
        }
    }
}
