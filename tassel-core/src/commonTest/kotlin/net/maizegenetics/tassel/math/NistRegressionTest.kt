package net.maizegenetics.tassel.math

import net.maizegenetics.tassel.testing.TestData
import kotlin.math.abs
import kotlin.math.log10
import kotlin.math.sqrt
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Fits the NIST StRD linear regression datasets (`testdata/oracle/nist/`, ADR-0029 kind 2) with [HouseholderQr] and
 * compares coefficients, standard errors, and the residual standard deviation with NIST's certified values.
 *
 * The fit is the one the model layer will use: pivoted QR, `b = solve(y)`, residual sum of squares from the tail of
 * `Q'y`, and standard errors `sigma * rInverseRowNorms`, with no normal equations and no inverse. Agreement is
 * counted in significant digits, `-log10(|actual - certified| / |certified|)`. Where NIST certifies zero (Wampler1
 * and Wampler2 have exact data, so a zero residual and zero standard errors) the comparison is absolute.
 *
 * The fit is [HouseholderQr.solveRefined] (Björck refinement with double-double residuals), because one plain pass
 * cannot reach the targets on the badly conditioned sets: LAPACK (numpy `lstsq` and `qr`, checked once) gets 10.9,
 * 8.0, and 6.0 digits on Longley, Wampler4, and Wampler5, and this kernel's plain `solve` gets 10.9, 8.6, and 6.5
 * ([plainSolveKeepsItsDigitsOnTheHardSets] pins that). Required digits, on coefficients, standard errors, and
 * residual SD together:
 *
 * - [LONGLEY_DIGITS] (12) for Longley, standard errors included (`rInverseRowNorms(a)` refines them too);
 * - [DIGITS] (9) for every other dataset except Filip, Wampler4 and Wampler5 included;
 * - Filip (degree-10 polynomial, condition number near 1e15) either reaches [FILIP_DIGITS] or reports rank
 *   deficiency; a full-rank answer with fewer digits fails.
 */
class NistRegressionTest {
    @Test
    fun longley() = check("Longley", LONGLEY_DIGITS)

    @Test
    fun norris() = check("Norris", DIGITS)

    @Test
    fun pontius() = check("Pontius", DIGITS)

    @Test
    fun noInt1() = check("NoInt1", DIGITS)

    @Test
    fun noInt2() = check("NoInt2", DIGITS)

    @Test
    fun wampler1() = check("Wampler1", DIGITS)

    @Test
    fun wampler2() = check("Wampler2", DIGITS)

    @Test
    fun wampler3() = check("Wampler3", DIGITS)

    @Test
    fun wampler4() = check("Wampler4", DIGITS)

    @Test
    fun wampler5() = check("Wampler5", DIGITS)

    @Test
    fun filipReachesSevenDigitsOrReportsRankDeficiency() {
        val fit = fit("Filip")
        if (fit.qr.rank < fit.parameters) {
            assertTrue(fit.qr.aliasedColumns.isNotEmpty(), "rank ${fit.qr.rank} but no aliased column is reported")
            return
        }
        val digits = fit.digits()
        assertTrue(digits >= FILIP_DIGITS, "Filip is full rank but agrees to only $digits digits: ${fit.report()}")
    }

    @Test
    fun plainSolveKeepsItsDigitsOnTheHardSets() {
        // A regression guard for the unrefined kernel, just under what it reaches today.
        for ((name, floor) in listOf("Longley" to 10.5, "Wampler4" to 8.5, "Wampler5" to 6.5)) {
            val digits = fit(name, refine = false).digits()
            assertTrue(digits >= floor, "$name plain solve agrees to $digits digits, floor $floor")
        }
    }

    private fun check(name: String, required: Double) {
        val fit = fit(name)
        assertEquals(fit.parameters, fit.qr.rank, "$name must be detected as full rank")
        val digits = fit.digits()
        assertTrue(digits >= required, "$name agrees to $digits digits, needs $required: ${fit.report()}")
    }

    private class Fit(
        val name: String,
        val parameters: Int,
        val qr: HouseholderQr,
        val coefficients: DoubleArray,
        val standardErrors: DoubleArray,
        val residualSd: Double,
        val certified: Map<String, String>,
        val intercept: Boolean,
    ) {
        fun comparisons(): List<Triple<String, Double, Double>> {
            val rows = mutableListOf<Triple<String, Double, Double>>()
            for (j in 0 until parameters) {
                val (estimate, error) = certified.getValue("B${if (intercept) j else j + 1}").split(' ').map {
                    it.toDouble()
                }
                rows += Triple("B$j", coefficients[j], estimate)
                rows += Triple("SE(B$j)", standardErrors[j], error)
            }
            rows += Triple("residual_sd", residualSd, certified.getValue("residual_sd").toDouble())
            return rows
        }

        fun digits(standardErrors: Boolean = true): Double = comparisons()
            .filter { standardErrors || !it.first.startsWith("SE") }
            .minOf { (_, actual, expected) -> digitsOf(actual, expected) }

        fun report(): String = comparisons().joinToString { (label, actual, expected) ->
            "$label=${digitsOf(actual, expected).toString().take(5)}"
        }
    }

    private fun fit(name: String, refine: Boolean = true): Fit {
        val certified = TestData.readOracleText("nist/$name.certified.txt").lines()
            .filter { " = " in it }
            .associate { it.substringBefore(" = ") to it.substringAfter(" = ").trim() }
        val table = TestData.readOracleText("nist/$name.data.tsv").lines().filter { it.isNotEmpty() }
        val header = table.first().split('\t')
        val rows = table.drop(1).map { line -> line.split('\t').map { it.toDouble() } }
        val n = rows.size
        assertEquals(certified.getValue("observations").toInt(), n)
        val p = certified.getValue("parameters").toInt()
        val intercept = certified.getValue("intercept") == "true"

        val design = DenseMatrix(n, p)
        val xColumns = header.size - 1
        for (i in 0 until n) {
            for (j in 0 until p) {
                design[i, j] = when {
                    // Several predictors (Longley): an intercept, then each x as given.
                    xColumns > 1 -> if (j == 0) 1.0 else rows[i][j]
                    // One predictor: powers of x, starting at x^0 with an intercept and x^1 without.
                    else -> powi(rows[i][1], if (intercept) j else j + 1)
                }
            }
        }
        val y = DoubleArray(n) { rows[it][0] }

        val qr = HouseholderQr.of(design)
        val coefficients: DoubleArray
        var rss = 0.0
        if (refine) {
            val refined = qr.solveRefinedWithResiduals(design, y)
            coefficients = refined.coefficients
            for (r in refined.residuals) rss += r * r
        } else {
            coefficients = qr.solve(y)
            val qty = y.copyOf().also { qr.applyQt(it) }
            for (i in qr.rank until n) rss += qty[i] * qty[i]
        }
        val sigma = sqrt(rss / (n - p))
        val norms = if (refine) qr.rInverseRowNorms(design) else qr.rInverseRowNorms()
        return Fit(name, p, qr, coefficients, DoubleArray(p) { sigma * norms[it] }, sigma, certified, intercept)
    }

    private fun powi(x: Double, power: Int): Double {
        var result = 1.0
        repeat(power) { result *= x }
        return result
    }

    companion object {
        const val DIGITS = 9.0
        const val LONGLEY_DIGITS = 12.0
        const val FILIP_DIGITS = 7.0

        /**
         * Significant digits of agreement, capped at 15. Where [expected] is zero the error is measured against 1,
         * that is, as an absolute error.
         */
        fun digitsOf(actual: Double, expected: Double): Double {
            if (actual == expected) return 15.0
            val scale = if (expected == 0.0) 1.0 else abs(expected)
            val relative = abs(actual - expected) / scale
            return if (!relative.isFinite()) 0.0 else minOf(15.0, -log10(relative))
        }
    }
}
