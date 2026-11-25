package net.maizegenetics.tassel.stats

import net.maizegenetics.tassel.testing.TestData
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.ln1p
import kotlin.math.log10
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Log-space tails ([ContinuousDistribution.logSf], [fTestLog10PValue]) against the mpmath grid in
 * `oracle/tails/log-sf.tsv` (ADR-0029 kind 2, rel 1e-12 in log space), closed forms where the tail is elementary, and
 * agreement with [fTestPValue] wherever that is representable.
 */
class LogTailTest {
    private fun relativeError(expected: Double, actual: Double): Double = abs(actual - expected) / abs(expected)

    @Test
    fun logSfMatchesTheMpmathGrid() {
        val rows = TestData.readOracleText("tails/log-sf.tsv").lines().filter { it.isNotEmpty() }.drop(1)
        assertEquals(1076, rows.size)
        val failures = mutableListOf<String>()
        for (row in rows) {
            val cells = row.split('\t')
            val x = cells[3].toDouble()
            val expected = cells[4].toDouble()
            val actual =
                when (cells[0]) {
                    "chisq" -> ChiSquaredDistribution(cells[1].toDouble()).logSf(x)
                    "t" -> StudentTDistribution(cells[1].toDouble()).logSf(x)
                    "f" -> FDistribution(cells[1].toDouble(), cells[2].toDouble()).logSf(x)
                    else -> error("unknown distribution in $row")
                }
            // The fixture's 80-digit mpmath rounds log(1 - 1e-85) to 0.0, so a 0.0 row only bounds the magnitude.
            val ok = if (expected == 0.0) abs(actual) < 1e-79 else relativeError(expected, actual) <= 1e-12
            if (!ok) failures += "$row -> $actual (rel ${relativeError(expected, actual)})"
        }
        assertTrue(
            failures.isEmpty(),
            "${failures.size} of ${rows.size} rows miss rel 1e-12:\n${failures.take(15).joinToString("\n")}",
        )
    }

    @Test
    fun log10SfIsLogSfOverLn10() {
        val f = FDistribution(3.0, 17.0)
        assertEquals(f.logSf(4.2) / ln(10.0), f.log10Sf(4.2), 1e-15)
        assertEquals(fTestLog10PValue(4.2, 3.0, 17.0), f.log10Sf(4.2), 1e-15)
    }

    @Test
    fun logPOfMinusTenThousandIsFinite() {
        // Chi-square with 2 df is exponential: P(X > x) = exp(-x / 2), so ln p = -x / 2 exactly.
        assertEquals(-1e4, ChiSquaredDistribution(2.0).logSf(2e4), 1e-9)
        assertEquals(0.0, ChiSquaredDistribution(2.0).sf(2e4))
        // F with df1 = 2: P(F > f) = (1 + 2 f / df2)^(-df2 / 2).
        val f = 2e4 // ln p = -(df2 / 2) ln1p(2 f / df2)
        val df2 = 8.0
        val expected = -(df2 / 2.0) * ln1p(2.0 * f / df2)
        assertEquals(expected, FDistribution(2.0, df2).logSf(f), 1e-12 * abs(expected))
        assertEquals(expected / ln(10.0), fTestLog10PValue(f, 2.0, df2), 1e-12 * abs(expected))
        // Student t with 1 df is Cauchy: P(T > t) = 1/2 - atan(t) / pi, which is 1 / (pi t) far out.
        assertEquals(-ln(kotlin.math.PI * 1e300), StudentTDistribution(1.0).logSf(1e300), 1e-12 * 690.0)
    }

    @Test
    fun pNearTheSmallestDoubleKeepsItsLog() {
        // exp(-x / 2) hits Double.MIN_VALUE (ln = -744.44) at x = 1488.9; the log stays exact either side of it.
        for (x in listOf(1400.0, 1480.0, 1488.0, 1490.0, 1600.0, 5000.0)) {
            val log = ChiSquaredDistribution(2.0).logSf(x)
            assertEquals(-x / 2.0, log, 1e-12 * abs(log))
        }
        // sf is already 0.0 (Colt's cut-off is e^-709.8, above the denormals) while the log is still exact.
        assertEquals(0.0, ChiSquaredDistribution(2.0).sf(1488.0))
        assertTrue(ChiSquaredDistribution(2.0).sf(1400.0) > 0.0)
        assertEquals(0.0, ChiSquaredDistribution(2.0).sf(1500.0))
        assertTrue(ChiSquaredDistribution(2.0).logSf(1500.0).isFinite())
    }

    @Test
    fun fTestLog10PValueAgreesWithFTestPValueWhereRepresentable() {
        val dfs = listOf(1.0, 2.0, 5.5, 10.0, 30.0, 200.0, 1000.0)
        val fs = listOf(1e-8, 1e-3, 0.1, 0.5, 1.0, 2.0, 5.0, 17.3, 100.0, 1e3, 1e5, 1e8)
        var checked = 0
        for (df1 in dfs) {
            for (df2 in dfs) {
                for (f in fs) {
                    val p = fTestPValue(f, df1, df2)
                    if (p < 1e-290) continue
                    val log10p = fTestLog10PValue(f, df1, df2)
                    // fTestPValue is 1 - (a tiny tail) for p near 1, so it carries 1e-16 / (1 - p) of error the log does not.
                    val tolerance = 1e-12 + 1e-15 / (1.0 - p).coerceAtLeast(1e-300)
                    assertEquals(log10(p), log10p, tolerance, "f=$f df1=$df1 df2=$df2 p=$p")
                    checked++
                }
            }
        }
        assertTrue(checked > 300, "only $checked points compared")
    }

    @Test
    fun fTestLog10PValueDomain() {
        assertEquals(0.0, fTestLog10PValue(0.0, 3.0, 12.0))
        assertEquals(Double.NEGATIVE_INFINITY, fTestLog10PValue(Double.POSITIVE_INFINITY, 3.0, 12.0))
        assertTrue(fTestLog10PValue(Double.NaN, 3.0, 12.0).isNaN())
        assertTrue(fTestLog10PValue(1.0, Double.NaN, 12.0).isNaN())
        assertFailsWith<IllegalArgumentException> { fTestLog10PValue(-1.0, 3.0, 12.0) }
        assertFailsWith<IllegalArgumentException> { fTestLog10PValue(1.0, 0.0, 12.0) }
        assertFailsWith<IllegalArgumentException> { fTestLog10PValue(1.0, 3.0, -2.0) }
        // fTestPValue underflows to 0.0 here; the log10 does not.
        assertEquals(0.0, fTestPValue(1e50, 10.0, 1000.0))
        assertTrue(fTestLog10PValue(1e50, 10.0, 1000.0) < -300.0)
        // An f so large that df1 * f overflows still has a finite log.
        val huge = fTestLog10PValue(1e308, 1000.0, 10.0)
        assertTrue(huge.isFinite() && huge < 0.0, "log10 p = $huge")
    }

    @Test
    fun nonIntegerDegreesOfFreedom() {
        // Wilks and Rao style F approximations: df1 = 2.5, df2 = 7.25. Compared with fTestPValue, and the F tail also
        // follows from the chi-square relation as df2 grows: df1 * F -> chi-square(df1).
        val p = fTestPValue(3.0, 2.5, 7.25)
        assertEquals(log10(p), fTestLog10PValue(3.0, 2.5, 7.25), 1e-13)
        val big = FDistribution(2.5, 1e7).logSf(3.0)
        val chi = ChiSquaredDistribution(2.5).logSf(7.5)
        assertEquals(chi, big, 1e-5 * abs(chi))
        // t and chi-square at non-integer df.
        assertEquals(ln(StudentTDistribution(4.5).sf(2.0)), StudentTDistribution(4.5).logSf(2.0), 1e-13)
        assertEquals(ln(ChiSquaredDistribution(7.5).sf(9.0)), ChiSquaredDistribution(7.5).logSf(9.0), 1e-13)
    }

    @Test
    fun boundariesAndSymmetry() {
        val chi = ChiSquaredDistribution(3.0)
        assertEquals(0.0, chi.logSf(0.0))
        assertEquals(0.0, chi.logSf(-5.0))
        assertEquals(Double.NEGATIVE_INFINITY, chi.logSf(Double.POSITIVE_INFINITY))
        assertTrue(chi.logSf(Double.NaN).isNaN())
        val f = FDistribution(3.0, 4.0)
        assertEquals(0.0, f.logSf(0.0))
        assertEquals(Double.NEGATIVE_INFINITY, f.logSf(Double.POSITIVE_INFINITY))
        assertTrue(f.logSf(Double.NaN).isNaN())
        val t = StudentTDistribution(5.0)
        assertEquals(-ln(2.0), t.logSf(0.0), 1e-15)
        assertEquals(0.0, t.logSf(Double.NEGATIVE_INFINITY))
        assertEquals(Double.NEGATIVE_INFINITY, t.logSf(Double.POSITIVE_INFINITY))
        assertTrue(t.logSf(Double.NaN).isNaN())
        // sf(-x) = 1 - sf(x), so the two logs satisfy exp(a) + exp(b) = 1.
        assertEquals(1.0, exp(t.logSf(1.7)) + exp(t.logSf(-1.7)), 1e-14)
        assertEquals(ln(t.sf(-3.0)), t.logSf(-3.0), 1e-14)
    }

    @Test
    fun logIncompleteFunctionsMatchTheirLinearForms() {
        for ((a, b, x) in listOf(
            Triple(2.0, 3.0, 0.4),
            Triple(0.5, 7.0, 0.9),
            Triple(30.0, 2.5, 0.99),
            Triple(1.0, 1.0, 0.3),
        )) {
            assertEquals(
                ln(SpecialFunctions.regularizedIncompleteBeta(a, b, x)),
                SpecialFunctions.logRegularizedIncompleteBeta(a, b, x),
                1e-13,
                "a=$a b=$b x=$x",
            )
        }
        assertEquals(Double.NEGATIVE_INFINITY, SpecialFunctions.logRegularizedIncompleteBeta(2.0, 3.0, 0.0))
        assertEquals(0.0, SpecialFunctions.logRegularizedIncompleteBeta(2.0, 3.0, 1.0))
        assertFailsWith<IllegalArgumentException> { SpecialFunctions.logRegularizedIncompleteBeta(0.0, 3.0, 0.5) }
        assertFailsWith<IllegalArgumentException> { SpecialFunctions.logRegularizedIncompleteBeta(2.0, 3.0, 1.5) }
        for ((a, x) in listOf(0.5 to 0.2, 3.0 to 2.0, 3.0 to 40.0, 100.0 to 90.0, 100.0 to 130.0)) {
            assertEquals(
                ln(SpecialFunctions.regularizedGammaQ(a, x)),
                SpecialFunctions.logRegularizedGammaQ(a, x),
                1e-12,
                "a=$a x=$x",
            )
        }
        assertEquals(0.0, SpecialFunctions.logRegularizedGammaQ(2.0, 0.0))
        assertEquals(0.0, SpecialFunctions.logRegularizedGammaQ(0.0, 2.0))
        assertTrue(SpecialFunctions.logRegularizedGammaQ(Double.NaN, 2.0).isNaN())
    }
}
