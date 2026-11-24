package net.maizegenetics.tassel.stats

import net.maizegenetics.tassel.testing.assertMatchesGolden
import net.maizegenetics.tassel.testing.goldenText
import kotlin.math.PI
import kotlin.math.abs
import kotlin.math.atan
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.sqrt
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * [FDistribution], [StudentTDistribution], and [ChiSquaredDistribution] against Colt 1.2.0 `Probability` tails
 * (`golden/stats/distributions.txt`, rel 1e-12) and commons-math3 3.4.1 quantiles
 * (`golden/stats/distribution-quantiles.txt`, abs 1e-9), then closed forms, validation, the support, and
 * `cdf(quantile(p)) = p` round trips.
 *
 * @see cern.jet.stat.Probability
 * @see org.apache.commons.math3.distribution.FDistribution
 */
class DistributionsTest {
    private fun rows(path: String): List<List<String>> = goldenText(path).trim().lines().map { it.split('\t') }

    private fun regenerate(path: String, value: (String, Double, Double?, Double) -> Double): String {
        val rows = rows(path)
        return buildString {
            append(rows[0].joinToString("\t")).append('\n')
            for (cells in rows.drop(1)) {
                val a = cells[1].toDouble()
                val b = cells[2].toDoubleOrNull()
                val x = cells[3].toDouble()
                append(cells.take(4).joinToString("\t")).append('\t').append(value(cells[0], a, b, x)).append('\n')
            }
        }
    }

    @Test
    fun cdfAndSfMatchColtOnTheGoldenGrid() {
        val path = "golden/stats/distributions.txt"
        val actual =
            regenerate(path) { function, a, b, x ->
                when (function) {
                    "chiSquareCdf" -> ChiSquaredDistribution(a).cdf(x)
                    "chiSquareSf" -> ChiSquaredDistribution(a).sf(x)
                    "tCdf" -> StudentTDistribution(a).cdf(x)
                    "tSf" -> StudentTDistribution(a).sf(x)
                    "fCdf" -> FDistribution(a, b!!).cdf(x)
                    "fSf" -> FDistribution(a, b!!).sf(x)
                    else -> error("unknown function $function")
                }
            }
        assertMatchesGolden(path, actual)
    }

    @Test
    fun quantilesMatchCommonsMathOnTheGoldenGrid() {
        val path = "golden/stats/distribution-quantiles.txt"
        val actual =
            regenerate(path) { function, a, b, p ->
                when (function) {
                    "chiSquareQuantile" -> ChiSquaredDistribution(a).quantile(p)
                    "tQuantile" -> StudentTDistribution(a).quantile(p)
                    "fQuantile" -> FDistribution(a, b!!).quantile(p)
                    else -> error("unknown function $function")
                }
            }
        assertMatchesGolden(path, actual)
    }

    @Test
    fun parametersAreValidated() {
        for (bad in listOf(0.0, -1.0, Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY)) {
            assertEquals(
                "df must be a finite number greater than 0, got $bad",
                assertFailsWith<IllegalArgumentException> { ChiSquaredDistribution(bad) }.message,
            )
            assertEquals(
                "df must be a finite number greater than 0, got $bad",
                assertFailsWith<IllegalArgumentException> { StudentTDistribution(bad) }.message,
            )
            assertEquals(
                "df1 must be a finite number greater than 0, got $bad",
                assertFailsWith<IllegalArgumentException> { FDistribution(bad, 3.0) }.message,
            )
            assertEquals(
                "df2 must be a finite number greater than 0, got $bad",
                assertFailsWith<IllegalArgumentException> { FDistribution(3.0, bad) }.message,
            )
        }
    }

    @Test
    fun probabilitiesOutsideZeroToOneAreErrors() {
        for (d in all()) {
            for (bad in listOf(-0.1, 1.1, Double.NaN, Double.POSITIVE_INFINITY)) {
                assertEquals(
                    "p must be in [0, 1], got $bad",
                    assertFailsWith<IllegalArgumentException> {
                        d.quantile(bad)
                    }.message,
                )
                assertEquals(
                    "q must be in [0, 1], got $bad",
                    assertFailsWith<IllegalArgumentException> {
                        d.inverseSf(bad)
                    }.message,
                )
            }
        }
    }

    @Test
    fun endpointsReturnTheSupportBounds() {
        for (d in listOf(ChiSquaredDistribution(3.0), FDistribution(3.0, 7.0))) {
            assertEquals(0.0, d.quantile(0.0))
            assertEquals(Double.POSITIVE_INFINITY, d.quantile(1.0))
            assertEquals(Double.POSITIVE_INFINITY, d.inverseSf(0.0))
            assertEquals(0.0, d.inverseSf(1.0))
        }
        val t = StudentTDistribution(5.0)
        assertEquals(Double.NEGATIVE_INFINITY, t.quantile(0.0))
        assertEquals(Double.POSITIVE_INFINITY, t.quantile(1.0))
        assertEquals(Double.POSITIVE_INFINITY, t.inverseSf(0.0))
        assertEquals(Double.NEGATIVE_INFINITY, t.inverseSf(1.0))
        assertEquals(0.0, t.quantile(0.5))
        assertEquals(0.0, t.inverseSf(0.5))
    }

    @Test
    fun xOutsideTheSupport() {
        for (d in listOf(ChiSquaredDistribution(3.0), FDistribution(3.0, 7.0))) {
            assertEquals(0.0, d.cdf(-1.0))
            assertEquals(1.0, d.sf(-1.0))
            assertEquals(0.0, d.pdf(-1.0))
            assertEquals(Double.NEGATIVE_INFINITY, d.logPdf(-1.0))
            assertEquals(0.0, d.cdf(0.0))
            assertEquals(1.0, d.sf(0.0))
            assertEquals(1.0, d.cdf(Double.POSITIVE_INFINITY))
            assertEquals(0.0, d.sf(Double.POSITIVE_INFINITY))
            assertEquals(0.0, d.pdf(Double.POSITIVE_INFINITY))
        }
        val t = StudentTDistribution(3.0)
        assertEquals(0.0, t.cdf(Double.NEGATIVE_INFINITY))
        assertEquals(1.0, t.cdf(Double.POSITIVE_INFINITY))
        assertEquals(0.0, t.sf(Double.POSITIVE_INFINITY))
        assertEquals(1.0, t.sf(Double.NEGATIVE_INFINITY))
        assertEquals(0.0, t.pdf(Double.POSITIVE_INFINITY))
        for (d in all()) {
            assertTrue(d.cdf(Double.NaN).isNaN())
            assertTrue(d.sf(Double.NaN).isNaN())
            assertTrue(d.pdf(Double.NaN).isNaN())
            assertTrue(d.logPdf(Double.NaN).isNaN())
        }
    }

    @Test
    fun densityAtTheLowerBound() {
        assertEquals(Double.POSITIVE_INFINITY, ChiSquaredDistribution(1.0).pdf(0.0))
        assertEquals(0.5, ChiSquaredDistribution(2.0).pdf(0.0), 1e-15)
        assertEquals(0.0, ChiSquaredDistribution(3.0).pdf(0.0))
        assertEquals(Double.POSITIVE_INFINITY, FDistribution(1.0, 4.0).pdf(0.0))
        assertEquals(1.0, FDistribution(2.0, 4.0).pdf(0.0), 1e-15)
        assertEquals(0.0, FDistribution(3.0, 4.0).pdf(0.0))
    }

    @Test
    fun oneDegreeOfFreedomClosedForms() {
        // Cauchy: cdf = 1/2 + atan(t) / pi, pdf = 1 / (pi (1 + t^2)).
        val t = StudentTDistribution(1.0)
        for (x in listOf(-50.0, -2.0, -0.3, 0.0, 0.7, 4.0, 1e3)) {
            assertEquals(0.5 + atan(x) / PI, t.cdf(x), 1e-15)
            assertEquals(1.0 / (PI * (1.0 + x * x)), t.pdf(x), 1e-15)
        }
        // The Cauchy upper tail 1/2 - atan(x) / pi = atan(1/x) / pi stays relatively accurate far out.
        assertRelative(atan(1e-150) / PI, t.sf(1e150), 1e-13)
        // chi-square with one df is Z^2: sf(x) = erfc(sqrt(x / 2)); at x = 1 that is 0.31731050786291415.
        assertRelative(0.31731050786291415, ChiSquaredDistribution(1.0).sf(1.0), 1e-14)
        // F(1, 1) is the ratio of two chi-square(1): cdf = 2 atan(sqrt(x)) / pi.
        val f = FDistribution(1.0, 1.0)
        for (x in listOf(0.01, 0.5, 1.0, 3.0, 400.0)) {
            assertEquals(2.0 * atan(sqrt(x)) / PI, f.cdf(x), 1e-15)
            assertEquals(1.0 / (PI * sqrt(x) * (1.0 + x)), f.pdf(x), 1e-14)
        }
    }

    @Test
    fun chiSquareTwoDegreesOfFreedomIsExponential() {
        val d = ChiSquaredDistribution(2.0)
        for (x in listOf(0.01, 1.0, 5.0, 100.0, 1000.0, 1380.0)) {
            assertRelative(exp(-x / 2.0), d.sf(x), 1e-13)
            assertRelative(exp(-x / 2.0) / 2.0, d.pdf(x), 1e-13)
        }
        assertRelative(-ln(2.0), d.logPdf(0.0), 1e-15)
        assertRelative(2.0 * -ln(1e-300), d.inverseSf(1e-300), 1e-13)
    }

    @Test
    fun largeDegreesOfFreedomApproachTheNormal() {
        val t = StudentTDistribution(1e7)
        // Phi(1.959963984540054) = 0.975 and Phi(-8) = 6.22096057427178e-16; the t differs by O(1/df).
        assertEquals(0.975, t.cdf(1.959963984540054), 1e-6)
        assertEquals(1.959963984540054, t.quantile(0.975), 1e-5)
        // Computed by integrating the density: 6.2216174e-16 at df = 1e7, 6.2869599e-16 at 1e5 (normal 6.2209606e-16).
        assertRelative(6.221617412938895e-16, t.sf(8.0), 1e-6)
        assertRelative(6.286959938016749e-16, StudentTDistribution(1e5).sf(8.0), 1e-8)
        assertRelative(1.0 / sqrt(2.0 * PI), t.pdf(0.0), 1e-7)
        // The log density keeps its accuracy at huge df, where lnGamma differences would cancel.
        assertRelative(1.0 / sqrt(2.0 * PI), StudentTDistribution(1e12).pdf(0.0), 1e-9)
        // F(df1, infinity) is chi-square(df1) / df1, and a huge df2 gets close.
        assertRelative(ChiSquaredDistribution(4.0).sf(8.0), FDistribution(4.0, 1e8).sf(2.0), 1e-6)
    }

    @Test
    fun studentTIsSymmetric() {
        val t = StudentTDistribution(7.5)
        for (x in listOf(0.0, 1e-9, 0.4, 2.0, 11.0, 1e4)) {
            assertEquals(t.cdf(-x), t.sf(x), 0.0)
            assertEquals(t.pdf(-x), t.pdf(x), 0.0)
            assertEquals(1.0, t.cdf(x) + t.cdf(-x), 1e-15)
        }
        assertEquals(0.5, t.cdf(0.0))
        // Both tails are accurate: the lower tail of a very negative t is as small as the upper tail of its mirror.
        assertTrue(t.cdf(-1e5) > 0.0 && t.cdf(-1e5) < 1e-20)
        assertEquals(t.sf(1e5), t.cdf(-1e5), 0.0)
        // Near the centre the small t is kept: cdf(1e-12) - 1/2 is pdf(0) * 1e-12.
        assertRelative(t.pdf(0.0) * 1e-12, t.cdf(1e-12) - 0.5, 1e-3)
    }

    @Test
    fun fSurvivalFunctionIsTheTasselFTestPValue() {
        val d = FDistribution(3.0, 278.0)
        for (f in listOf(0.01, 1.0, 2.7, 15.0, 400.0, 1e6)) {
            assertEquals(fTestPValue(f, 3.0, 278.0), d.sf(f), 0.0)
        }
    }

    @Test
    fun extremeTailsStayRelativelyAccurate() {
        // sf is the small side computed directly, so these are nowhere near 1 - cdf rounding to zero.
        // F(2, 2): sf = 1 / (1 + x).
        assertRelative(1e-298, FDistribution(2.0, 2.0).sf(1e298), 1e-12)
        assertRelative(1.0 / (1.0 + 1e100), FDistribution(2.0, 2.0).sf(1e100), 1e-12)
        val t = StudentTDistribution(2.0)
        // t(2): sf = 1/2 - x / (2 sqrt(2 + x^2)), which is 1 / (2 x^2) to leading order.
        assertRelative(1.0 / (2.0 * 1e300), t.sf(1e150), 1e-12)
        assertEquals(1.0, ChiSquaredDistribution(5.0).cdf(1500.0))
        assertTrue(ChiSquaredDistribution(5.0).sf(1000.0) in 1e-230..1e-200)
    }

    @Test
    fun quantileSolvesTheFarTails() {
        // chi-square(2): sf(x) = exp(-x / 2), so inverseSf(q) = -2 ln q.
        assertRelative(-2.0 * ln(1e-300), ChiSquaredDistribution(2.0).inverseSf(1e-300), 1e-13)
        assertRelative(-2.0 * ln(TWO_POW_M53), ChiSquaredDistribution(2.0).quantile(1.0 - TWO_POW_M53), 1e-13)
        // F(2, 2): sf = 1 / (1 + x), so inverseSf(q) = 1 / q - 1.
        assertRelative(1e300, FDistribution(2.0, 2.0).inverseSf(1e-300), 1e-13)
        // t(2): sf(x) is about 1 / (2 x^2) far out, so inverseSf(1e-300) is about 1 / sqrt(2e-300) = 7.07e149.
        assertRelative(sqrt(0.5e300), StudentTDistribution(2.0).inverseSf(1e-300), 1e-12)
    }

    @Test
    fun quantileRoundTrips() {
        val ps =
            (-300..-1 step 7).map { exp10(it) } +
                listOf(0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0 - 1e-8, 1.0 - 1e-12, 1.0 - 1e-16)
        for (d in all()) {
            for (p in ps) {
                val x = d.quantile(p)
                // A root below about 1e-290 is where the lower tail itself underflows, and one beyond the doubles is
                // infinite; neither can round-trip.
                if (!x.isFinite() || abs(x) < 1e-290) continue
                val back = if (p < 0.5) d.cdf(x) else 1.0 - d.sf(x)
                if (p <
                    0.5
                ) {
                    assertRelative(p, back, 1e-12, "$d cdf(quantile($p)), x = $x")
                } else {
                    assertEquals(p, back, 1e-12)
                }
                // The same tail from the other end: sf(inverseSf(q)) = q.
                val q = p
                val y = d.inverseSf(q)
                if (!y.isFinite() || abs(y) < 1e-290) continue
                if (q <=
                    0.5
                ) {
                    assertRelative(q, d.sf(y), 1e-12, "$d sf(inverseSf($q)), y = $y")
                } else {
                    assertEquals(1.0 - q, d.cdf(y), 1e-12)
                }
            }
        }
    }

    @Test
    fun quantileAndInverseSfAreMirrorImages() {
        for (d in all()) {
            // 1 - p is exact for these p, so the two solves target the same tail.
            for (p in listOf(TWO_POW_M30, TWO_POW_M10, 0.25, 0.5, 0.75, 1.0 - TWO_POW_M10)) {
                assertRelative(d.quantile(p), d.inverseSf(1.0 - p), 1e-12, "$d p=$p")
            }
        }
    }

    private fun all(): List<ContinuousDistribution> = listOf(
        ChiSquaredDistribution(1.0),
        ChiSquaredDistribution(2.5),
        ChiSquaredDistribution(300.0),
        StudentTDistribution(1.0),
        StudentTDistribution(4.5),
        StudentTDistribution(1000.0),
        FDistribution(1.0, 1.0),
        FDistribution(3.0, 278.0),
        FDistribution(12.0, 2.5),
        FDistribution(500.0, 500.0),
    )

    private companion object {
        const val TWO_POW_M10 = 1.0 / 1024.0
        const val TWO_POW_M30 = 1.0 / 1073741824.0
        const val TWO_POW_M53 = 1.0 / 9007199254740992.0
    }

    private fun exp10(power: Int): Double = "1e$power".toDouble()

    private fun assertRelative(expected: Double, actual: Double, tolerance: Double, message: String? = null) {
        assertTrue(
            abs(actual - expected) <= tolerance * abs(expected),
            "${message ?: "values differ"}: expected $expected, got $actual",
        )
    }
}
