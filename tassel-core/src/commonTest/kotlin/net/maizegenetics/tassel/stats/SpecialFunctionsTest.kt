package net.maizegenetics.tassel.stats

import net.maizegenetics.tassel.testing.assertMatchesGolden
import net.maizegenetics.tassel.testing.goldenText
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * [SpecialFunctions] and [fTestPValue] against Colt 1.2.0 `cern.jet.stat.Gamma` values from the pinned TASSEL 5
 * standalone (`golden/stats/special-functions.txt`, rel 1e-13, ADR-0009 and ADR-0020), plus the control-flow branches
 * the grid reaches.
 *
 * @see net.maizegenetics.stats.linearmodels.LinearModelUtils.Ftest
 */
class SpecialFunctionsTest {
    @Test
    fun matchesColtOnTheGoldenGrid() {
        val rows = goldenText("golden/stats/special-functions.txt").trim().lines()
        val actual =
            buildString {
                append(rows[0]).append('\n')
                for (row in rows.drop(1)) {
                    val cells = row.split('\t')
                    val a = cells[1].toDoubleOrNull()
                    val b = cells[2].toDoubleOrNull()
                    val x = cells[3].toDoubleOrNull()
                    val value =
                        try {
                            when (val function = cells[0]) {
                                "incompleteBeta" -> SpecialFunctions.regularizedIncompleteBeta(a!!, b!!, x!!)
                                "incompleteGamma" -> SpecialFunctions.regularizedGammaP(a!!, x!!)
                                "incompleteGammaComplement" -> SpecialFunctions.regularizedGammaQ(a!!, x!!)
                                "gamma" -> SpecialFunctions.gamma(a!!)
                                "logGamma" -> SpecialFunctions.logGamma(a!!)
                                else -> error("unknown function $function")
                            }.toString()
                        } catch (e: IllegalArgumentException) {
                            "error"
                        } catch (e: ArithmeticException) {
                            "error"
                        }
                    append(cells.take(4).joinToString("\t")).append('\t').append(value).append('\n')
                }
            }
        assertMatchesGolden("golden/stats/special-functions.txt", actual)
    }

    @Test
    fun betaPowerSeriesBranch() {
        // b * x <= 1 and x <= 0.95: I_x(2, 1) = x^2.
        assertEquals(0.25, SpecialFunctions.regularizedIncompleteBeta(2.0, 1.0, 0.5), 1e-15)
        // I_x(1, 1) = x.
        assertEquals(0.3, SpecialFunctions.regularizedIncompleteBeta(1.0, 1.0, 0.3), 1e-15)
    }

    @Test
    fun betaSwapBranchAndCap() {
        // x above the mean a / (a + b) swaps a and b; I_x(a, b) = 1 - I_(1-x)(b, a).
        val direct = SpecialFunctions.regularizedIncompleteBeta(2.0, 3.0, 0.8)
        val swapped = 1.0 - SpecialFunctions.regularizedIncompleteBeta(3.0, 2.0, 0.2)
        assertEquals(swapped, direct, 1e-14)
        // A swapped tail at or below machine epsilon is held at 1 - 2^-53, never 1.
        assertEquals(0.9999999999999999, SpecialFunctions.regularizedIncompleteBeta(1.0, 1000.0, 0.5))
    }

    @Test
    fun betaUnderflowsToZeroInLogSpace() {
        assertEquals(0.0, SpecialFunctions.regularizedIncompleteBeta(1000.0, 1000.0, 1e-100))
        assertEquals(0.0, SpecialFunctions.regularizedIncompleteBeta(1e7, 1e7, 1e-5))
    }

    @Test
    fun betaEndpoints() {
        assertEquals(0.0, SpecialFunctions.regularizedIncompleteBeta(2.0, 3.0, 0.0))
        assertEquals(1.0, SpecialFunctions.regularizedIncompleteBeta(2.0, 3.0, 1.0))
    }

    @Test
    fun betaDomainErrorsThrow() {
        assertFailsWith<IllegalArgumentException> { SpecialFunctions.regularizedIncompleteBeta(0.0, 1.0, 0.5) }
        assertFailsWith<IllegalArgumentException> { SpecialFunctions.regularizedIncompleteBeta(1.0, -1.0, 0.5) }
        assertFailsWith<IllegalArgumentException> { SpecialFunctions.regularizedIncompleteBeta(1.0, 1.0, -0.1) }
        assertFailsWith<IllegalArgumentException> { SpecialFunctions.regularizedIncompleteBeta(1.0, 1.0, 1.1) }
    }

    @Test
    fun nanInputGivesNaN() {
        val nan = Double.NaN
        assertTrue(SpecialFunctions.regularizedIncompleteBeta(nan, 1.0, 0.5).isNaN())
        assertTrue(SpecialFunctions.regularizedIncompleteBeta(1.0, nan, 0.5).isNaN())
        assertTrue(SpecialFunctions.regularizedIncompleteBeta(1.0, 1.0, nan).isNaN())
        assertTrue(SpecialFunctions.regularizedGammaP(nan, 1.0).isNaN())
        assertTrue(SpecialFunctions.regularizedGammaP(1.0, nan).isNaN())
        assertTrue(SpecialFunctions.regularizedGammaQ(nan, 1.0).isNaN())
        assertTrue(SpecialFunctions.regularizedGammaQ(1.0, nan).isNaN())
        assertTrue(SpecialFunctions.gamma(nan).isNaN())
        assertTrue(SpecialFunctions.logGamma(nan).isNaN())
        assertTrue(fTestPValue(nan, 2.0, 10.0).isNaN())
    }

    @Test
    fun gammaOfIntegersIsFactorial() {
        assertEquals(1.0, SpecialFunctions.gamma(1.0))
        assertEquals(24.0, SpecialFunctions.gamma(5.0), 1e-12)
        assertEquals(1.329340388179137, SpecialFunctions.gamma(2.5), 1e-12)
        assertFailsWith<ArithmeticException> { SpecialFunctions.gamma(0.0) }
        assertFailsWith<ArithmeticException> { SpecialFunctions.gamma(-2.0) }
    }

    @Test
    fun incompleteGammaIdentities() {
        // P(1, x) = 1 - exp(-x); P + Q = 1.
        assertEquals(1.0 - kotlin.math.exp(-2.0), SpecialFunctions.regularizedGammaP(1.0, 2.0), 1e-15)
        assertEquals(
            1.0,
            SpecialFunctions.regularizedGammaP(3.5, 10.0) + SpecialFunctions.regularizedGammaQ(3.5, 10.0),
            1e-15,
        )
        assertEquals(0.0, SpecialFunctions.regularizedGammaP(2.0, 0.0))
        assertEquals(1.0, SpecialFunctions.regularizedGammaQ(2.0, 0.0))
        assertEquals(0.0, SpecialFunctions.regularizedGammaQ(1.0, 1e8))
    }

    @Test
    fun fTestPValueEdges() {
        assertEquals(1.0, fTestPValue(0.0, 3.0, 20.0))
        assertEquals(0.0, fTestPValue(Double.POSITIVE_INFINITY, 3.0, 20.0))
        assertFailsWith<IllegalArgumentException> { fTestPValue(-1.0, 3.0, 20.0) }
    }

    @Test
    fun fTestPValueKnownValues() {
        // F(1, 1) upper tail at 1 is exactly 0.5; F(2, d2) has the closed form (d2 / (d2 + 2 f))^(d2 / 2).
        assertEquals(0.5, fTestPValue(1.0, 1.0, 1.0), 1e-15)
        assertEquals(kotlin.math.exp(5.0 * kotlin.math.ln(10.0 / 16.0)), fTestPValue(3.0, 2.0, 10.0), 1e-14)
    }

    @Test
    fun fTestPValueAcceptsNonIntegerDegreesOfFreedom() {
        val p = fTestPValue(2.7, 3.25, 21.0)
        assertTrue(p > fTestPValue(2.7, 3.0, 21.0).coerceAtMost(fTestPValue(2.7, 4.0, 21.0)) - 1e-12)
        assertTrue(p < fTestPValue(2.7, 3.0, 21.0).coerceAtLeast(fTestPValue(2.7, 4.0, 21.0)) + 1e-12)
        assertTrue(abs(p) in 0.0..1.0)
    }
}
