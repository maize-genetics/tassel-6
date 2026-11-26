package net.maizegenetics.tassel.distance

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertFalse
import kotlin.test.assertIs
import kotlin.test.assertNotNull
import kotlin.test.assertNull
import kotlin.test.assertTrue

/** [KinshipValidation]. */
class KinshipValidationTest {
    private fun matrix(vararg rows: FloatArray): DistanceMatrix {
        val values = rows.flatMap { it.toList() }.toFloatArray()
        return DistanceMatrix(TaxaList.of((0 until rows.size).map { "T$it" }), values)
    }

    @Test
    fun acceptsAPositiveDefiniteMatrix() = runTest {
        val report = KinshipValidation.check(matrix(floatArrayOf(2f), floatArrayOf(1f, 2f)))
        assertTrue(report.isValid)
        assertEquals(3.0, report.maxEigenvalue!!, 1e-12)
        assertEquals(1.0, report.minEigenvalue!!, 1e-12)
        assertEquals(0.0, report.maxAsymmetry)
        report.requireValid()
    }

    @Test
    fun acceptsASingularButPositiveSemidefiniteMatrix() = runTest {
        // Rank one: v v' with v = (1, 2, 3). Its smallest eigenvalue is zero up to rounding.
        val report = KinshipValidation.check(
            matrix(floatArrayOf(1f), floatArrayOf(2f, 4f), floatArrayOf(3f, 6f, 9f)),
        )
        assertTrue(report.isValid, report.problems.toString())
        assertEquals(0.0, report.minEigenvalue!!, 1e-10)
        assertEquals(14.0, report.maxEigenvalue!!, 1e-10)
    }

    @Test
    fun acceptsAZeroMatrix() = runTest {
        assertTrue(KinshipValidation.check(matrix(floatArrayOf(0f), floatArrayOf(0f, 0f))).isValid)
    }

    @Test
    fun rejectsAnIndefiniteMatrix() = runTest {
        // Eigenvalues 3 and -1.
        val report = KinshipValidation.check(matrix(floatArrayOf(1f), floatArrayOf(2f, 1f)))
        assertFalse(report.isValid)
        val problem = assertIs<KinshipProblem.NotPositiveSemidefinite>(report.problems.single())
        assertEquals(-1.0, problem.minEigenvalue, 1e-12)
        assertEquals(-3e-6, problem.threshold, 1e-18)
        assertContains(
            assertFailsWith<IllegalArgumentException> {
                report.requireValid()
            }.message.orEmpty(),
            "positive semidefinite",
        )
    }

    @Test
    fun toleranceDecidesHowNegativeIsTooNegative() = runTest {
        // Eigenvalues 1 + 1e-3 and -1e-3 relative to... scaled: 2.001 and -0.001, ratio about -5e-4.
        val nearly = matrix(floatArrayOf(1f), floatArrayOf(1.001f, 1f))
        assertFalse(KinshipValidation.check(nearly, tolerance = 1e-6).isValid)
        assertTrue(KinshipValidation.check(nearly, tolerance = 1e-3).isValid)
    }

    @Test
    fun rejectsANegativeDefiniteMatrix() = runTest {
        assertFalse(KinshipValidation.check(matrix(floatArrayOf(-1f))).isValid)
    }

    @Test
    fun reportsNaNAndInfinityAndSkipsTheEigenvalues() = runTest {
        val report = KinshipValidation.check(matrix(floatArrayOf(1f), floatArrayOf(Float.NaN, 1f)))
        assertFalse(report.isValid)
        val problem = assertIs<KinshipProblem.NonFinite>(report.problems.single())
        assertEquals(2, problem.count)
        assertEquals(0, problem.row)
        assertEquals(1, problem.column)
        assertNull(report.minEigenvalue)
        assertFalse(KinshipValidation.check(matrix(floatArrayOf(Float.POSITIVE_INFINITY))).isValid)
    }

    @Test
    fun reportsAnAsymmetricSquareOfValues() = runTest {
        val square = arrayOf(doubleArrayOf(1.0, 0.5), doubleArrayOf(0.2, 1.0))
        val report = KinshipValidation.check(2) { row, column -> square[row][column] }
        assertFalse(report.isValid)
        val problem = assertIs<KinshipProblem.Asymmetric>(report.problems.single())
        assertEquals(0.3, problem.difference, 1e-12)
        assertEquals(0.3, report.maxAsymmetry, 1e-12)
        assertNull(report.eigenvalueThreshold)
        // Within tolerance it is accepted and the eigenvalues are reported.
        assertNotNull(KinshipValidation.check(2, tolerance = 0.5) { row, column -> square[row][column] }.minEigenvalue)
    }

    @Test
    fun reportsEveryProblemFoundTogether() = runTest {
        val report = KinshipValidation.check(2) { row, column ->
            if (row == 1 && column == 1) {
                Double.NaN
            } else if (row > column) {
                1.0
            } else {
                2.0
            }
        }
        assertEquals(
            listOf(KinshipProblem.NonFinite::class, KinshipProblem.Asymmetric::class),
            report.problems.map {
                it::class
            },
        )
    }

    @Test
    fun anEmptyMatrixIsValid() = runTest {
        assertTrue(KinshipValidation.check(0) { _, _ -> 0.0 }.isValid)
    }

    @Test
    fun rejectsABadTolerance() = runTest {
        assertFailsWith<IllegalArgumentException> {
            KinshipValidation.check(matrix(floatArrayOf(1f)), tolerance = -1.0)
        }
    }

    @Test
    fun finishesProgressAtTheTotal() = runTest {
        val reports = mutableListOf<Progress>()
        KinshipValidation.check(matrix(floatArrayOf(2f), floatArrayOf(1f, 2f)), progress = { reports += it })
        assertEquals(reports.last().total, reports.last().completed)
    }
}
