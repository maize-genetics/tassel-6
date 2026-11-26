package net.maizegenetics.tassel.math

import kotlinx.coroutines.test.runTest
import kotlinx.coroutines.yield
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/** The progress-reporting overloads of [SymmetricEigenDecomposition.of] run the same code as the blocking ones. */
class SymmetricEigenProgressTest {
    private fun symmetric(size: Int, seed: Int): Array<DoubleArray> {
        val random = Random(seed)
        val matrix = Array(size) { DoubleArray(size) }
        for (i in 0 until size) {
            for (j in 0..i) {
                val value = random.nextDouble(-1.0, 1.0)
                matrix[i][j] = value
                matrix[j][i] = value
            }
        }
        return matrix
    }

    private fun assertIdentical(expected: SymmetricEigenDecomposition, actual: SymmetricEigenDecomposition) {
        assertEquals(expected.size, actual.size)
        assertContentEquals(expected.eigenvalues(), actual.eigenvalues())
        for (index in 0 until expected.size) assertContentEquals(expected.eigenvector(index), actual.eigenvector(index))
    }

    @Test
    fun resultsAreBitIdenticalToTheBlockingVersion() = runTest {
        for (size in listOf(1, 2, 3, 17, 40)) {
            val matrix = symmetric(size, size)
            val blocking = SymmetricEigenDecomposition.of(matrix)
            assertIdentical(blocking, SymmetricEigenDecomposition.of(matrix) { })
            val flat = DoubleArray(size * size) { matrix[it / size][it % size] }
            assertIdentical(blocking, SymmetricEigenDecomposition.of(size, flat) { })
            assertIdentical(blocking, SymmetricEigenDecomposition.of(size, { }) { r, c -> matrix[r][c] })
        }
    }

    @Test
    fun progressIsMonotoneStartsAtZeroAndEndsAtOne() = runTest {
        val size = 30
        val reports = mutableListOf<Double>()
        SymmetricEigenDecomposition.of(symmetric(size, 5)) { reports += it }
        assertEquals(0.0, reports.first())
        assertEquals(1.0, reports.last())
        assertTrue(reports.zipWithNext().all { (a, b) -> a <= b }, "progress must not go backwards: $reports")
        assertTrue(reports.all { it in 0.0..1.0 })
        // At least one report per Householder column and one per QL eigenvalue.
        assertTrue(reports.size >= 2 * size, "only ${reports.size} reports")
    }

    @Test
    fun theCallbackIsASuspensionPointThatMayYieldOrFail() = runTest {
        val matrix = symmetric(20, 6)
        var yields = 0
        val result = SymmetricEigenDecomposition.of(matrix) {
            yield()
            yields++
        }
        assertTrue(yields > 0)
        assertIdentical(SymmetricEigenDecomposition.of(matrix), result)
        assertFailsWith<IllegalStateException> {
            SymmetricEigenDecomposition.of(matrix) { fraction ->
                if (fraction >
                    0.6
                ) {
                    throw IllegalStateException("stop")
                }
            }
        }
    }

    @Test
    fun validationMatchesTheBlockingOverloads() = runTest {
        assertFailsWith<IllegalArgumentException> {
            SymmetricEigenDecomposition.of(arrayOf(doubleArrayOf(1.0, 2.0))) { }
        }
        assertFailsWith<IllegalArgumentException> { SymmetricEigenDecomposition.of(2, DoubleArray(3)) { } }
        assertFailsWith<IllegalArgumentException> {
            SymmetricEigenDecomposition.of(2, { }) { _, _ -> Double.NaN }
        }
    }
}
