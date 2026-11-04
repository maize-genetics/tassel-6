package net.maizegenetics.tassel.math

import kotlin.math.PI
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.sqrt
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class SymmetricEigenDecompositionTest {
    @Test
    fun oneByOneIsTheElementItself() {
        val eigen = SymmetricEigenDecomposition.of(arrayOf(doubleArrayOf(-3.5)))
        assertEquals(1, eigen.size)
        assertEquals(-3.5, eigen.eigenvalue(0), 0.0)
        assertEquals(1.0, eigen[0, 0], 0.0)
    }

    @Test
    fun diagonalMatrixGivesSortedDiagonalAndBasisVectors() {
        val diagonal = doubleArrayOf(2.0, -7.0, 5.0, 0.0)
        val eigen = SymmetricEigenDecomposition.of(4) { row, column -> if (row == column) diagonal[row] else 0.0 }

        assertArrayEquals(doubleArrayOf(5.0, 2.0, 0.0, -7.0), eigen.eigenvalues(), 0.0)
        // Columns are the basis vectors for the sorted eigenvalues, made positive by the sign convention.
        val expectedColumns = intArrayOf(2, 0, 3, 1)
        for (column in 0 until 4) {
            for (row in 0 until 4) {
                assertEquals(if (row == expectedColumns[column]) 1.0 else 0.0, eigen[row, column], 0.0)
            }
        }
    }

    @Test
    fun twoByTwoMatchesClosedFormIncludingSigns() {
        val eigen = SymmetricEigenDecomposition.of(arrayOf(doubleArrayOf(2.0, 1.0), doubleArrayOf(1.0, 2.0)))
        val half = sqrt(0.5)

        assertEquals(3.0, eigen.eigenvalue(0), 1e-15)
        assertEquals(1.0, eigen.eigenvalue(1), 1e-15)
        // Both columns tie in absolute value, so the sign convention makes row 0 positive.
        assertArrayEquals(doubleArrayOf(half, half), eigen.eigenvector(0), 1e-15)
        assertArrayEquals(doubleArrayOf(half, -half), eigen.eigenvector(1), 1e-15)
    }

    @Test
    fun eigenvectorEntryWithLargestMagnitudeIsPositive() {
        val eigen = SymmetricEigenDecomposition.of(gramMatrix(randomMatrix(9, 9, Random(7))))
        for (column in 0 until eigen.size) {
            val vector = eigen.eigenvector(column)
            var pivot = 0
            for (row in vector.indices) if (abs(vector[row]) > abs(vector[pivot])) pivot = row
            assertTrue(vector[pivot] > 0.0, "column $column has a negative leading entry ${vector[pivot]}")
        }
    }

    @Test
    fun matchesClosedFormOfTheToeplitzTridiagonalMatrix() {
        // The n x n matrix with 2 on the diagonal and -1 next to it has eigenvalues 4 sin^2(k pi / (2n + 2)).
        val n = 20
        val eigen =
            SymmetricEigenDecomposition.of(n) { row, column ->
                when (row - column) {
                    0 -> 2.0
                    1, -1 -> -1.0
                    else -> 0.0
                }
            }
        val expected = (1..n).map { 2.0 - 2.0 * cos(it * PI / (n + 1)) }.sortedDescending()
        expected.forEachIndexed { index, value ->
            assertEquals(value, eigen.eigenvalue(index), 1e-13)
        }
        assertOrthonormal(eigen)
    }

    @Test
    fun repeatedEigenvaluesStayOrthonormal() {
        // The all-ones matrix has eigenvalues n, 0, 0, ... so all but the first column are a degenerate subspace.
        val n = 7
        val eigen = SymmetricEigenDecomposition.of(n) { _, _ -> 1.0 }
        assertEquals(n.toDouble(), eigen.eigenvalue(0), 1e-13)
        for (index in 1 until n) assertEquals(0.0, eigen.eigenvalue(index), 1e-13)
        assertOrthonormal(eigen)
        assertReconstructs(eigen) { _, _ -> 1.0 }
    }

    @Test
    fun reconstructsRandomSymmetricMatrices() {
        for (n in intArrayOf(2, 3, 17, 64, 120)) {
            val matrix = symmetricMatrix(n, Random(n.toLong() * 31))
            val eigen = SymmetricEigenDecomposition.of(matrix)
            assertOrthonormal(eigen)
            assertReconstructs(eigen) { row, column -> matrix[row][column] }
            // The trace is the sum of the eigenvalues.
            val trace = (0 until n).sumOf { matrix[it][it] }
            assertEquals(trace, eigen.eigenvalues().sum(), 1e-10 * n)
        }
    }

    @Test
    fun reconstructsAGramMatrixWithNearZeroEigenvalues() {
        // 40 observations of 12 variables: rank 12, so 28 eigenvalues are zero. This is the shape #19 decomposes.
        val gram = gramMatrix(randomMatrix(40, 12, Random(2026)))
        val eigen = SymmetricEigenDecomposition.of(gram)
        assertOrthonormal(eigen)
        assertReconstructs(eigen) { row, column -> gram[row][column] }
        for (index in 12 until 40) assertEquals(0.0, eigen.eigenvalue(index), 1e-11)
        for (index in 0 until 40) {
            assertTrue(
                eigen.eigenvalue(index) <= eigen.eigenvalue(maxOf(index - 1, 0)) + 1e-15,
                "eigenvalues are not descending at $index",
            )
        }
    }

    @Test
    fun ignoresTheUpperTriangle() {
        val n = 15
        val matrix = symmetricMatrix(n, Random(99))
        val expected = SymmetricEigenDecomposition.of(matrix)

        val vandalized = Array(n) { row -> DoubleArray(n) { column -> if (column > row) 1e6 else matrix[row][column] } }
        val actual = SymmetricEigenDecomposition.of(vandalized)

        assertArrayEquals(expected.eigenvalues(), actual.eigenvalues(), 0.0)
        for (column in 0 until n) assertArrayEquals(expected.eigenvector(column), actual.eigenvector(column), 0.0)
    }

    @Test
    fun theThreeFactoriesAgree() {
        val n = 11
        val matrix = symmetricMatrix(n, Random(5))
        val flat = DoubleArray(n * n) { matrix[it / n][it % n] }

        val fromRows = SymmetricEigenDecomposition.of(matrix)
        val fromFlat = SymmetricEigenDecomposition.of(n, flat)
        val fromLambda = SymmetricEigenDecomposition.of(n) { row, column -> matrix[row][column] }

        assertArrayEquals(fromRows.eigenvalues(), fromFlat.eigenvalues(), 0.0)
        assertArrayEquals(fromRows.eigenvalues(), fromLambda.eigenvalues(), 0.0)
        for (column in 0 until n) {
            assertArrayEquals(fromRows.eigenvector(column), fromFlat.eigenvector(column), 0.0)
            assertArrayEquals(fromRows.eigenvector(column), fromLambda.eigenvector(column), 0.0)
        }
    }

    @Test
    fun rejectsBadInput() {
        assertFailsWith<IllegalArgumentException> { SymmetricEigenDecomposition.of(emptyArray<DoubleArray>()) }
        assertFailsWith<IllegalArgumentException> {
            SymmetricEigenDecomposition.of(arrayOf(doubleArrayOf(1.0, 2.0)))
        }
        assertFailsWith<IllegalArgumentException> { SymmetricEigenDecomposition.of(2, DoubleArray(3)) }
        assertFailsWith<IllegalArgumentException> { SymmetricEigenDecomposition.of(0) { _, _ -> 1.0 } }
        assertFailsWith<IllegalArgumentException> {
            SymmetricEigenDecomposition.of(2) { row, column -> if (row == column) Double.NaN else 1.0 }
        }
        val eigen = SymmetricEigenDecomposition.of(arrayOf(doubleArrayOf(1.0)))
        assertFailsWith<IllegalArgumentException> { eigen.eigenvalue(1) }
        assertFailsWith<IllegalArgumentException> { eigen.eigenvector(-1) }
        assertFailsWith<IllegalArgumentException> { eigen[0, 1] }
    }

    private fun assertOrthonormal(eigen: SymmetricEigenDecomposition) {
        val n = eigen.size
        for (i in 0 until n) {
            for (j in i until n) {
                var dot = 0.0
                for (row in 0 until n) dot += eigen[row, i] * eigen[row, j]
                assertEquals(if (i == j) 1.0 else 0.0, dot, 1e-12, "V'V is not the identity at ($i, $j)")
            }
        }
    }

    /** Checks `A v = lambda v` column by column, scaled by the size of the matrix. */
    private fun assertReconstructs(eigen: SymmetricEigenDecomposition, element: (row: Int, column: Int) -> Double) {
        val n = eigen.size
        var scale = 0.0
        for (row in 0 until n) {
            for (column in 0 until n) scale = maxOf(scale, abs(element(row, column)))
        }
        val tolerance = 1e-12 * n * maxOf(scale, 1.0)
        for (column in 0 until n) {
            val lambda = eigen.eigenvalue(column)
            for (row in 0 until n) {
                var product = 0.0
                for (k in 0 until n) product += element(row, k) * eigen[k, column]
                assertEquals(lambda * eigen[row, column], product, tolerance, "A v != lambda v at ($row, $column)")
            }
        }
    }

    private fun assertArrayEquals(expected: DoubleArray, actual: DoubleArray, tolerance: Double) {
        assertEquals(expected.size, actual.size, "sizes differ")
        expected.indices.forEach { assertEquals(expected[it], actual[it], tolerance, "differ at $it") }
    }

    private fun symmetricMatrix(n: Int, random: Random): Array<DoubleArray> {
        val matrix = Array(n) { DoubleArray(n) }
        for (row in 0 until n) {
            for (column in 0..row) {
                val value = random.nextDouble(-5.0, 5.0)
                matrix[row][column] = value
                matrix[column][row] = value
            }
        }
        return matrix
    }

    private fun randomMatrix(rows: Int, columns: Int, random: Random): Array<DoubleArray> =
        Array(rows) { DoubleArray(columns) { random.nextDouble(-2.0, 2.0) } }

    /** `matrix * matrix'`, a positive-semidefinite matrix of side `matrix.size`. */
    private fun gramMatrix(matrix: Array<DoubleArray>): Array<DoubleArray> {
        val n = matrix.size
        val gram = Array(n) { DoubleArray(n) }
        for (row in 0 until n) {
            for (column in 0..row) {
                var sum = 0.0
                for (k in matrix[row].indices) sum += matrix[row][k] * matrix[column][k]
                gram[row][column] = sum
                gram[column][row] = sum
            }
        }
        return gram
    }
}
