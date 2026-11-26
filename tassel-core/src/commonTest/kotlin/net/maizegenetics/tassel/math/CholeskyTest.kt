package net.maizegenetics.tassel.math

import kotlin.math.abs
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertIs
import kotlin.test.assertTrue

class CholeskyTest {
    @Test
    fun factorsAndSolvesASymmetricPositiveDefiniteSystem() {
        val a = DenseMatrix.ofRows(
            doubleArrayOf(4.0, 12.0, -16.0),
            doubleArrayOf(12.0, 37.0, -43.0),
            doubleArrayOf(-16.0, -43.0, 98.0),
        )
        val factor = assertIs<CholeskyResult.Success>(Cholesky.of(a)).factor
        val lower = factor.lower()
        // The textbook factor of this matrix.
        val expected = DenseMatrix.ofRows(
            doubleArrayOf(2.0, 0.0, 0.0),
            doubleArrayOf(6.0, 1.0, 0.0),
            doubleArrayOf(-8.0, 5.0, 3.0),
        )
        for (i in 0 until 3) for (j in 0 until 3) assertEquals(expected[i, j], lower[i, j], "L[$i,$j]")
        val x = doubleArrayOf(1.0, -2.0, 0.5)
        val b = DoubleArray(3) { i -> (0 until 3).sumOf { a[i, it] * x[it] } }
        val solved = factor.solve(b)
        for (i in 0 until 3) assertTrue(abs(solved[i] - x[i]) < 1e-12, "x[$i] = ${solved[i]}")
    }

    @Test
    fun readsOnlyTheLowerTriangle() {
        val clean = DenseMatrix.ofRows(doubleArrayOf(4.0, 2.0), doubleArrayOf(2.0, 3.0))
        val dirty = clean.copy().also { it[0, 1] = 1e9 }
        val a = assertIs<CholeskyResult.Success>(Cholesky.of(clean)).factor.lower()
        val b = assertIs<CholeskyResult.Success>(Cholesky.of(dirty)).factor.lower()
        assertEquals(a.data.toList(), b.data.toList())
    }

    @Test
    fun reportsTheFailingPivotAsAValue() {
        // The second pivot is 1 - 2 * 2 / 1 = -3.
        val indefinite = DenseMatrix.ofRows(doubleArrayOf(1.0, 2.0), doubleArrayOf(2.0, 1.0))
        assertEquals(1, assertIs<CholeskyResult.NotPositiveDefinite>(Cholesky.of(indefinite)).pivot)
        assertEquals(0, assertIs<CholeskyResult.NotPositiveDefinite>(Cholesky.of(DenseMatrix(1, 1))).pivot)
        // A positive semi-definite matrix: the rank-one matrix [1 1; 1 1] has a zero second pivot.
        val semi = DenseMatrix.ofRows(doubleArrayOf(1.0, 1.0), doubleArrayOf(1.0, 1.0))
        assertEquals(1, assertIs<CholeskyResult.NotPositiveDefinite>(Cholesky.of(semi)).pivot)
        val nan = DenseMatrix.ofRows(doubleArrayOf(1.0, 0.0), doubleArrayOf(0.0, Double.NaN))
        assertEquals(1, assertIs<CholeskyResult.NotPositiveDefinite>(Cholesky.of(nan)).pivot)
        // Failure at a later pivot of a larger matrix.
        val random = Random(1)
        val big = DenseMatrix(6, 6)
        for (i in 0 until 6) big[i, i] = 1.0 + random.nextDouble()
        big[4, 4] = -1.0
        assertEquals(4, assertIs<CholeskyResult.NotPositiveDefinite>(Cholesky.of(big)).pivot)
    }

    @Test
    fun agreesWithQrOnAGramMatrix() {
        val random = Random(2)
        val x = DenseMatrix(30, 5, DoubleArray(150) { random.nextDouble(-1.0, 1.0) })
        val y = DoubleArray(30) { random.nextDouble(-1.0, 1.0) }
        val gram = DenseMatrix(5, 5)
        val rhs = DoubleArray(5)
        for (i in 0 until 5) {
            for (j in 0 until 5) gram[i, j] = (0 until 30).sumOf { x[it, i] * x[it, j] }
            rhs[i] = (0 until 30).sumOf { x[it, i] * y[it] }
        }
        val viaCholesky = assertIs<CholeskyResult.Success>(Cholesky.of(gram)).factor.solve(rhs)
        val viaQr = HouseholderQr.of(x).solve(y)
        for (j in 0 until 5) assertTrue(abs(viaCholesky[j] - viaQr[j]) < 1e-10)
    }

    @Test
    fun rejectsNonSquareInputAndWrongRightHandSide() {
        assertFailsWith<IllegalArgumentException> { Cholesky.of(DenseMatrix(2, 3)) }
        val factor = assertIs<CholeskyResult.Success>(Cholesky.of(DenseMatrix.identity(3))).factor
        assertFailsWith<IllegalArgumentException> { factor.solve(DoubleArray(2)) }
    }
}
