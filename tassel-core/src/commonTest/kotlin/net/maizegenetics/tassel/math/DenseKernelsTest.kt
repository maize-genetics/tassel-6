package net.maizegenetics.tassel.math

import kotlinx.coroutines.test.runTest
import kotlin.math.abs
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/** The pure-Kotlin kernels against naive triple loops, and their independence from the worker count. */
class DenseKernelsTest {
    private val kernels = PureKotlinDenseKernels { 1 }

    @Test
    fun dotAndAxpyHonourOffsets() {
        val x = doubleArrayOf(9.0, 1.0, 2.0, 3.0)
        val y = doubleArrayOf(4.0, 5.0, 6.0, 9.0)
        assertEquals(32.0, kernels.dot(3, x, 1, y, 0))
        kernels.axpy(3, 2.0, x, 1, y, 0)
        assertContentEquals(doubleArrayOf(6.0, 9.0, 12.0, 9.0), y)
        assertEquals(0.0, kernels.dot(0, x, 0, y, 0))
    }

    @Test
    fun gemvMatchesNaive() {
        val a = random(5, 3, 1)
        val x = Random(2).let { r -> DoubleArray(3) { r.nextDouble(-1.0, 1.0) } }
        val z = Random(3).let { r -> DoubleArray(5) { r.nextDouble(-1.0, 1.0) } }
        val y = DoubleArray(5) { 0.5 }
        kernels.gemv(false, a, 2.0, x, 3.0, y)
        for (i in 0 until 5) {
            val expected = 3.0 * 0.5 + 2.0 * (0 until 3).sumOf { a[i, it] * x[it] }
            assertClose(expected, y[i])
        }
        val t = DoubleArray(3) { 1.5 }
        kernels.gemv(true, a, -1.0, z, 1.0, t)
        for (j in 0 until 3) assertClose(1.5 - (0 until 5).sumOf { a[it, j] * z[it] }, t[j])
    }

    @Test
    fun betaZeroOverwritesGarbage() = runTest {
        val a = random(3, 2, 4)
        val y = DoubleArray(3) { Double.NaN }
        kernels.gemv(false, a, 1.0, doubleArrayOf(1.0, 1.0), 0.0, y)
        assertTrue(y.all { it.isFinite() })
        val c = DenseMatrix(3, 3, DoubleArray(9) { Double.NaN })
        kernels.gemm(false, true, 1.0, a, a, 0.0, c)
        assertTrue(c.data.all { it.isFinite() })
        val s = DenseMatrix(3, 3, DoubleArray(9) { Double.NaN })
        kernels.syrk(false, 1.0, a, 0.0, s)
        for (j in 0 until 3) for (i in j until 3) assertTrue(s[i, j].isFinite())
    }

    @Test
    fun gemmMatchesNaiveForEveryTransposeCombination() = runTest {
        val m = 6
        val k = 4
        val n = 5
        for (transposeA in listOf(false, true)) {
            for (transposeB in listOf(false, true)) {
                val a = if (transposeA) random(k, m, 5) else random(m, k, 5)
                val b = if (transposeB) random(n, k, 6) else random(k, n, 6)
                val c = random(m, n, 7)
                val before = c.copy()
                kernels.gemm(transposeA, transposeB, 1.5, a, b, -0.5, c)
                for (i in 0 until m) {
                    for (j in 0 until n) {
                        var sum = 0.0
                        for (l in 0 until k) {
                            sum += (if (transposeA) a[l, i] else a[i, l]) * (if (transposeB) b[j, l] else b[l, j])
                        }
                        assertClose(-0.5 * before[i, j] + 1.5 * sum, c[i, j], "A'=$transposeA B'=$transposeB ($i,$j)")
                    }
                }
            }
        }
    }

    @Test
    fun syrkFillsOnlyTheLowerTriangle() = runTest {
        val n = 5
        for (transpose in listOf(false, true)) {
            val a = if (transpose) random(7, n, 8) else random(n, 7, 8)
            val c = random(n, n, 9)
            val before = c.copy()
            kernels.syrk(transpose, 2.0, a, 0.5, c)
            for (j in 0 until n) {
                for (i in 0 until n) {
                    if (i < j) {
                        assertEquals(before[i, j], c[i, j], "upper triangle ($i,$j) must be untouched")
                    } else {
                        val dot = (0 until 7).sumOf { l ->
                            if (transpose) a[l, i] * a[l, j] else a[i, l] * a[j, l]
                        }
                        assertClose(0.5 * before[i, j] + 2.0 * dot, c[i, j], "transpose=$transpose ($i,$j)")
                    }
                }
            }
        }
    }

    @Test
    fun gemmAndSyrkAreBitIdenticalForAnyNumberOfWorkers() = runTest {
        val a = random(37, 23, 10)
        val b = random(23, 41, 11)
        val reference = PureKotlinDenseKernels { 1 }
        val gemmExpected = DenseMatrix(37, 41).also { reference.gemm(false, false, 1.0, a, b, 0.0, it) }
        val gemmTExpected = DenseMatrix(23, 23).also { reference.gemm(true, false, 1.0, a, a, 0.0, it) }
        val syrkExpected = DenseMatrix(37, 37).also { reference.syrk(false, 1.0, a, 0.0, it) }
        val syrkTExpected = DenseMatrix(23, 23).also { reference.syrk(true, 1.0, a, 0.0, it) }
        for (workers in listOf(2, 3, 4, 7, 64)) {
            val parallel = PureKotlinDenseKernels { workers }
            val g = DenseMatrix(37, 41).also { parallel.gemm(false, false, 1.0, a, b, 0.0, it) }
            val gt = DenseMatrix(23, 23).also { parallel.gemm(true, false, 1.0, a, a, 0.0, it) }
            val s = DenseMatrix(37, 37).also { parallel.syrk(false, 1.0, a, 0.0, it) }
            val st = DenseMatrix(23, 23).also { parallel.syrk(true, 1.0, a, 0.0, it) }
            assertContentEquals(gemmExpected.data, g.data, "gemm workers=$workers")
            assertContentEquals(gemmTExpected.data, gt.data, "gemm A'B workers=$workers")
            assertContentEquals(syrkExpected.data, s.data, "syrk workers=$workers")
            assertContentEquals(syrkTExpected.data, st.data, "syrk' workers=$workers")
        }
    }

    @Test
    fun evenBandsCoverEveryIndexOnceInOrder() {
        for (count in listOf(0, 1, 5, 10, 101)) {
            for (workers in listOf(1, 2, 3, 8, 200)) {
                val bands = evenBands(count, workers)
                assertEquals((0 until count).toList(), bands.flatMap { it.toList() })
                assertTrue(bands.size <= maxOf(1, workers))
            }
        }
    }

    @Test
    fun rejectsMismatchedSizes() = runTest {
        val a = DenseMatrix(3, 2)
        assertFailsWith<IllegalArgumentException> { kernels.gemv(false, a, 1.0, DoubleArray(3), 0.0, DoubleArray(3)) }
        assertFailsWith<IllegalArgumentException> {
            kernels.gemm(false, false, 1.0, a, DenseMatrix(3, 2), 0.0, DenseMatrix(3, 2))
        }
        assertFailsWith<IllegalArgumentException> { kernels.syrk(false, 1.0, a, 0.0, DenseMatrix(2, 2)) }
    }

    @Test
    fun theDefaultHookIsThePureKotlinImplementation() {
        assertTrue(denseKernels is PureKotlinDenseKernels)
    }

    private fun random(rows: Int, columns: Int, seed: Int): DenseMatrix {
        val random = Random(seed)
        return DenseMatrix(rows, columns, DoubleArray(rows * columns) { random.nextDouble(-1.0, 1.0) })
    }

    private fun assertClose(expected: Double, actual: Double, message: String = "") {
        assertTrue(
            abs(expected - actual) <= 1e-12 * maxOf(1.0, abs(expected)),
            "$message expected $expected got $actual",
        )
    }
}
