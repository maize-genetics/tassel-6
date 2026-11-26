package net.maizegenetics.tassel.math

import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.coroutineScope
import net.maizegenetics.tassel.distance.kinshipWorkerCount
import net.maizegenetics.tassel.distance.pairBands

/**
 * BLAS-style level-1, level-2, and level-3 kernels on column-major [DenseMatrix] storage.
 *
 * Everything the decompositions and the model layer multiply goes through [denseKernels], so a later platform
 * backend (native BLAS, FFM, a GPU) can replace the product kernels and leave callers alone, the way
 * `kinshipTriangleUpdate` does for kinship (ADR-0020, ADR-0025). The default is [PureKotlinDenseKernels].
 *
 * An implementation must not change results with the number of workers: the order in which each output element is
 * accumulated is fixed by the arguments alone, and workers only decide which elements a thread computes. Different
 * backends may of course round differently from each other.
 *
 * This is an `interface` rather than the `fun interface` of the kinship hook because it has five operations; a
 * `fun interface` allows exactly one abstract member.
 *
 * Conventions follow BLAS: `beta = 0` overwrites the destination without reading it, so a `NaN` or garbage left in a
 * reused workspace cannot leak into the result; sizes are validated and a mismatch throws
 * [IllegalArgumentException].
 */
internal interface DenseKernels {
    /** `sum_i x[xOffset + i] * y[yOffset + i]` for `i` in `0 until n`. */
    fun dot(n: Int, x: DoubleArray, xOffset: Int, y: DoubleArray, yOffset: Int): Double

    /** `y[yOffset + i] += alpha * x[xOffset + i]` for `i` in `0 until n`. */
    fun axpy(n: Int, alpha: Double, x: DoubleArray, xOffset: Int, y: DoubleArray, yOffset: Int)

    /**
     * `y = alpha * op(a) * x + beta * y`, where `op(a)` is `a'` when [transpose] is true. [x] has `op(a).columns`
     * values and [y] has `op(a).rows` values.
     */
    fun gemv(transpose: Boolean, a: DenseMatrix, alpha: Double, x: DoubleArray, beta: Double, y: DoubleArray)

    /**
     * `c = alpha * op(a) * op(b) + beta * c`, where [transposeA] and [transposeB] select `op`. [c] is
     * `op(a).rows x op(b).columns`.
     */
    suspend fun gemm(
        transposeA: Boolean,
        transposeB: Boolean,
        alpha: Double,
        a: DenseMatrix,
        b: DenseMatrix,
        beta: Double,
        c: DenseMatrix,
    )

    /**
     * Symmetric rank-k update of the lower triangle of [c], which is `n x n`: `c = alpha * a * a' + beta * c` when
     * [transpose] is false (`a` is `n x k`), `c = alpha * a' * a + beta * c` when it is true (`a` is `k x n`). The
     * strict upper triangle of [c] is neither read nor written.
     */
    suspend fun syrk(transpose: Boolean, alpha: Double, a: DenseMatrix, beta: Double, c: DenseMatrix)
}

/** The kernels every dense routine calls. Off-by-default spikes and later backends assign a different one. */
internal var denseKernels: DenseKernels = PureKotlinDenseKernels()

/**
 * Plain Kotlin kernels. [gemm] and [syrk] split the columns of the result across up to `workers()` coroutines; each
 * column is computed by one coroutine in a fixed order, so the bits do not depend on the split. The coroutines
 * inherit the caller's dispatcher: the JVM and native engines run them on `Dispatchers.Default`, a test dispatcher
 * and the browser's single thread run them one after another.
 */
internal class PureKotlinDenseKernels(private val workers: () -> Int = ::kinshipWorkerCount) : DenseKernels {
    override fun dot(n: Int, x: DoubleArray, xOffset: Int, y: DoubleArray, yOffset: Int): Double {
        var sum = 0.0
        for (i in 0 until n) sum += x[xOffset + i] * y[yOffset + i]
        return sum
    }

    override fun axpy(n: Int, alpha: Double, x: DoubleArray, xOffset: Int, y: DoubleArray, yOffset: Int) {
        if (alpha == 0.0) return
        for (i in 0 until n) y[yOffset + i] += alpha * x[xOffset + i]
    }

    override fun gemv(transpose: Boolean, a: DenseMatrix, alpha: Double, x: DoubleArray, beta: Double, y: DoubleArray) {
        val outputs = if (transpose) a.columns else a.rows
        val inputs = if (transpose) a.rows else a.columns
        require(x.size == inputs) { "x has ${x.size} values, expected $inputs" }
        require(y.size == outputs) { "y has ${y.size} values, expected $outputs" }
        scale(y, 0, outputs, beta)
        if (transpose) {
            for (j in 0 until a.columns) y[j] += alpha * dot(a.rows, a.data, j * a.rows, x, 0)
        } else {
            for (j in 0 until a.columns) axpy(a.rows, alpha * x[j], a.data, j * a.rows, y, 0)
        }
    }

    override suspend fun gemm(
        transposeA: Boolean,
        transposeB: Boolean,
        alpha: Double,
        a: DenseMatrix,
        b: DenseMatrix,
        beta: Double,
        c: DenseMatrix,
    ) {
        val m = if (transposeA) a.columns else a.rows
        val k = if (transposeA) a.rows else a.columns
        val kb = if (transposeB) b.columns else b.rows
        val n = if (transposeB) b.rows else b.columns
        require(k == kb) { "inner dimensions differ: op(a) has $k columns, op(b) has $kb rows" }
        require(c.rows == m && c.columns == n) { "c is ${c.rows} x ${c.columns}, expected $m x $n" }
        val bands = evenBands(n, workers())
        runBands(bands) { columns ->
            val gathered = if (transposeA) DoubleArray(k) else DoubleArray(0)
            for (j in columns) {
                val out = j * m
                scale(c.data, out, m, beta)
                if (transposeA) {
                    // c[i, j] += alpha * a[:, i] . op(b)[:, j]; gather op(b)[:, j] so it is contiguous.
                    for (l in 0 until k) {
                        gathered[l] =
                            if (transposeB) b.data[l * b.rows + j] else b.data[j * b.rows + l]
                    }
                    for (i in 0 until m) c.data[out + i] += alpha * dot(k, a.data, i * a.rows, gathered, 0)
                } else {
                    for (l in 0 until k) {
                        val scalar = alpha * (if (transposeB) b.data[l * b.rows + j] else b.data[j * b.rows + l])
                        axpy(m, scalar, a.data, l * a.rows, c.data, out)
                    }
                }
            }
        }
    }

    override suspend fun syrk(transpose: Boolean, alpha: Double, a: DenseMatrix, beta: Double, c: DenseMatrix) {
        val n = if (transpose) a.columns else a.rows
        val k = if (transpose) a.rows else a.columns
        require(c.rows == n && c.columns == n) { "c is ${c.rows} x ${c.columns}, expected $n x $n" }
        if (n == 0) return
        // Column j of the lower triangle holds n - j values, the same shape as row j of the kinship triangle.
        runBands(pairBands(n, workers())) { columns ->
            for (j in columns) {
                val out = j * n
                scale(c.data, out + j, n - j, beta)
                if (transpose) {
                    for (i in j until n) c.data[out + i] += alpha * dot(k, a.data, i * a.rows, a.data, j * a.rows)
                } else {
                    for (l in 0 until k) axpy(n - j, alpha * a.data[l * n + j], a.data, l * n + j, c.data, out + j)
                }
            }
        }
    }

    private fun scale(array: DoubleArray, offset: Int, length: Int, factor: Double) {
        if (factor == 1.0) return
        if (factor == 0.0) {
            array.fill(0.0, offset, offset + length)
        } else {
            for (i in offset until offset + length) array[i] *= factor
        }
    }

    private suspend fun runBands(bands: List<IntRange>, body: (IntRange) -> Unit) {
        if (bands.size <= 1) {
            for (band in bands) body(band)
            return
        }
        coroutineScope { bands.map { band -> async { body(band) } }.awaitAll() }
    }
}

/** Splits `0 until count` into at most [workers] contiguous bands of nearly equal size, in order. */
internal fun evenBands(count: Int, workers: Int): List<IntRange> {
    if (count <= 0) return emptyList()
    val bands = workers.coerceIn(1, count)
    return List(bands) { index ->
        (index.toLong() * count / bands).toInt() until
            ((index + 1).toLong() * count / bands).toInt()
    }
}
