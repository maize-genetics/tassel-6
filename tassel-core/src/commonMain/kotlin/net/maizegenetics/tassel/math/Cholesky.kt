package net.maizegenetics.tassel.math

import kotlin.math.sqrt

/** The outcome of [Cholesky.of]: a factor, or the pivot where the matrix stopped being positive definite. */
sealed interface CholeskyResult {
    /** `A = L L'` succeeded; [factor] solves systems with `A`. */
    class Success(val factor: Cholesky) : CholeskyResult

    /**
     * The leading `(pivot + 1) x (pivot + 1)` block is not positive definite: the pivot `pivot` (zero-based) was zero,
     * negative, or not a number once the earlier columns were eliminated. Reported as a value so a caller scanning
     * thousands of matrices can count and skip failures without exception handling in a loop.
     */
    class NotPositiveDefinite(val pivot: Int) : CholeskyResult
}

/**
 * The Cholesky factor `L` of a symmetric positive-definite matrix, `A = L L'`, with `L` lower triangular.
 *
 * Only the lower triangle of the input, diagonal included, is read, the same convention as
 * [SymmetricEigenDecomposition]. The factorization is left-looking and goes down columns, so its level-1 work is the
 * [denseKernels] `axpy`. No inverse is formed; [solve] does two triangular solves.
 */
class Cholesky private constructor(
    /** Order of the factored matrix. */
    val size: Int,
    private val lowerFactor: DenseMatrix,
) {
    /** A copy of `L`; entries above the diagonal are zero. */
    fun lower(): DenseMatrix = lowerFactor.copy()

    /** Solves `A x = b`, returning a new array. */
    fun solve(b: DoubleArray): DoubleArray = b.copyOf().also { solveInPlace(it) }

    /** Solves `A x = b`, overwriting [b]. */
    fun solveInPlace(b: DoubleArray) {
        require(b.size == size) { "b has ${b.size} values, expected $size" }
        Triangular.solveLower(lowerFactor, size, b)
        Triangular.solveLowerTransposed(lowerFactor, size, b)
    }

    companion object {
        /** Factors [matrix], which must be square, reading its lower triangle. */
        fun of(matrix: DenseMatrix): CholeskyResult {
            require(matrix.rows == matrix.columns) { "matrix must be square, got ${matrix.rows} x ${matrix.columns}" }
            val n = matrix.rows
            val kernels = denseKernels
            val l = DenseMatrix(n, n)
            for (j in 0 until n) {
                // Copy the lower part of column j, then subtract the contribution of every earlier column.
                matrix.data.copyInto(l.data, j * n + j, j * n + j, j * n + n)
                for (k in 0 until j) {
                    kernels.axpy(n - j, -l.data[k * n + j], l.data, k * n + j, l.data, j * n + j)
                }
                val pivot = l.data[j * n + j]
                // `!(pivot > 0)` also catches NaN.
                if (!(pivot > 0.0) || !pivot.isFinite()) return CholeskyResult.NotPositiveDefinite(j)
                val root = sqrt(pivot)
                l.data[j * n + j] = root
                for (i in j + 1 until n) l.data[j * n + i] /= root
            }
            return CholeskyResult.Success(Cholesky(n, l))
        }
    }
}
