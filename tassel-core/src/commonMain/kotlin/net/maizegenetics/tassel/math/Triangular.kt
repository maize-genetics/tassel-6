package net.maizegenetics.tassel.math

/**
 * In-place solves with a triangular matrix stored in the leading `n x n` block of a [DenseMatrix]. Only the named
 * triangle is read, so a QR factor (R above the diagonal, reflectors below) or a Cholesky workspace can be passed as
 * it is. Each solve overwrites `b[0 until n]` with the solution.
 *
 * A zero diagonal entry divides by zero and produces infinities or `NaN`; callers decide rank first (see
 * [HouseholderQr.rank]) and pass only the leading block they trust. Level-1 work goes through [denseKernels].
 */
object Triangular {
    /** Solves `R x = b` for upper-triangular `R` (back substitution, one column at a time). */
    fun solveUpper(r: DenseMatrix, n: Int, b: DoubleArray) {
        check(r, n, b)
        val kernels = denseKernels
        for (j in n - 1 downTo 0) {
            b[j] /= r.data[j * r.rows + j]
            kernels.axpy(j, -b[j], r.data, j * r.rows, b, 0)
        }
    }

    /** Solves `R' x = b` for upper-triangular `R` (forward substitution on the transpose). */
    fun solveUpperTransposed(r: DenseMatrix, n: Int, b: DoubleArray) {
        check(r, n, b)
        val kernels = denseKernels
        for (j in 0 until n) {
            b[j] = (b[j] - kernels.dot(j, r.data, j * r.rows, b, 0)) / r.data[j * r.rows + j]
        }
    }

    /** Solves `L x = b` for lower-triangular `L` (forward substitution, one column at a time). */
    fun solveLower(l: DenseMatrix, n: Int, b: DoubleArray) {
        check(l, n, b)
        val kernels = denseKernels
        for (j in 0 until n) {
            b[j] /= l.data[j * l.rows + j]
            kernels.axpy(n - j - 1, -b[j], l.data, j * l.rows + j + 1, b, j + 1)
        }
    }

    /** Solves `L' x = b` for lower-triangular `L` (back substitution on the transpose). */
    fun solveLowerTransposed(l: DenseMatrix, n: Int, b: DoubleArray) {
        check(l, n, b)
        val kernels = denseKernels
        for (j in n - 1 downTo 0) {
            b[j] = (b[j] - kernels.dot(n - j - 1, l.data, j * l.rows + j + 1, b, j + 1)) / l.data[j * l.rows + j]
        }
    }

    private fun check(matrix: DenseMatrix, n: Int, b: DoubleArray) {
        require(n >= 0 && n <= matrix.rows && n <= matrix.columns) {
            "a $n x $n block does not fit in a ${matrix.rows} x ${matrix.columns} matrix"
        }
        require(b.size >= n) { "b has ${b.size} values, expected at least $n" }
    }
}
