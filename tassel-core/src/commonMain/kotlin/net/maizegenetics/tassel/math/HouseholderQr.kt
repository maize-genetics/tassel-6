package net.maizegenetics.tassel.math

import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt

/** The output of [HouseholderQr.solveRefinedWithResiduals]: coefficients (aliased columns zero) and `b - A x`. */
class RefinedSolution(val coefficients: DoubleArray, val residuals: DoubleArray)

/**
 * Householder QR factorization with optional column pivoting and rank detection: `A P = Q R` for an `m x n` matrix
 * `A`, a permutation `P`, an orthogonal `Q`, and an upper-trapezoidal `R`.
 *
 * `Q` is never formed. The reflectors live below the diagonal of the factor, LAPACK style (`H = I - tau v v'` with
 * `v[0] = 1` implied), and [applyQt] and [applyQ] apply them to vectors or matrices. Least squares and standard
 * errors are triangular solves against `R`, never an inverse of `A'A` (ADR-0025). The pivoting follows LAPACK
 * `dgeqp3`/`dlaqp2`: at each step the remaining column with the largest norm moves to the front, with norms updated
 * by the usual downdating formula and recomputed when cancellation makes it unreliable. Ties go to the lowest column
 * index, so the result is a deterministic function of the input on every target.
 *
 * ## Rank
 *
 * The numerical rank is the number of leading diagonal entries of `R` with `|r[k,k]| > tol`, where
 * `tol = max(m, n) * eps * |r[0,0]|` and `eps = 2^-52` (the `dgeqp3` convention, ADR-0025). Columns past the rank,
 * in permuted order, are *aliased*: they are numerically combinations of the earlier ones. Rank deficiency is
 * reported through [rank], [aliasedColumns], and [isAliased], never thrown, and [solve] gives aliased columns a
 * coefficient of zero (a basic solution; the minimum-norm solution is not computed). A matrix with more columns than
 * rows has rank at most `m`; a zero matrix has rank 0.
 *
 * With `pivot = false` the columns stay in order and the rank counts only the leading diagonal entries above the
 * tolerance, stopping at the first one that is not. Without pivoting the diagonal of `R` is not ordered, so a
 * numerically dependent column can be followed by columns that look fine; use `pivot = true` whenever rank matters.
 *
 * @see net.maizegenetics.matrixalgebra.decomposition.QRDecomposition
 * @see net.maizegenetics.stats.linearmodels.SolveByOrthogonalizing
 */
class HouseholderQr private constructor(
    /** Number of rows of the factored matrix. */
    val rows: Int,
    /** Number of columns of the factored matrix. */
    val columns: Int,
    /** `R` on and above the diagonal, the reflector vectors below it. `rows x columns`, column-major. */
    private val factor: DenseMatrix,
    private val tau: DoubleArray,
    private val permutation: IntArray,
    /** Numerical rank, from the tolerance described on the class. */
    val rank: Int,
) {
    /** Number of Householder reflectors, `min(rows, columns)`. */
    private val reflectors: Int get() = tau.size

    /** `permutation[k]` is the original index of the column in position `k` of `R`. A copy. */
    fun permutation(): IntArray = permutation.copyOf()

    /** The diagonal of `R`, `min(rows, columns)` values; their magnitudes fall off with the pivoted rank. A copy. */
    val rDiagonal: DoubleArray
        get() = DoubleArray(reflectors) { factor.data[it * rows + it] }

    /** Original indices of the aliased columns, ascending. Empty when the matrix has full column rank. */
    val aliasedColumns: IntArray
        get() = permutation.copyOfRange(rank, columns).also { it.sort() }

    /** True when the column at [column], an original index, is numerically dependent on the others (not estimable). */
    fun isAliased(column: Int): Boolean {
        require(column in 0 until columns) { "column $column is outside 0 until $columns" }
        for (k in rank until columns) if (permutation[k] == column) return true
        return false
    }

    /** `R` as a new `min(rows, columns) x columns` upper-trapezoidal matrix, columns in permuted order. */
    fun r(): DenseMatrix {
        val result = DenseMatrix(reflectors, columns)
        for (c in 0 until columns) {
            for (r in 0..min(c, reflectors - 1)) result.data[c * reflectors + r] = factor.data[c * rows + r]
        }
        return result
    }

    /** Replaces [b], a vector of [rows] values, with `Q' b`. */
    fun applyQt(b: DoubleArray) {
        require(b.size == rows) { "b has ${b.size} values, expected $rows" }
        applyQt(b, 0)
    }

    /** Replaces every column of [b], which has [rows] rows, with `Q'` times it. */
    fun applyQt(b: DenseMatrix) {
        require(b.rows == rows) { "b has ${b.rows} rows, expected $rows" }
        for (c in 0 until b.columns) applyQt(b.data, c * rows)
    }

    /** Replaces [b], a vector of [rows] values, with `Q b`. */
    fun applyQ(b: DoubleArray) {
        require(b.size == rows) { "b has ${b.size} values, expected $rows" }
        applyQ(b, 0)
    }

    /** Replaces every column of [b], which has [rows] rows, with `Q` times it. */
    fun applyQ(b: DenseMatrix) {
        require(b.rows == rows) { "b has ${b.rows} rows, expected $rows" }
        for (c in 0 until b.columns) applyQ(b.data, c * rows)
    }

    /**
     * Least-squares solution of `A x = b`: the [rank] basic coefficients come from `R11 z = (Q' b)[0 until rank]`
     * and land at their original column positions, and aliased columns get `0.0`. For a full-rank `A` with
     * `rows >= columns` this minimizes `|A x - b|`. The residual sum of squares is the squared norm of
     * `(Q' b)[rank until rows]`, available from [applyQt].
     */
    fun solve(b: DoubleArray): DoubleArray {
        require(b.size == rows) { "b has ${b.size} values, expected $rows" }
        val work = b.copyOf()
        applyQt(work, 0)
        Triangular.solveUpper(factor, rank, work)
        val coefficients = DoubleArray(columns)
        for (k in 0 until rank) coefficients[permutation[k]] = work[k]
        return coefficients
    }

    /**
     * [solve] followed by iterative refinement against the original matrix [a], for badly conditioned problems where
     * one pass of Householder QR keeps fewer digits than the data deserves (Longley, the Wampler polynomials).
     *
     * [a] must be the matrix this factorization came from. The refinement is Björck's (1967) method for least
     * squares: it refines the coefficients and the residual together through the augmented system, reusing this
     * factorization for every correction, so it works when the residual is large (Wampler5) and not only when it is
     * small. The residuals of the augmented system, `b - r - A x` and `A' r`, are accumulated in double-double
     * arithmetic ([CompensatedSum]); that is what lets the corrections carry information the working precision lost.
     * Each step costs two applications of `Q` or `Q'` and two passes over `a`.
     *
     * It stops after [maxIterations] steps or once a correction is below `2^-52` of the coefficients. With
     * `maxIterations = 0` it equals [solve]. Aliased columns stay at zero. Like [solve] it is a deterministic
     * function of its inputs on every target.
     */
    fun solveRefined(a: DenseMatrix, b: DoubleArray, maxIterations: Int = DEFAULT_REFINEMENTS): DoubleArray =
        solveRefinedWithResiduals(a, b, maxIterations).coefficients

    /**
     * As [solveRefined], also returning the refined residual vector `b - A x`, whose squared norm is the residual
     * sum of squares. The refinement carries the residual along with the coefficients, so this is more accurate than
     * recomputing `b - A x` in working precision, and it is what a standard error's `sigma` should come from on a
     * badly conditioned design.
     */
    fun solveRefinedWithResiduals(
        a: DenseMatrix,
        b: DoubleArray,
        maxIterations: Int = DEFAULT_REFINEMENTS,
    ): RefinedSolution {
        require(a.rows == rows && a.columns == columns) {
            "a is ${a.rows} x ${a.columns}, but this factorization is of a $rows x $columns matrix"
        }
        require(b.size == rows) { "b has ${b.size} values, expected $rows" }
        require(maxIterations >= 0) { "maxIterations must not be negative, got $maxIterations" }
        val coefficients = DoubleArray(columns)
        if (rank == 0) return RefinedSolution(coefficients, b.copyOf())

        // Start from the plain solution and its residual, Q [0; (Q'b)[rank until rows]].
        val qtb = b.copyOf()
        applyQt(qtb, 0)
        val x = qtb.copyOf(rank)
        Triangular.solveUpper(factor, rank, x)
        val residual = DoubleArray(rows)
        qtb.copyInto(residual, rank, rank, rows)
        applyQ(residual, 0)

        val work = RefinementWork(rows, rank)
        val zero = DoubleArray(rank)
        repeat(maxIterations) {
            if (refineStep(a, b, zero, residual, x, work)) return RefinedSolution(scatter(x, coefficients), residual)
        }
        return RefinedSolution(scatter(x, coefficients), residual)
    }

    private class RefinementWork(rows: Int, rank: Int) {
        val sum = CompensatedSum()
        val f = DoubleArray(rows)
        val h = DoubleArray(rank)
        val dx = DoubleArray(rank)
    }

    /**
     * One Björck step on the augmented system `[I B; B' 0] [u; v] = [f; g]`, where `B` is the basic columns of [a] in
     * pivot order: updates [u] (length [rows]) and [v] (length [rank], pivot order) in place. Least squares is
     * `f = b`, `g = 0` (`u` the residual, `v` the coefficients); the minimum-norm solution of `B' w = e` is
     * `f = 0`, `g = e` (`u = w`). Returns true once the correction to [v] is below `2^-52` of `v`.
     */
    private fun refineStep(
        a: DenseMatrix,
        f: DoubleArray,
        g: DoubleArray,
        u: DoubleArray,
        v: DoubleArray,
        work: RefinementWork,
    ): Boolean {
        val sum = work.sum
        val fr = work.f
        val h = work.h
        // fr = f - u - B v and h = R^-T (g - B' u), both with extended-precision accumulation.
        for (i in 0 until rows) {
            sum.reset()
            sum.add(f[i])
            sum.add(-u[i])
            for (k in 0 until rank) sum.addProduct(-a.data[permutation[k] * rows + i], v[k])
            fr[i] = sum.value
        }
        for (k in 0 until rank) {
            sum.reset()
            sum.add(g[k])
            val column = permutation[k] * rows
            for (i in 0 until rows) sum.addProduct(-a.data[column + i], u[i])
            h[k] = sum.value
        }
        Triangular.solveUpperTransposed(factor, rank, h)
        applyQt(fr, 0)
        // dv = R^-1 (fr[0 until rank] - h); du = Q [h; fr[rank until rows]].
        val dv = work.dx
        for (k in 0 until rank) dv[k] = fr[k] - h[k]
        Triangular.solveUpper(factor, rank, dv)
        h.copyInto(fr, 0, 0, rank)
        applyQ(fr, 0)
        var largestStep = 0.0
        var largestValue = 0.0
        for (k in 0 until rank) {
            v[k] += dv[k]
            largestStep = max(largestStep, abs(dv[k]))
            largestValue = max(largestValue, abs(v[k]))
        }
        for (i in 0 until rows) u[i] += fr[i]
        return largestStep <= EPSILON * largestValue
    }

    private fun scatter(basic: DoubleArray, coefficients: DoubleArray): DoubleArray {
        for (k in 0 until rank) coefficients[permutation[k]] = basic[k]
        return coefficients
    }

    /**
     * For each original column `j`, the Euclidean norm of the row of `R11^-1` that belongs to it, which is
     * `sqrt(((A'A)^-1)[j, j])` for the basic columns. The standard error of coefficient `j` in a least-squares fit is
     * therefore `sigma * rInverseRowNorms()[j]`. Aliased columns get `NaN`.
     *
     * `R^-1` is never stored: it is solved one column at a time into a vector of length [rank], and the squares are
     * accumulated per row, about `rank^3 / 3` multiplications and `O(rank)` memory.
     */
    fun rInverseRowNorms(): DoubleArray {
        val sumsOfSquares = DoubleArray(rank)
        val column = DoubleArray(rank)
        for (j in 0 until rank) {
            column.fill(0.0, 0, j + 1)
            column[j] = 1.0
            // Column j of R^-1 is zero below row j, so only the leading (j + 1) block is involved.
            Triangular.solveUpper(factor, j + 1, column)
            for (i in 0..j) sumsOfSquares[i] += column[i] * column[i]
        }
        val norms = DoubleArray(columns) { Double.NaN }
        for (k in 0 until rank) norms[permutation[k]] = sqrt(sumsOfSquares[k])
        return norms
    }

    /**
     * [rInverseRowNorms] with each norm refined against the original matrix [a], for designs where `R` alone is only
     * good to `cond * eps` (Longley's standard errors keep 11.8 digits otherwise).
     *
     * Column `j`'s norm is `|w|` for the minimum-norm solution `w` of `B' w = e_j`, because
     * `((B'B)^-1)[j, j] = |B^+' e_j|^2`. Starting from the one from `R`, Björck's refinement (the step behind
     * [solveRefined], on the augmented system with `f = 0`, `g = e_j`) improves `w` with double-double residuals, so
     * the answer is limited by the data, not by the factorization. Cost is [maxIterations] steps per estimable
     * column, each `O(rows * rank)`. Aliased columns get `NaN`.
     */
    fun rInverseRowNorms(a: DenseMatrix, maxIterations: Int = DEFAULT_REFINEMENTS): DoubleArray {
        require(a.rows == rows && a.columns == columns) {
            "a is ${a.rows} x ${a.columns}, but this factorization is of a $rows x $columns matrix"
        }
        require(maxIterations >= 0) { "maxIterations must not be negative, got $maxIterations" }
        val norms = DoubleArray(columns) { Double.NaN }
        val work = RefinementWork(rows, rank)
        val zero = DoubleArray(rows)
        val unit = DoubleArray(rank)
        val u = DoubleArray(rows)
        val v = DoubleArray(rank)
        for (j in 0 until rank) {
            unit.fill(0.0)
            unit[j] = 1.0
            // Unrefined start: w = Q [R^-T e_j; 0], v = -R^-1 R^-T e_j.
            u.fill(0.0)
            unit.copyInto(u, 0, 0, rank)
            Triangular.solveUpperTransposed(factor, rank, u)
            u.copyInto(v, 0, 0, rank)
            Triangular.solveUpper(factor, rank, v)
            for (k in 0 until rank) v[k] = -v[k]
            applyQ(u, 0)
            for (step in 0 until maxIterations) {
                if (refineStep(a, zero, unit, u, v, work)) break
            }
            norms[permutation[j]] = norm2(u, 0, rows)
        }
        return norms
    }

    private fun applyQt(b: DoubleArray, offset: Int) {
        for (k in 0 until reflectors) reflect(k, b, offset)
    }

    private fun applyQ(b: DoubleArray, offset: Int) {
        for (k in reflectors - 1 downTo 0) reflect(k, b, offset)
    }

    /** `b = H_k b` for the `k`-th reflector, on the vector starting at [offset] in [b]. */
    private fun reflect(k: Int, b: DoubleArray, offset: Int) {
        val t = tau[k]
        if (t == 0.0) return
        val tail = rows - k - 1
        val kernels = denseKernels
        val w = t * (b[offset + k] + kernels.dot(tail, factor.data, k * rows + k + 1, b, offset + k + 1))
        b[offset + k] -= w
        kernels.axpy(tail, -w, factor.data, k * rows + k + 1, b, offset + k + 1)
    }

    companion object {
        // 2^-52, the gap between 1.0 and the next Double. A literal keeps it identical on every target.
        private const val EPSILON = 2.220446049250313E-16

        /** Refinement steps [solveRefined] takes by default; each gains digits until the data's own precision. */
        const val DEFAULT_REFINEMENTS = 5

        // sqrt(EPSILON): below this a downdated column norm has lost about half its digits and is recomputed.
        private const val RECOMPUTE_THRESHOLD = 1.4901161193847656E-8

        /**
         * Factors [matrix], leaving it untouched. With [pivot] true (the default) columns are pivoted by remaining
         * norm, which is what makes the rank reliable.
         */
        fun of(matrix: DenseMatrix, pivot: Boolean = true): HouseholderQr {
            val m = matrix.rows
            val n = matrix.columns
            val steps = min(m, n)
            val a = matrix.copy()
            val data = a.data
            val tau = DoubleArray(steps)
            val permutation = IntArray(n) { it }
            val kernels = denseKernels

            // Partial column norms: current (below the active row) and original, for the downdating test.
            val current = DoubleArray(n) { norm2(data, it * m, m) }
            val original = current.copyOf()

            for (k in 0 until steps) {
                if (pivot) {
                    var best = k
                    for (j in k + 1 until n) if (current[j] > current[best]) best = j
                    if (best != k) {
                        for (r in 0 until m) {
                            val swap = data[k * m + r]
                            data[k * m + r] = data[best * m + r]
                            data[best * m + r] = swap
                        }
                        val index = permutation[k]
                        permutation[k] = permutation[best]
                        permutation[best] = index
                        current[best] = current[k]
                        original[best] = original[k]
                    }
                }

                // Generate the reflector that zeroes column k below the diagonal (LAPACK dlarfg).
                val tail = m - k - 1
                val alpha = data[k * m + k]
                val tailNorm = norm2(data, k * m + k + 1, tail)
                if (tailNorm == 0.0) {
                    tau[k] = 0.0
                } else {
                    val beta = if (alpha >= 0.0) -pythag(alpha, tailNorm) else pythag(alpha, tailNorm)
                    tau[k] = (beta - alpha) / beta
                    val scale = 1.0 / (alpha - beta)
                    for (r in k + 1 until m) data[k * m + r] *= scale
                    data[k * m + k] = beta
                }

                // Apply it to the trailing columns and downdate their norms.
                for (j in k + 1 until n) {
                    val t = tau[k]
                    if (t != 0.0) {
                        val w = t * (data[j * m + k] + kernels.dot(tail, data, k * m + k + 1, data, j * m + k + 1))
                        data[j * m + k] -= w
                        kernels.axpy(tail, -w, data, k * m + k + 1, data, j * m + k + 1)
                    }
                    if (pivot && current[j] != 0.0) {
                        val ratio = abs(data[j * m + k]) / current[j]
                        val remaining = max(0.0, 1.0 - ratio * ratio)
                        val drift = current[j] / original[j]
                        if (remaining * drift * drift <= RECOMPUTE_THRESHOLD) {
                            current[j] = norm2(data, j * m + k + 1, tail)
                            original[j] = current[j]
                        } else {
                            current[j] *= sqrt(remaining)
                        }
                    }
                }
            }

            var rank = 0
            if (steps > 0) {
                val tolerance = max(m, n) * EPSILON * abs(data[0])
                while (rank < steps && abs(data[rank * m + rank]) > tolerance) rank++
            }
            return HouseholderQr(m, n, a, tau, permutation, rank)
        }
    }
}
