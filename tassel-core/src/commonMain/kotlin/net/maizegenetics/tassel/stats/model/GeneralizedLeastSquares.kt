package net.maizegenetics.tassel.stats.model

import net.maizegenetics.tassel.math.Cholesky
import net.maizegenetics.tassel.math.CholeskyResult
import net.maizegenetics.tassel.math.DenseMatrix
import net.maizegenetics.tassel.math.Triangular
import kotlin.math.sqrt

/**
 * A whitening operator `W` with `W' W = V^-1` for a covariance `V` (up to a variance scale): applied to the response
 * and every design column it turns generalized least squares into ordinary least squares. A covariance model
 * supplies it (a kinship-based one for MLM, ADR-0028); [fromCovariance] and [diagonal] are the simple cases.
 */
interface Whitener {
    /** Number of observations the operator acts on. */
    val size: Int

    /** Replaces [vector], of [size] values, with `W vector`. */
    fun whiten(vector: DoubleArray)

    companion object {
        /** `W = L^-1` for the Cholesky factor `V = L L'`: a forward solve, no inverse. [covariance] must be positive definite. */
        fun fromCovariance(covariance: DenseMatrix): Whitener {
            val factor = when (val result = Cholesky.of(covariance)) {
                is CholeskyResult.Success -> result.factor.lower()
                is CholeskyResult.NotPositiveDefinite ->
                    throw IllegalArgumentException("The covariance is not positive definite (pivot ${result.pivot})")
            }
            return object : Whitener {
                override val size: Int get() = factor.rows

                override fun whiten(vector: DoubleArray) {
                    require(vector.size == size) { "vector has ${vector.size} values, expected $size" }
                    Triangular.solveLower(factor, size, vector)
                }
            }
        }

        /** `V = diag(1 / weights)`, the weighted-least-squares case: each value is multiplied by `sqrt(weight)`. */
        fun diagonal(weights: DoubleArray): Whitener {
            weights.forEachIndexed { index, w ->
                require(w > 0.0 && w.isFinite()) { "weight $w at $index must be positive" }
            }
            val roots = DoubleArray(weights.size) { sqrt(weights[it]) }
            return object : Whitener {
                override val size: Int get() = roots.size

                override fun whiten(vector: DoubleArray) {
                    require(vector.size == size) { "vector has ${vector.size} values, expected $size" }
                    for (i in vector.indices) vector[i] *= roots[i]
                }
            }
        }
    }
}

/**
 * Generalized least squares: the same fit, tests, and sums of squares as OLS, after a [Whitener] transforms the
 * response and design. [whiten] returns the transformed [DesignMatrix] (same terms and columns), which
 * [LinearModelFit], [Anova], [FTest], [WaldTest], and [TTest] accept as they are. MLM with known variance
 * components (tassel-6#197) uses it.
 */
object GeneralizedLeastSquares {
    /** The design and response with [whitener] applied to the response and to every column. */
    fun whiten(design: DesignMatrix, whitener: Whitener): DesignMatrix {
        require(whitener.size == design.rows) {
            "The whitener acts on ${whitener.size} values, the design has ${design.rows} rows"
        }
        val matrix = design.matrix.copy()
        val column = DoubleArray(design.rows)
        for (c in 0 until matrix.columns) {
            matrix.data.copyInto(column, 0, c * design.rows, (c + 1) * design.rows)
            whitener.whiten(column)
            column.copyInto(matrix.data, c * design.rows)
        }
        val response = design.response.copyOf().also { whitener.whiten(it) }
        return design.withData(matrix, response)
    }

    /**
     * Fits [design] by generalized least squares. Coefficients and standard errors are the GLS ones (standard errors
     * with the residual variance of the whitened fit); residuals and fitted values are in the whitened scale.
     */
    fun fit(design: DesignMatrix, whitener: Whitener, refine: Boolean = false): LinearModelFit =
        LinearModelFit.of(whiten(design, whitener), refine)
}
