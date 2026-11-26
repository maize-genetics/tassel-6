package net.maizegenetics.tassel.stats.model

import kotlin.math.sqrt

/**
 * An ordinary least-squares fit of a [DesignMatrix], from the unpivoted, alias-aware Householder QR of the design
 * ([DesignMatrix.rank], ADR-0025). Nothing is inverted: coefficients solve against `R`, and standard errors are
 * `sigma` times the row norms of `R^-1`.
 *
 * A collinear design does not throw. The effective model degrees of freedom are the [rank], never the column count;
 * [aliasedColumns] lists the columns that are combinations of earlier ones, and they have `NaN` in [coefficients] and
 * [standardErrors] (R's `NA`). Residuals and fitted values are the fit of the estimable columns, so they are the
 * same whichever of several dependent columns is dropped.
 *
 * For a whitened design from [GeneralizedLeastSquares], residuals and fitted values are in the whitened scale, and
 * every test on the fit is the generalized-least-squares test.
 *
 * @see net.maizegenetics.stats.linearmodels.SweepFastLinearModel
 */
class LinearModelFit private constructor(
    /** The design that was fitted. */
    val design: DesignMatrix,
    /** Numerical rank: the model's degrees of freedom. */
    val rank: Int,
    /** One per design column; `NaN` where the column is aliased. */
    val coefficients: DoubleArray,
    /** One per design column; `NaN` where aliased, and everywhere when there are no residual degrees of freedom. */
    val standardErrors: DoubleArray,
    /** Response minus fitted, one per row. */
    val residuals: DoubleArray,
    /** Fitted values, one per row. */
    val fitted: DoubleArray,
    /** Residual sum of squares. */
    val residualSumOfSquares: Double,
) {
    /** Number of observations. */
    val observations: Int get() = design.rows

    /** Residual degrees of freedom, `n - rank`. */
    val residualDf: Int get() = design.rows - rank

    /** Residual mean square, `RSS / residualDf`; `NaN` with no residual degrees of freedom. */
    val residualVariance: Double
        get() = if (residualDf > 0) residualSumOfSquares / residualDf else Double.NaN

    /** Original indices of the aliased columns, ascending. */
    val aliasedColumns: IntArray get() = design.aliasedColumns

    /** Whether the coefficient of [column] is estimable. */
    fun isEstimable(column: Int): Boolean = !coefficients[column].isNaN()

    /** Labels of the terms whose columns are all aliased. */
    val aliasedTerms: List<String> get() = design.aliasedTerms

    companion object {
        /**
         * Fits [design]. With [refine] true the solution is polished with [net.maizegenetics.tassel.math.HouseholderQr.solveRefined]
         * (double-double residuals), for badly conditioned designs; it costs a few passes over the matrix.
         */
        fun of(design: DesignMatrix, refine: Boolean = false): LinearModelFit {
            val n = design.rows
            val p = design.columnCount
            val decomposition = design.decomposition
            val y = design.response
            val k = decomposition.rank
            val coefficients = DoubleArray(p) { Double.NaN }
            val standardErrors = DoubleArray(p) { Double.NaN }
            val qr = decomposition.qr
            if (qr == null) {
                val rss = y.sumOfSquares()
                return LinearModelFit(design, 0, coefficients, standardErrors, y.copyOf(), DoubleArray(n), rss)
            }

            val scaledCoefficients: DoubleArray
            val residuals: DoubleArray
            if (refine) {
                val refined = qr.solveRefinedWithResiduals(decomposition.scaled, y)
                scaledCoefficients = refined.coefficients
                residuals = refined.residuals
            } else {
                scaledCoefficients = qr.solve(y)
                val tail = y.copyOf().also { qr.applyQt(it) }
                tail.fill(0.0, 0, k)
                qr.applyQ(tail)
                residuals = tail
            }
            val rss = residuals.sumOfSquares()
            val fitted = DoubleArray(n) { y[it] - residuals[it] }
            val rowNorms = if (refine) qr.rInverseRowNorms(decomposition.scaled) else qr.rInverseRowNorms()
            val sigma = if (n - k > 0) sqrt(rss / (n - k)) else Double.NaN
            decomposition.kept.forEachIndexed { index, column ->
                coefficients[column] = scaledCoefficients[index] * decomposition.scale[column]
                standardErrors[column] = sigma * rowNorms[index] * decomposition.scale[column]
            }
            return LinearModelFit(design, k, coefficients, standardErrors, residuals, fitted, rss)
        }

        private fun DoubleArray.sumOfSquares(): Double {
            var sum = 0.0
            for (v in this) sum += v * v
            return sum
        }
    }
}
