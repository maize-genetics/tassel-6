package net.maizegenetics.tassel.distance

import net.maizegenetics.tassel.math.SymmetricEigenDecomposition
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import kotlin.math.abs
import kotlin.math.max

/** One reason a matrix cannot serve as a covariance (ADR-0026). */
sealed interface KinshipProblem {
    /** A readable description naming the offending cells or values. */
    val message: String

    /** [value] at [row] and [column] is NaN or infinite. [count] cells are, counting this one. */
    data class NonFinite(val row: Int, val column: Int, val value: Double, val count: Int) : KinshipProblem {
        override val message: String
            get() = "$count values are not finite; the first is $value at row $row, column $column"
    }

    /** Cells ([row], [column]) and ([column], [row]) differ by [difference], more than the tolerance. */
    data class Asymmetric(val row: Int, val column: Int, val difference: Double) : KinshipProblem {
        override val message: String
            get() = "the matrix is not symmetric; rows $row and $column differ by $difference at those columns"
    }

    /** The smallest eigenvalue [minEigenvalue] is below [threshold], so the matrix is not positive semidefinite. */
    data class NotPositiveSemidefinite(val minEigenvalue: Double, val threshold: Double) : KinshipProblem {
        override val message: String
            get() = "the matrix is not positive semidefinite; " +
                "its smallest eigenvalue $minEigenvalue is below $threshold"
    }
}

/**
 * The outcome of [KinshipValidation]. The eigenvalues are `null` when the matrix was not symmetric
 * and finite, since they would mean nothing then.
 *
 * @property maxAsymmetry The largest absolute difference between mirrored cells.
 * @property minEigenvalue The smallest eigenvalue.
 * @property maxEigenvalue The largest eigenvalue.
 * @property eigenvalueThreshold `-tolerance * maxEigenvalue`; the smallest eigenvalue must not fall below it.
 */
data class KinshipReport(
    val size: Int,
    val tolerance: Double,
    val maxAsymmetry: Double,
    val minEigenvalue: Double?,
    val maxEigenvalue: Double?,
    val eigenvalueThreshold: Double?,
    val problems: List<KinshipProblem>,
) {
    /** True when nothing stops the matrix being used as a covariance. */
    val isValid: Boolean get() = problems.isEmpty()

    /** Throws if [isValid] is false, with every problem in the message. */
    fun requireValid() {
        require(isValid) { "Kinship matrix is not usable: " + problems.joinToString("; ") { it.message } }
    }
}

/**
 * Checks that a matrix can be a covariance: finite, symmetric within `tolerance`, and positive
 * semidefinite, meaning its smallest eigenvalue is at least `-tolerance * maxEigenvalue`. A
 * singular but positive semidefinite matrix passes. Nothing is repaired; bending or adding a ridge
 * is out of scope (ADR-0026).
 *
 * Unique labels are guaranteed by [net.maizegenetics.tassel.taxa.TaxaList], and squareness by
 * [DistanceMatrix].
 */
object KinshipValidation {
    /** Default relative tolerance on asymmetry (absolute) and on the smallest eigenvalue (relative to the largest). */
    const val DEFAULT_TOLERANCE: Double = 1e-6

    /**
     * Validates [matrix]. A [DistanceMatrix] stores one triangle, so it is symmetric by
     * construction; the asymmetry check matters for [check] on raw values.
     */
    suspend fun check(
        matrix: DistanceMatrix,
        tolerance: Double = DEFAULT_TOLERANCE,
        progress: ProgressReporter = ProgressReporter.None,
    ): KinshipReport =
        check(matrix.numberOfTaxa, tolerance, progress) { row, column -> matrix.distance(row, column).toDouble() }

    /**
     * Validates the `size` by `size` matrix whose entry at `(row, column)` is `element(row, column)`,
     * which is called for the whole square.
     */
    suspend fun check(
        size: Int,
        tolerance: Double = DEFAULT_TOLERANCE,
        progress: ProgressReporter = ProgressReporter.None,
        element: (row: Int, column: Int) -> Double,
    ): KinshipReport {
        require(tolerance >= 0.0 && tolerance.isFinite()) {
            "tolerance must be finite and not negative, not $tolerance"
        }
        val problems = ArrayList<KinshipProblem>()

        var nonFiniteCount = 0
        var firstNonFinite: KinshipProblem.NonFinite? = null
        var maxAsymmetry = 0.0
        var worstAsymmetry: KinshipProblem.Asymmetric? = null
        for (row in 0 until size) {
            for (column in 0 until size) {
                val value = element(row, column)
                if (!value.isFinite()) {
                    nonFiniteCount++
                    if (firstNonFinite == null) firstNonFinite = KinshipProblem.NonFinite(row, column, value, 0)
                } else if (column < row) {
                    val mirror = element(column, row)
                    if (mirror.isFinite()) {
                        val difference = abs(value - mirror)
                        if (difference > maxAsymmetry) {
                            maxAsymmetry = difference
                            worstAsymmetry = KinshipProblem.Asymmetric(row, column, difference)
                        }
                    }
                }
            }
        }
        firstNonFinite?.let { problems += it.copy(count = nonFiniteCount) }
        worstAsymmetry?.takeIf { maxAsymmetry > tolerance }?.let { problems += it }
        progress.report(0, 1, "Checking kinship matrix")

        var minEigenvalue: Double? = null
        var maxEigenvalue: Double? = null
        var threshold: Double? = null
        if (problems.isEmpty() && size > 0) {
            val eigen = SymmetricEigenDecomposition.of(size, { fraction ->
                progress.report((fraction * PROGRESS_UNITS).toLong(), PROGRESS_UNITS, "Checking kinship matrix")
            }, element)
            val values = eigen.eigenvalues()
            val largest = values.max()
            val smallest = values.min()
            maxEigenvalue = largest
            minEigenvalue = smallest
            // A matrix whose largest eigenvalue is not positive has no scale; only zero is semidefinite then.
            threshold = -tolerance * max(largest, 0.0)
            if (smallest < threshold) problems += KinshipProblem.NotPositiveSemidefinite(smallest, threshold)
        }
        progress.report(1, 1, "Checking kinship matrix")
        return KinshipReport(size, tolerance, maxAsymmetry, minEigenvalue, maxEigenvalue, threshold, problems)
    }

    private const val PROGRESS_UNITS = 1000L
}
