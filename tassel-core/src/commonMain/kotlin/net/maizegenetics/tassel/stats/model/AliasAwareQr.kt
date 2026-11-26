package net.maizegenetics.tassel.stats.model

import net.maizegenetics.tassel.math.DenseMatrix
import net.maizegenetics.tassel.math.HouseholderQr
import kotlin.math.sqrt

/**
 * An unpivoted Householder QR of the design's usable columns, found in column order the way R's `lm` finds them: a
 * column that is numerically a combination of the columns before it is *aliased* and left out, and later columns are
 * judged without it. (A pivoted QR also detects rank, but chooses which of several dependent columns to drop by
 * norm; a model needs the later term to be the aliased one.)
 *
 * Each column is first scaled to unit length, so the rank tolerance of [HouseholderQr] (`max(n, p) * eps` against a
 * leading diagonal of about one) is a *relative* test on every column and the factorization is better conditioned.
 * The scaling is undone in the coefficients and standard errors. An all-zero column is aliased outright. After the
 * loop, [qr] is full column rank over [kept].
 */
internal class AliasAwareQr private constructor(
    val rows: Int,
    val columns: Int,
    /** Original indices of the columns the factorization uses, ascending. */
    val kept: IntArray,
    /** `1 / norm` of each original column; zero for an all-zero column. */
    val scale: DoubleArray,
    /** The scaled kept columns, `rows x kept.size`. */
    val scaled: DenseMatrix,
    /** Null when no column is kept. */
    val qr: HouseholderQr?,
    /** Original indices of aliased columns, ascending. */
    val aliased: IntArray,
) {
    val rank: Int get() = kept.size

    companion object {
        /** [relativeTolerance] replaces the default rank tolerance of [HouseholderQr] (see its `of`). */
        fun of(matrix: DenseMatrix, relativeTolerance: Double? = null): AliasAwareQr {
            val m = matrix.rows
            val scale = DoubleArray(matrix.columns) {
                val norm = columnNorm(matrix, it)
                if (norm > 0.0 && norm.isFinite()) 1.0 / norm else 0.0
            }
            val kept = (0 until matrix.columns).filter { scale[it] != 0.0 }.toMutableList()
            val aliased = (0 until matrix.columns).filter { scale[it] == 0.0 }.toMutableList()
            while (true) {
                val scaled = scaledColumns(matrix, kept, scale)
                if (kept.isEmpty()) {
                    return AliasAwareQr(
                        m,
                        matrix.columns,
                        IntArray(0),
                        scale,
                        scaled,
                        null,
                        aliased.toIntArray(),
                    )
                }
                val qr = HouseholderQr.of(scaled, pivot = false, relativeTolerance = relativeTolerance)
                if (qr.rank == kept.size) {
                    return AliasAwareQr(
                        m,
                        matrix.columns,
                        kept.toIntArray(),
                        scale,
                        scaled,
                        qr,
                        aliased.sorted().toIntArray(),
                    )
                }
                // The first column whose diagonal fell below the tolerance depends on those before it.
                aliased += kept.removeAt(qr.rank)
            }
        }

        private fun columnNorm(matrix: DenseMatrix, column: Int): Double {
            var largest = 0.0
            for (i in 0 until matrix.rows) largest = maxOf(largest, kotlin.math.abs(matrix[i, column]))
            if (largest == 0.0 || !largest.isFinite()) return largest
            var sum = 0.0
            for (i in 0 until matrix.rows) (matrix[i, column] / largest).let { sum += it * it }
            return largest * sqrt(sum)
        }

        private fun scaledColumns(matrix: DenseMatrix, kept: List<Int>, scale: DoubleArray): DenseMatrix {
            val result = DenseMatrix(matrix.rows, kept.size)
            kept.forEachIndexed { target, source ->
                for (i in 0 until matrix.rows) result.data[target * matrix.rows + i] = matrix[i, source] * scale[source]
            }
            return result
        }
    }
}
