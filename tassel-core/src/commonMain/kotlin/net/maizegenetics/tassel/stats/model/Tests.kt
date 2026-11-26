package net.maizegenetics.tassel.stats.model

import net.maizegenetics.tassel.math.DenseMatrix
import net.maizegenetics.tassel.math.HouseholderQr
import net.maizegenetics.tassel.math.Triangular
import net.maizegenetics.tassel.stats.ChiSquaredDistribution
import net.maizegenetics.tassel.stats.FDistribution
import net.maizegenetics.tassel.stats.StudentTDistribution
import kotlin.math.abs
import kotlin.math.min
import kotlin.math.sqrt

/** Which distribution a [TestResult]'s statistic follows under the null. */
enum class TestKind { F, T, ChiSquare }

/**
 * The outcome of a test: the [statistic], its degrees of freedom, and the upper-tail [p] with `log10(p)`, which
 * stays finite when `p` underflows to zero (ADR-0025).
 *
 * [df] is the numerator (or only) degrees of freedom; [residualDf] is the denominator of an F or t test and `NaN`
 * for a chi-square test. A test that cannot be carried out (aliased coefficients, no numerator or residual degrees
 * of freedom) has [estimable] false and `NaN` for the statistic and p-values.
 */
class TestResult(
    val kind: TestKind,
    val statistic: Double,
    val df: Double,
    val residualDf: Double,
    val p: Double,
    val log10P: Double,
    val estimable: Boolean = true,
) {
    override fun toString(): String =
        "TestResult($kind, statistic=$statistic, df=$df, residualDf=$residualDf, p=$p, log10P=$log10P)"

    companion object {
        internal fun notEstimable(kind: TestKind, df: Double = Double.NaN, residualDf: Double = Double.NaN) =
            TestResult(kind, Double.NaN, df, residualDf, Double.NaN, Double.NaN, estimable = false)

        internal fun fromF(f: Double, df1: Int, df2: Int): TestResult {
            if (df1 <= 0 || df2 <= 0 || f.isNaN()) return notEstimable(TestKind.F, df1.toDouble(), df2.toDouble())
            val distribution = FDistribution(df1.toDouble(), df2.toDouble())
            return TestResult(
                TestKind.F,
                f,
                df1.toDouble(),
                df2.toDouble(),
                distribution.sf(f),
                distribution.log10Sf(f),
            )
        }
    }
}

/**
 * F test of nested models: does [full] explain more than [reduced] by more than chance? Numerator degrees of freedom
 * are the difference of ranks (any size), the denominator is [full]'s residual mean square.
 *
 * @see net.maizegenetics.stats.linearmodels.LinearModelUtils.Ftest
 */
object FTest {
    /**
     * Both fits must be of the same observations in the same order, and [reduced] must be [full] without some
     * columns (this is not checked beyond shape and rank). `F = ((RSS0 - RSS1) / (r1 - r0)) / (RSS1 / (n - r1))`.
     */
    fun compare(reduced: LinearModelFit, full: LinearModelFit): TestResult {
        require(reduced.design.observations.contentEquals(full.design.observations)) {
            "The two models were fitted to different observations"
        }
        val numerator = full.rank - reduced.rank
        val denominator = full.residualDf
        if (numerator <= 0 ||
            denominator <= 0
        ) {
            return TestResult.notEstimable(TestKind.F, numerator.toDouble(), denominator.toDouble())
        }
        val extra = (reduced.residualSumOfSquares - full.residualSumOfSquares).coerceAtLeast(0.0)
        return TestResult.fromF((extra / numerator) / (full.residualSumOfSquares / denominator), numerator, denominator)
    }
}

/** t test of one coefficient or one linear combination, two-sided. */
object TTest {
    /** Tests that the coefficient of design column [column] is zero. */
    fun coefficient(fit: LinearModelFit, column: Int): TestResult {
        require(column in 0 until fit.design.columnCount) {
            "column $column is outside 0 until ${fit.design.columnCount}"
        }
        if (!fit.isEstimable(column)) return TestResult.notEstimable(TestKind.T, 1.0, fit.residualDf.toDouble())
        return twoSided(fit.coefficients[column] / fit.standardErrors[column], fit.residualDf)
    }

    /** Tests `l' beta = value`, with [l] one weight per design column. Aliased columns must have weight zero. */
    fun contrast(fit: LinearModelFit, l: DoubleArray, value: Double = 0.0): TestResult {
        val estimate = ContrastGeometry.of(fit, DenseMatrix.ofRows(l), doubleArrayOf(value))
            ?: return TestResult.notEstimable(TestKind.T, 1.0, fit.residualDf.toDouble())
        if (estimate.rank != 1) return TestResult.notEstimable(TestKind.T, 1.0, fit.residualDf.toDouble())
        val t = (estimate.difference[0]) / (sqrt(fit.residualVariance) * estimate.scale)
        return twoSided(t, fit.residualDf)
    }

    private fun twoSided(t: Double, df: Int): TestResult {
        if (df <= 0 || t.isNaN()) return TestResult.notEstimable(TestKind.T, 1.0, df.toDouble())
        val distribution = StudentTDistribution(df.toDouble())
        val a = abs(t)
        val p = min(1.0, 2.0 * distribution.sf(a))
        val log10P = min(0.0, (LN_2 + distribution.logSf(a)) / LN_10)
        return TestResult(TestKind.T, t, 1.0, df.toDouble(), p, log10P)
    }
}

/**
 * Wald test of the linear hypothesis `L beta = c`, with `L` of `q` rows, one weight per design column. The
 * statistic is `(L b - c)' [L (X'X)^-1 L']^-1 (L b - c)`, computed with triangular solves and a QR of the `q x q`
 * system, never an inverse. Dependent rows of `L` reduce the numerator degrees of freedom to the rank of the
 * hypothesis; a row that involves an aliased column is not estimable.
 */
object WaldTest {
    /**
     * The F form `W / (q sigma^2)` against `F(q, residualDf)`, with `sigma^2` estimated from the fit; or, when
     * [knownVariance] is given, the chi-square form `W / sigma^2` against `chi-square(q)`.
     */
    fun test(
        fit: LinearModelFit,
        l: DenseMatrix,
        value: DoubleArray = DoubleArray(l.rows),
        knownVariance: Double? = null,
    ): TestResult {
        val kind = if (knownVariance == null) TestKind.F else TestKind.ChiSquare
        val geometry = ContrastGeometry.of(fit, l, value) ?: return TestResult.notEstimable(kind)
        val q = geometry.rank
        if (q == 0) return TestResult.notEstimable(kind, 0.0, fit.residualDf.toDouble())
        if (knownVariance == null) {
            return TestResult.fromF(geometry.quadraticForm / (q * fit.residualVariance), q, fit.residualDf)
        }
        require(knownVariance > 0.0) { "knownVariance must be positive, got $knownVariance" }
        val statistic = geometry.quadraticForm / knownVariance
        val distribution = ChiSquaredDistribution(q.toDouble())
        return TestResult(
            kind,
            statistic,
            q.toDouble(),
            Double.NaN,
            distribution.sf(statistic),
            distribution.log10Sf(statistic),
        )
    }
}

/**
 * What the Wald and contrast tests share. With `X = Q R` (columns scaled by `s`, kept columns only), `Var(L b) =
 * sigma^2 M'M` for `M = R^-T (L S)'`, so the quadratic form of `d = L b - c` is `|R2^-T P' d|^2` for a pivoted QR
 * `M P = Q2 R2`: triangular solves only. [scale] is `|M|` for a single row.
 */
private class ContrastGeometry(
    val rank: Int,
    val quadraticForm: Double,
    val difference: DoubleArray,
    val scale: Double,
) {
    companion object {
        fun of(fit: LinearModelFit, l: DenseMatrix, value: DoubleArray): ContrastGeometry? {
            val design = fit.design
            require(l.columns == design.columnCount) {
                "L has ${l.columns} columns, the design has ${design.columnCount}"
            }
            require(value.size == l.rows) { "the right-hand side has ${value.size} values, L has ${l.rows} rows" }
            val q = l.rows
            val decomposition = design.decomposition
            for (row in 0 until q) {
                for (column in decomposition.aliased) if (l[row, column] != 0.0) return null
            }
            val k = fit.rank
            val difference = DoubleArray(q) { row ->
                var sum = -value[row]
                for (column in decomposition.kept) sum += l[row, column] * fit.coefficients[column]
                sum
            }
            if (k == 0) return ContrastGeometry(0, 0.0, difference, Double.NaN)
            val r = decomposition.qr!!.r()
            // M is k x q: column `row` is R^-T (scaled L row).
            val m = DenseMatrix(k, q)
            for (row in 0 until q) {
                val work = DoubleArray(k) { index ->
                    decomposition.kept[index].let {
                        l[row, it] *
                            decomposition.scale[it]
                    }
                }
                Triangular.solveUpperTransposed(r, k, work)
                work.copyInto(m.data, row * k)
            }
            val qrM = HouseholderQr.of(m, pivot = true)
            val rank = qrM.rank
            val permutation = qrM.permutation()
            val z = DoubleArray(rank) { difference[permutation[it]] }
            Triangular.solveUpperTransposed(qrM.r(), rank, z)
            var form = 0.0
            for (v in z) form += v * v
            val scale = if (q == 1) {
                var sum = 0.0
                for (i in 0 until k) sum += m.data[i] * m.data[i]
                sqrt(sum)
            } else {
                Double.NaN
            }
            return ContrastGeometry(rank, form, difference, scale)
        }
    }
}

private const val LN_2 = 0.6931471805599453
private const val LN_10 = 2.302585092994046
