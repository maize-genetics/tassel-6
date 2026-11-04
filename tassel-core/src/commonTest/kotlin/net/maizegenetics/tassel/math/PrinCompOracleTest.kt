package net.maizegenetics.tassel.math

import kotlin.math.abs
import kotlin.math.max
import kotlin.math.sqrt
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Checks [SymmetricEigenDecomposition] against the R `prcomp` values in TASSEL 5's `PrinCompTest`, through the Gram
 * matrix route that the `pca` operation (#19) will take.
 *
 * TASSEL 5's `PrinComp` takes the SVD of `A`, the centered (and for `corr`, scaled) data matrix divided by
 * `sqrt(nRows - 1)`. With `A = U S V'` this test instead eigendecomposes the small `nRows x nRows` Gram matrix
 * `G = A A'`, whose eigenvalues are `S^2`, and recovers what TASSEL 5 prints:
 *
 * - eigenvalues: the eigenvalues of `G`;
 * - eigenvectors in variable space: `V = A' U S^-1`;
 * - principal components: `A V * sqrt(nRows - 1)`, which equals `sqrt(nRows - 1) * U * S`.
 *
 * R's values have about seven significant digits, and forming `G` squares the condition number, so each value is
 * compared with a tolerance of `max(1e-9, |expected| * 1e-6)`; the relative part is the same one `PrinCompTest`
 * uses. Columns whose eigenvalue is numerically zero are skipped, as in `PrinCompTest`, because their eigenvectors
 * are arbitrary. Eigenvectors are compared up to sign (ADR-0009), with one sign per column shared by the
 * eigenvector and the principal component, so a sign flip cannot hide a disagreement between them.
 *
 * @see net.maizegenetics.stats.PCA.PrinComp
 */
class PrinCompOracleTest {
    @Test
    fun covarianceDataSetMatchesRPrcomp() {
        val data =
            arrayOf(
                doubleArrayOf(1.0, 2.0, 3.0, 2.0, 3.0, 4.0),
                doubleArrayOf(3.0, 2.0, 3.0, 4.0, 5.0, 5.0),
                doubleArrayOf(4.0, 3.0, 4.0, 5.0, 6.0, 5.0),
                doubleArrayOf(4.0, 8.0, 7.0, 3.0, 4.0, 5.0),
            )
        val expectedEigenvalues =
            doubleArrayOf(1.265282e+01, 4.716365e+00, 4.747914e-02, 2.501067e-31)
        val expectedEigenvectors =
            arrayOf(
                doubleArrayOf(-0.26040315, 0.49115810, 0.2994106, -0.34767622),
                doubleArrayOf(-0.80318313, -0.13582600, 0.1136607, -0.36391773),
                doubleArrayOf(-0.53113010, -0.03333491, -0.4290802, 0.64855531),
                doubleArrayOf(0.01164988, 0.59364919, -0.2433303, 0.03151932),
                doubleArrayOf(0.01164988, 0.59364919, -0.2433303, 0.03151932),
                doubleArrayOf(-0.06867629, 0.18537818, 0.7713011, 0.56927516),
            )
        val expectedComponents =
            arrayOf(
                doubleArrayOf(2.6068470, -2.6229333, -0.10986215, 3.885781e-16),
                doubleArrayOf(2.0639639, 0.9193579, 0.28693917, -3.053113e-16),
                doubleArrayOf(0.4925473, 2.4286535, -0.21573033, -2.775558e-17),
                doubleArrayOf(-5.1633582, -0.7250780, 0.03865331, 2.775558e-17),
            )

        assertMatchesPrcomp(data, scaleColumns = false, expectedEigenvalues, expectedEigenvectors, expectedComponents)
    }

    @Test
    fun correlationDataSetMatchesRPrcomp() {
        val data =
            arrayOf(
                doubleArrayOf(9.0, 5.0, 7.0, 9.0, 4.0, 7.0, 1.0),
                doubleArrayOf(2.0, 1.0, 1.0, 1.0, 1.0, 7.0, 6.0),
                doubleArrayOf(7.0, 6.0, 9.0, 8.0, 4.0, 2.0, 7.0),
                doubleArrayOf(2.0, 6.0, 8.0, 7.0, 9.0, 8.0, 4.0),
                doubleArrayOf(3.0, 7.0, 4.0, 7.0, 4.0, 4.0, 4.0),
            )
        val expectedEigenvalues =
            doubleArrayOf(3.537200e+00, 1.688439e+00, 1.278101e+00, 4.962601e-01, 3.970544e-32)
        val expectedEigenvectors =
            arrayOf(
                doubleArrayOf(0.3123845, -0.19660085, -0.64877339, 0.3224638, -0.11369759),
                doubleArrayOf(0.4622523, -0.07987140, 0.26331346, -0.5401578, 0.45748534),
                doubleArrayOf(0.4798065, -0.06589826, 0.12565918, 0.5645608, 0.51010183),
                doubleArrayOf(0.5247024, -0.01732956, -0.12095154, -0.1184232, -0.43041095),
                doubleArrayOf(0.3413618, 0.37243242, 0.51229863, 0.1914193, -0.50807637),
                doubleArrayOf(-0.1426675, 0.72445284, -0.06051857, 0.2737018, 0.26207596),
                doubleArrayOf(-0.2152440, -0.53552586, 0.46176532, 0.3999245, -0.07414278),
            )
        val expectedComponents =
            arrayOf(
                doubleArrayOf(1.2309943, 0.83516558, -1.7306678, 0.08630054, -1.249001e-16),
                doubleArrayOf(-3.2826861, 0.01445715, -0.2166971, 0.24060123, -1.228184e-15),
                doubleArrayOf(1.0824049, -1.94998024, 0.2254010, 0.53510557, -2.053913e-15),
                doubleArrayOf(0.8133263, 1.45798609, 1.3790800, 0.36330020, -3.042705e-15),
                doubleArrayOf(0.1559606, -0.35762858, 0.3428839, -1.22530754, 7.088080e-15),
            )

        assertMatchesPrcomp(data, scaleColumns = true, expectedEigenvalues, expectedEigenvectors, expectedComponents)
    }

    private fun assertMatchesPrcomp(
        data: Array<DoubleArray>,
        scaleColumns: Boolean,
        expectedEigenvalues: DoubleArray,
        expectedEigenvectors: Array<DoubleArray>,
        expectedComponents: Array<DoubleArray>,
    ) {
        val rows = data.size
        val columns = data[0].size
        val a = scaledCenteredMatrix(data, scaleColumns)
        val eigen = SymmetricEigenDecomposition.of(rows) { i, j -> (0 until columns).sumOf { a[i][it] * a[j][it] } }

        expectedEigenvalues.forEachIndexed { index, expected ->
            assertEquals(expected, eigen.eigenvalue(index), tolerance(expected), "eigenvalue $index")
        }

        for (component in expectedEigenvalues.indices) {
            if (abs(eigen.eigenvalue(component)) <= ZERO_EIGENVALUE) continue
            val singularValue = sqrt(eigen.eigenvalue(component))

            // V = A' U S^-1, the eigenvectors in variable space that TASSEL 5 prints.
            val eigenvector =
                DoubleArray(columns) { variable ->
                    (0 until rows).sumOf { a[it][variable] * eigen[it, component] } / singularValue
                }
            // The principal components, sqrt(rows - 1) * U * S.
            val principalComponent =
                DoubleArray(rows) { row -> sqrt(rows - 1.0) * singularValue * eigen[row, component] }

            // One sign per component, taken from the largest entry of the expected eigenvector (ADR-0009).
            var pivot = 0
            for (variable in 0 until columns) {
                if (abs(expectedEigenvectors[variable][component]) > abs(expectedEigenvectors[pivot][component])) {
                    pivot = variable
                }
            }
            val sign = if (expectedEigenvectors[pivot][component] * eigenvector[pivot] < 0.0) -1.0 else 1.0

            for (variable in 0 until columns) {
                val expected = expectedEigenvectors[variable][component]
                assertEquals(
                    expected,
                    sign * eigenvector[variable],
                    tolerance(expected),
                    "eigenvector ($variable, $component)",
                )
            }
            for (row in 0 until rows) {
                val expected = expectedComponents[row][component]
                assertEquals(
                    expected,
                    sign * principalComponent[row],
                    tolerance(expected),
                    "principal component ($row, $component)",
                )
            }
        }

        // The eigenvectors of the Gram matrix are an orthonormal basis whatever the data looked like.
        val leadingNorm = eigen.eigenvector(0).sumOf { it * it }
        assertTrue(abs(leadingNorm - 1.0) < 1e-12, "the leading eigenvector is not unit length: $leadingNorm")
    }

    /** `PrinComp.centerCols` plus, when [scaleColumns] is set, `scaleCenteredMatrix`, then `1 / sqrt(rows - 1)`. */
    private fun scaledCenteredMatrix(data: Array<DoubleArray>, scaleColumns: Boolean): Array<DoubleArray> {
        val rows = data.size
        val columns = data[0].size
        val result = Array(rows) { data[it].copyOf() }
        for (column in 0 until columns) {
            val mean = (0 until rows).sumOf { result[it][column] } / rows
            for (row in 0 until rows) result[row][column] -= mean
            if (scaleColumns) {
                val sumOfSquares = (0 until rows).sumOf { result[it][column] * result[it][column] }
                val standardDeviation = sqrt(sumOfSquares / (rows - 1))
                for (row in 0 until rows) result[row][column] /= standardDeviation
            }
        }
        val multiplier = 1.0 / sqrt(rows - 1.0)
        for (row in 0 until rows) {
            for (column in 0 until columns) result[row][column] *= multiplier
        }
        return result
    }

    private fun tolerance(expected: Double): Double = max(1e-9, abs(expected) * 1e-6)

    private companion object {
        /** `PrinCompTest`'s cutoff for "this eigenvalue is zero, so its eigenvector is arbitrary". */
        const val ZERO_EIGENVALUE = 1e-12
    }
}
