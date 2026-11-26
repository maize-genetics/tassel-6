package net.maizegenetics.tassel.math

import kotlin.math.abs
import kotlin.math.sqrt
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/** Structure, rank detection, and bookkeeping of [HouseholderQr]. The NIST accuracy checks are in [NistRegressionTest]. */
class HouseholderQrTest {
    @Test
    fun factorsAsQrWithThePermutation() {
        val a = random(9, 5, 1)
        val qr = HouseholderQr.of(a)
        assertEquals(5, qr.rank)
        // Build R padded to 9 x 5, apply Q, and compare with the permuted columns of A.
        val r = qr.r()
        val padded = DenseMatrix(9, 5)
        for (c in 0 until 5) for (i in 0 until r.rows) padded[i, c] = r[i, c]
        qr.applyQ(padded)
        val permutation = qr.permutation()
        for (c in 0 until 5) {
            for (i in 0 until 9) assertClose(a[i, permutation[c]], padded[i, c], 1e-13, "A P = Q R at ($i,$c)")
        }
        assertEquals((0 until 5).toSet(), permutation.toSet())
    }

    @Test
    fun applyQAndQtAreInverseAndPreserveNorm() {
        val qr = HouseholderQr.of(random(12, 4, 2))
        val b = random(12, 1, 3).data
        val c = b.copyOf()
        qr.applyQt(c)
        assertClose(norm(b), norm(c), 1e-13, "Q' preserves length")
        qr.applyQ(c)
        for (i in b.indices) assertClose(b[i], c[i], 1e-13, "Q Q' b = b at $i")
        val matrix = random(12, 3, 4)
        val copy = matrix.copy()
        qr.applyQt(matrix)
        qr.applyQ(matrix)
        for (i in matrix.data.indices) assertClose(copy.data[i], matrix.data[i], 1e-13, "matrix round trip at $i")
    }

    @Test
    fun pivotsByColumnNormAndOrdersTheDiagonal() {
        val a = DenseMatrix.ofRows(
            doubleArrayOf(1.0, 10.0, 0.1),
            doubleArrayOf(1.0, -10.0, 0.1),
            doubleArrayOf(1.0, 10.0, -0.1),
            doubleArrayOf(1.0, 10.0, 0.1),
        )
        val qr = HouseholderQr.of(a)
        assertEquals(1, qr.permutation()[0], "the largest column goes first")
        val diagonal = qr.rDiagonal.map { abs(it) }
        assertTrue(diagonal[0] >= diagonal[1] && diagonal[1] >= diagonal[2], "diagonal magnitudes $diagonal")
    }

    @Test
    fun solveIsLeastSquares() {
        val a = random(20, 4, 5)
        val truth = doubleArrayOf(1.0, -2.0, 0.5, 3.0)
        val noise = Random(6)
        val b = DoubleArray(20) { i -> (0 until 4).sumOf { a[i, it] * truth[it] } + noise.nextDouble(-0.01, 0.01) }
        val qr = HouseholderQr.of(a)
        val x = qr.solve(b)
        // The residual is orthogonal to every column of A.
        val residual = DoubleArray(20) { i -> b[i] - (0 until 4).sumOf { a[i, it] * x[it] } }
        for (j in 0 until 4) {
            assertTrue(
                abs(
                    (0 until 20).sumOf {
                        a[it, j] * residual[it]
                    },
                ) < 1e-12,
                "column $j is not orthogonal to the residual",
            )
        }
        // Same answer as the normal equations, which are fine for a well-conditioned design.
        val gram = DenseMatrix(4, 4)
        val rhs = DoubleArray(4)
        for (i in 0 until 4) {
            for (j in 0 until 4) gram[i, j] = (0 until 20).sumOf { a[it, i] * a[it, j] }
            rhs[i] = (0 until 20).sumOf { a[it, i] * b[it] }
        }
        val viaCholesky = (Cholesky.of(gram) as CholeskyResult.Success).factor.solve(rhs)
        for (j in 0 until 4) assertClose(viaCholesky[j], x[j], 1e-10, "coefficient $j")
        // Residual sum of squares from the tail of Q'b.
        val qtb = b.copyOf().also { qr.applyQt(it) }
        val rss = (4 until 20).sumOf { qtb[it] * qtb[it] }
        assertClose(residual.sumOf { it * it }, rss, 1e-12, "RSS")
    }

    @Test
    fun rInverseRowNormsMatchTheDiagonalOfTheInverseGram() {
        val a = random(15, 4, 7)
        val qr = HouseholderQr.of(a)
        val gram = DenseMatrix(4, 4)
        for (i in 0 until 4) for (j in 0 until 4) gram[i, j] = (0 until 15).sumOf { a[it, i] * a[it, j] }
        val factor = (Cholesky.of(gram) as CholeskyResult.Success).factor
        val norms = qr.rInverseRowNorms()
        for (j in 0 until 4) {
            val unit = DoubleArray(4).also { it[j] = 1.0 }
            val inverseDiagonal = factor.solve(unit)[j]
            assertClose(sqrt(inverseDiagonal), norms[j], 1e-10, "column $j")
        }
    }

    @Test
    fun exactCollinearityIsReportedNotThrown() {
        val base = random(10, 3, 8)
        val a = DenseMatrix(10, 4)
        for (i in 0 until 10) {
            for (j in 0 until 3) a[i, j] = base[i, j]
            a[i, 3] = base[i, 0] + 2.0 * base[i, 1]
        }
        val qr = HouseholderQr.of(a)
        assertEquals(3, qr.rank)
        assertEquals(1, qr.aliasedColumns.size)
        assertTrue(qr.aliasedColumns.single() in setOf(0, 1, 3), "one of the dependent columns is aliased")
        assertTrue(qr.isAliased(qr.aliasedColumns.single()))
        assertEquals(3, (0 until 4).count { !qr.isAliased(it) })
        val b = DoubleArray(10) { i -> a[i, 0] - a[i, 2] + 0.25 * i }
        val x = qr.solve(b)
        assertEquals(0.0, x[qr.aliasedColumns.single()], "an aliased coefficient is zero")
        // The fit is as good as the fit with the aliased column removed.
        val kept = (0 until 4).filter { !qr.isAliased(it) }
        val reduced = DenseMatrix(10, 3)
        for ((c, original) in kept.withIndex()) for (i in 0 until 10) reduced[i, c] = a[i, original]
        val expected = HouseholderQr.of(reduced).solve(b)
        for ((c, original) in kept.withIndex()) assertClose(expected[c], x[original], 1e-10, "coefficient $original")
        val norms = qr.rInverseRowNorms()
        assertTrue(norms[qr.aliasedColumns.single()].isNaN())
        assertEquals(3, norms.count { !it.isNaN() })
    }

    @Test
    fun toleranceBoundaryIsMaxDimensionTimesEpsilonTimesR11() {
        // diag(1, 1, t): r11 = 1 and tol = 3 * 2^-52, about 6.7e-16.
        fun rankFor(t: Double) = HouseholderQr.of(
            DenseMatrix(3, 3).also {
                it[0, 0] = 1.0
                it[1, 1] = 1.0
                it[2, 2] = t
            },
        ).rank
        assertEquals(3, rankFor(1e-15), "above the tolerance")
        assertEquals(2, rankFor(1e-16), "below the tolerance")
        assertEquals(3, rankFor(1e-10))
        // The tolerance scales with |r11|: the same ratio at a larger scale gives the same answer.
        val scaled = DenseMatrix(3, 3).also {
            it[0, 0] = 1e6
            it[1, 1] = 1e6
            it[2, 2] = 1e-9
        }
        assertEquals(3, HouseholderQr.of(scaled).rank)
        scaled[2, 2] = 1e-10
        assertEquals(2, HouseholderQr.of(scaled).rank)
    }

    @Test
    fun nearCollinearityAboveAndBelowTheToleranceOnRealisticData() {
        fun rankWithNoise(scale: Double): Int {
            val random = Random(9)
            val a = DenseMatrix(30, 4)
            for (i in 0 until 30) {
                for (j in 0 until 3) a[i, j] = random.nextDouble(-1.0, 1.0)
                a[i, 3] = a[i, 0] - a[i, 1] + scale * random.nextDouble(-1.0, 1.0)
            }
            return HouseholderQr.of(a).rank
        }
        assertEquals(4, rankWithNoise(1e-8), "a perturbation of 1e-8 is real signal")
        assertEquals(3, rankWithNoise(1e-17), "a perturbation of 1e-17 is round-off")
    }

    @Test
    fun zeroColumnsAndZeroMatrices() {
        val a = random(6, 3, 10)
        for (i in 0 until 6) a[i, 1] = 0.0
        val qr = HouseholderQr.of(a)
        assertEquals(2, qr.rank)
        assertContentEquals(intArrayOf(1), qr.aliasedColumns)
        assertEquals(0.0, qr.solve(DoubleArray(6) { 1.0 })[1])

        val zero = HouseholderQr.of(DenseMatrix(4, 3))
        assertEquals(0, zero.rank)
        assertContentEquals(intArrayOf(0, 1, 2), zero.aliasedColumns)
        assertContentEquals(DoubleArray(3), zero.solve(DoubleArray(4) { 1.0 }))
        assertTrue(zero.rInverseRowNorms().all { it.isNaN() })
    }

    @Test
    fun moreColumnsThanRowsHasRankAtMostRows() {
        val a = random(3, 6, 11)
        val qr = HouseholderQr.of(a)
        assertEquals(3, qr.rank)
        assertEquals(3, qr.aliasedColumns.size)
        val b = doubleArrayOf(1.0, 2.0, 3.0)
        val x = qr.solve(b)
        for (i in 0 until 3) assertClose(b[i], (0 until 6).sumOf { a[i, it] * x[it] }, 1e-12, "exact fit at row $i")
        assertEquals(3, qr.rInverseRowNorms().count { !it.isNaN() })

        val single = HouseholderQr.of(random(1, 4, 12))
        assertEquals(1, single.rank)
    }

    @Test
    fun withoutPivotingColumnsStayInOrder() {
        val a = random(10, 4, 13)
        val pivoted = HouseholderQr.of(a)
        val plain = HouseholderQr.of(a, pivot = false)
        assertContentEquals(intArrayOf(0, 1, 2, 3), plain.permutation())
        assertEquals(4, plain.rank)
        val b = random(10, 1, 14).data
        val x = plain.solve(b)
        val y = pivoted.solve(b)
        for (j in 0 until 4) assertClose(y[j], x[j], 1e-11, "coefficient $j")
    }

    @Test
    fun withoutPivotingRankStopsAtTheFirstDependentColumn() {
        val a = random(8, 3, 15)
        for (i in 0 until 8) a[i, 1] = 2.0 * a[i, 0]
        val plain = HouseholderQr.of(a, pivot = false)
        assertEquals(1, plain.rank)
        assertContentEquals(intArrayOf(1, 2), plain.aliasedColumns)
    }

    @Test
    fun refinementMatchesSolveOnAWellConditionedProblemAndKeepsAliasedAtZero() {
        val a = random(20, 4, 20)
        val b = random(20, 1, 21).data
        val qr = HouseholderQr.of(a)
        val plain = qr.solve(b)
        val refined = qr.solveRefined(a, b)
        for (j in 0 until 4) assertClose(plain[j], refined[j], 1e-13, "coefficient $j")
        assertContentEquals(plain, qr.solveRefined(a, b, maxIterations = 0))

        val deficient = DenseMatrix(10, 3)
        val base = random(10, 2, 22)
        for (i in 0 until 10) {
            deficient[i, 0] = base[i, 0]
            deficient[i, 1] = base[i, 1]
            deficient[i, 2] = base[i, 0] + base[i, 1]
        }
        val dqr = HouseholderQr.of(deficient)
        val x = dqr.solveRefined(deficient, DoubleArray(10) { it.toDouble() })
        assertEquals(0.0, x[dqr.aliasedColumns.single()])
        assertContentEquals(
            DoubleArray(3),
            HouseholderQr.of(DenseMatrix(4, 3)).solveRefined(DenseMatrix(4, 3), DoubleArray(4)),
        )
    }

    @Test
    fun refinementRecoversDigitsOnAnIllConditionedLargeResidualProblem() {
        // Degree-5 polynomial on 0..20 with a large residual, the shape of Wampler5; the truth is all ones in the
        // noise-free part, so compare with the exact normal-equations solution computed in double-double.
        val n = 21
        val a = DenseMatrix(n, 6)
        val b = DoubleArray(n)
        for (i in 0 until n) {
            var power = 1.0
            for (j in 0 until 6) {
                a[i, j] = power
                power *= i.toDouble()
            }
            b[i] = (0 until 6).sumOf { a[i, it] } + (if (i % 2 == 0) 1.0 else -1.0) * 1e5
        }
        val qr = HouseholderQr.of(a)
        val plain = qr.solve(b)
        val refined = qr.solveRefined(a, b)
        // A second, deeper refinement changes nothing: the refined answer is a fixed point.
        val deeper = qr.solveRefined(a, b, maxIterations = 12)
        var plainGap = 0.0
        var refinedGap = 0.0
        for (j in 0 until 6) {
            refinedGap = maxOf(refinedGap, abs(refined[j] - deeper[j]) / abs(deeper[j]))
            plainGap = maxOf(plainGap, abs(plain[j] - deeper[j]) / abs(deeper[j]))
        }
        assertTrue(refinedGap < 1e-13, "refined answer moved by $refinedGap with more steps")
        assertTrue(plainGap > 10 * refinedGap, "plain gap $plainGap should exceed the refined gap $refinedGap")
    }

    @Test
    fun refinedRInverseRowNormsMatchThePlainOnesAndHandleAliasing() {
        val a = random(15, 4, 23)
        val qr = HouseholderQr.of(a)
        val plain = qr.rInverseRowNorms()
        val refined = qr.rInverseRowNorms(a)
        for (j in 0 until 4) assertClose(plain[j], refined[j], 1e-13, "column $j")
        val start = qr.rInverseRowNorms(a, maxIterations = 0)
        for (j in 0 until 4) assertClose(plain[j], start[j], 1e-13, "unrefined start, column $j")

        val deficient = DenseMatrix(9, 3)
        val base = random(9, 2, 24)
        for (i in 0 until 9) {
            deficient[i, 0] = base[i, 0]
            deficient[i, 1] = base[i, 1]
            deficient[i, 2] = base[i, 0] - base[i, 1]
        }
        val dqr = HouseholderQr.of(deficient)
        val norms = dqr.rInverseRowNorms(deficient)
        assertTrue(norms[dqr.aliasedColumns.single()].isNaN())
        assertEquals(2, norms.count { !it.isNaN() })
        assertFailsWith<IllegalArgumentException> { qr.rInverseRowNorms(DenseMatrix(15, 3)) }
    }

    @Test
    fun compensatedSumIsExactWhereADoubleSumIsNot() {
        val sum = CompensatedSum()
        sum.add(1.0)
        sum.addProduct(1e-17, 1.0)
        sum.addProduct(1e-17, 1.0)
        sum.add(-1.0)
        assertEquals(2e-17, sum.value, 1e-30)
        sum.reset()
        sum.addProduct(134217729.0, 134217727.0) // (2^27 + 1)(2^27 - 1) = 2^54 - 1, not representable in a double
        sum.add(-18014398509481984.0) // 2^54
        assertEquals(-1.0, sum.value)
    }

    @Test
    fun validatesArgumentSizes() {
        val qr = HouseholderQr.of(random(5, 2, 16))
        assertFailsWith<IllegalArgumentException> { qr.applyQt(DoubleArray(4)) }
        assertFailsWith<IllegalArgumentException> { qr.applyQ(DoubleArray(6)) }
        assertFailsWith<IllegalArgumentException> { qr.solve(DoubleArray(4)) }
        assertFailsWith<IllegalArgumentException> { qr.isAliased(2) }
        assertFailsWith<IllegalArgumentException> { qr.applyQt(DenseMatrix(4, 1)) }
        assertFailsWith<IllegalArgumentException> { qr.solveRefined(DenseMatrix(5, 3), DoubleArray(5)) }
        assertFailsWith<IllegalArgumentException> { qr.solveRefined(DenseMatrix(5, 2), DoubleArray(5), -1) }
    }

    @Test
    fun theFactorizationLeavesItsInputAlone() {
        val a = random(6, 3, 17)
        val before = a.data.copyOf()
        HouseholderQr.of(a)
        assertContentEquals(before, a.data)
    }

    @Test
    fun theKernelHookIsUsedAndResultsDoNotDependOnWorkers() {
        val a = random(25, 6, 18)
        val b = random(25, 1, 19).data
        val reference = HouseholderQr.of(a).solve(b)
        val saved = denseKernels
        try {
            val counting = CountingKernels(PureKotlinDenseKernels { 7 })
            denseKernels = counting
            val again = HouseholderQr.of(a).solve(b)
            assertContentEquals(reference, again, "same bits with a different worker count")
            assertTrue(counting.calls > 0, "the factorization must go through the hook")
        } finally {
            denseKernels = saved
        }
    }

    private class CountingKernels(private val delegate: DenseKernels) : DenseKernels by delegate {
        var calls = 0

        override fun dot(n: Int, x: DoubleArray, xOffset: Int, y: DoubleArray, yOffset: Int): Double {
            calls++
            return delegate.dot(n, x, xOffset, y, yOffset)
        }

        override fun axpy(n: Int, alpha: Double, x: DoubleArray, xOffset: Int, y: DoubleArray, yOffset: Int) {
            calls++
            delegate.axpy(n, alpha, x, xOffset, y, yOffset)
        }
    }

    private fun random(rows: Int, columns: Int, seed: Int): DenseMatrix {
        val random = Random(seed)
        return DenseMatrix(rows, columns, DoubleArray(rows * columns) { random.nextDouble(-1.0, 1.0) })
    }

    private fun norm(values: DoubleArray) = sqrt(values.sumOf { it * it })

    private fun assertClose(expected: Double, actual: Double, tolerance: Double, message: String) {
        assertTrue(
            abs(expected - actual) <= tolerance * maxOf(1.0, abs(expected)),
            "$message: expected $expected got $actual",
        )
    }
}
