package net.maizegenetics.tassel.stats.model

import net.maizegenetics.tassel.math.DenseMatrix
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class ModelLayerTest {
    private val y = floatArrayOf(5.1f, 6.0f, 5.5f, 7.9f, 8.2f, 7.4f, 9.8f, Float.NaN, 10.1f, 6.3f, 8.8f, 11.6f)
    private val x = floatArrayOf(1f, 2f, 1.5f, 3f, 3.5f, 2.5f, 4f, 5f, 4.5f, 2.2f, 3.8f, 5.5f)
    private val group = listOf("A", "A", "A", "B", "B", "B", "C", "C", "C", "A", "B", "C")
    private val phenotype = PhenotypeTable(
        listOf<PhenotypeAttribute>(
            PhenotypeTable.taxa(List(12) { "T$it" }),
            NumericAttribute("y", net.maizegenetics.tassel.phenotype.AttributeType.data, y),
            NumericAttribute("x", x),
            FactorAttribute.of("group", group),
        ),
    )

    private fun spec(contrast: Contrast = Contrast.treatment) = ModelSpec(
        "y",
        ModelTerm.Intercept,
        ModelTerm.Covariate("x"),
        ModelTerm.Factor("group", contrast),
    )

    @Test
    fun treatmentCodingDropsTheLastLevelAndSumCodingSubtractsIt() {
        val rows = intArrayOf(0, 3, 6)
        val treatment = DesignMatrix.build(spec(), phenotype, rows)
        assertEquals(listOf("(Intercept)", "x", "group=A", "group=B"), treatment.columns.map { it.name })
        assertContentEquals(doubleArrayOf(1.0, 0.0, 0.0), treatment.matrix.columnCopy(2))
        assertContentEquals(doubleArrayOf(0.0, 1.0, 0.0), treatment.matrix.columnCopy(3))
        val sum = DesignMatrix.build(spec(Contrast.sum), phenotype, rows)
        assertContentEquals(doubleArrayOf(1.0, 0.0, -1.0), sum.matrix.columnCopy(2))
        assertContentEquals(doubleArrayOf(0.0, 1.0, -1.0), sum.matrix.columnCopy(3))
        assertContentEquals(intArrayOf(2, 3), sum.termColumns(2))
    }

    @Test
    fun completeObservationsSkipsMissingValuesAndBuildRefusesThem() {
        val complete = DesignMatrix.completeObservations(spec(), phenotype)
        assertContentEquals(intArrayOf(0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11), complete)
        assertFailsWith<IllegalArgumentException> { DesignMatrix.build(spec(), phenotype, intArrayOf(0, 7)) }
        assertFailsWith<IllegalArgumentException> {
            DesignMatrix.build(ModelSpec("nope", ModelTerm.Intercept), phenotype)
        }
        assertFailsWith<IllegalArgumentException> {
            DesignMatrix.build(ModelSpec("y", ModelTerm.Factor("x")), phenotype, complete)
        }
        assertFailsWith<UnsupportedOperationException> {
            DesignMatrix.build(ModelSpec("y", ModelTerm.Intercept, ModelTerm.MarkerTerm()), phenotype, complete)
        }
        assertFailsWith<IllegalArgumentException> { ModelSpec("y", ModelTerm.Intercept, ModelTerm.Intercept) }
    }

    @Test
    fun fitResidualsAreOrthogonalToTheDesignAndSumsOfSquaresAdd() {
        val design = DesignMatrix.build(spec(), phenotype)
        val fit = LinearModelFit.of(design)
        assertEquals(4, fit.rank)
        assertEquals(design.rows - 4, fit.residualDf)
        for (c in 0 until design.columnCount) {
            var dot = 0.0
            for (r in 0 until design.rows) dot += design.matrix[r, c] * fit.residuals[r]
            assertTrue(abs(dot) < 1e-10, "column $c: $dot")
        }
        for (i in 0 until design.rows) assertEquals(design.response[i], fit.fitted[i] + fit.residuals[i], 1e-12)
        val table = Anova.sequential(design)
        val total = design.response.sumOf { it * it }
        // With an intercept row dropped, the rows plus the intercept's share plus the residual make up y'y.
        val interceptShare = total - LinearModelFit.of(design.subset(listOf(0))).residualSumOfSquares
        assertEquals(total, interceptShare + table.rows.sumOf { it.sumOfSquares } + table.residualSumOfSquares, 1e-9)
        val refined = LinearModelFit.of(design, refine = true)
        for (j in fit.coefficients.indices) assertEquals(fit.coefficients[j], refined.coefficients[j], 1e-10)
    }

    @Test
    fun nestedFTestEqualsTheSequentialRowAndTheWaldTestOnTheTermColumns() {
        val design = DesignMatrix.build(spec(), phenotype)
        val full = LinearModelFit.of(design)
        val reduced = LinearModelFit.of(design.subset(listOf(0, 1)))
        val nested = FTest.compare(reduced, full)
        val row = Anova.sequential(design).row("group")!!
        assertEquals(2.0, nested.df)
        assertEquals(row.f, nested.statistic, 1e-12 * row.f)
        assertEquals(row.log10P, nested.log10P, 1e-12 * abs(row.log10P))
        // Group is the last term, so its marginal test is the sequential one, and Wald on its two columns agrees.
        val l = DenseMatrix(2, design.columnCount)
        l[0, 2] = 1.0
        l[1, 3] = 1.0
        val wald = WaldTest.test(full, l)
        assertEquals(nested.statistic, wald.statistic, 1e-9 * nested.statistic)
        assertEquals(2.0, wald.df)
        val chi = WaldTest.test(full, l, knownVariance = full.residualVariance)
        assertEquals(TestKind.ChiSquare, chi.kind)
        assertEquals(wald.statistic * 2, chi.statistic, 1e-9 * chi.statistic)
        // A duplicated row of L does not add a degree of freedom.
        val dependent = DenseMatrix.ofRows(doubleArrayOf(0.0, 0.0, 1.0, 0.0), doubleArrayOf(0.0, 0.0, 2.0, 0.0))
        assertEquals(1.0, WaldTest.test(full, dependent).df)
    }

    @Test
    fun tTestOfACoefficientMatchesWaldAndTheContrastFormsAndSquaresToF() {
        val design = DesignMatrix.build(spec(), phenotype)
        val fit = LinearModelFit.of(design)
        for (column in 0 until design.columnCount) {
            val t = TTest.coefficient(fit, column)
            val unit = DoubleArray(design.columnCount).also { it[column] = 1.0 }
            val viaContrast = TTest.contrast(fit, unit)
            assertEquals(t.statistic, viaContrast.statistic, 1e-9 * abs(t.statistic))
            val wald = WaldTest.test(fit, DenseMatrix.ofRows(unit))
            assertEquals(t.statistic * t.statistic, wald.statistic, 1e-9 * wald.statistic)
            assertEquals(t.p, wald.p, 1e-9 * t.p)
        }
        // The difference of two level effects: group A minus group B.
        val difference = TTest.contrast(fit, doubleArrayOf(0.0, 0.0, 1.0, -1.0))
        val sameAsWald = WaldTest.test(fit, DenseMatrix.ofRows(doubleArrayOf(0.0, 0.0, 1.0, -1.0)))
        assertEquals(difference.statistic * difference.statistic, sameAsWald.statistic, 1e-9 * sameAsWald.statistic)
    }

    @Test
    fun aliasedCoefficientsAreNotEstimable() {
        val duplicated = PhenotypeTable(
            listOf<PhenotypeAttribute>(
                PhenotypeTable.taxa(List(6) { "T$it" }),
                NumericAttribute(
                    "y",
                    net.maizegenetics.tassel.phenotype.AttributeType.data,
                    floatArrayOf(1f, 3f, 2f, 5f, 4f, 7f),
                ),
                NumericAttribute("a", floatArrayOf(1f, 2f, 3f, 4f, 5f, 6f)),
                NumericAttribute("b", floatArrayOf(2f, 4f, 6f, 8f, 10f, 12f)),
                NumericAttribute("zero", floatArrayOf(0f, 0f, 0f, 0f, 0f, 0f)),
            ),
        )
        val design = DesignMatrix.build(
            ModelSpec(
                "y",
                ModelTerm.Intercept,
                ModelTerm.Covariate("a"),
                ModelTerm.Covariate("b"),
                ModelTerm.Covariate("zero"),
            ),
            duplicated,
        )
        assertEquals(2, design.rank)
        assertContentEquals(intArrayOf(2, 3), design.aliasedColumns)
        val fit = LinearModelFit.of(design)
        assertTrue(fit.coefficients[2].isNaN() && fit.coefficients[3].isNaN())
        assertTrue(!TTest.coefficient(fit, 2).estimable)
        val l = DenseMatrix.ofRows(doubleArrayOf(0.0, 0.0, 1.0, 0.0))
        assertTrue(!WaldTest.test(fit, l).estimable)
        val rows = Anova.sequential(design).rows
        assertEquals(listOf(1, 0, 0), rows.map { it.df })
    }

    @Test
    fun fullCovarianceWhiteningMatchesDiagonalWeights() {
        val design = DesignMatrix.build(spec(), phenotype)
        val weights = DoubleArray(design.rows) { 1.0 + 0.25 * it }
        val covariance = DenseMatrix(design.rows, design.rows)
        for (i in 0 until design.rows) covariance[i, i] = 1.0 / weights[i]
        val full = GeneralizedLeastSquares.fit(design, Whitener.fromCovariance(covariance))
        val diagonal = GeneralizedLeastSquares.fit(design, Whitener.diagonal(weights))
        for (j in full.coefficients.indices) {
            assertEquals(diagonal.coefficients[j], full.coefficients[j], 1e-12 * abs(diagonal.coefficients[j]) + 1e-14)
            assertEquals(diagonal.standardErrors[j], full.standardErrors[j], 1e-12 * diagonal.standardErrors[j])
        }
        covariance[0, 0] = -1.0
        assertFailsWith<IllegalArgumentException> { Whitener.fromCovariance(covariance) }
    }
}
