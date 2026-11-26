package net.maizegenetics.tassel.stats.model

import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.testing.OracleManifest
import net.maizegenetics.tassel.testing.TestData
import kotlin.math.abs
import kotlin.math.max
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * The model layer against R's `lm()`, `anova()`, and `drop1()` (`testdata/oracle/r/`, ADR-0029 kind 2): coefficients,
 * standard errors, t, F, sums of squares, and `log10 p` at the manifest's tolerance (`rel 1e-10`), degrees of
 * freedom exactly. Rows are matched by position for coefficients (R and TASSEL 6 name sum-coded columns differently)
 * and by term name for the tables.
 */
class RLinearModelTest {
    @Test fun covariatesOnly() = check("lm-covariates")

    @Test fun factorsOnly() = check("lm-factors")

    @Test fun mixedTreatmentCoding() = check("lm-mixed")

    @Test fun mixedSumCoding() = check("lm-mixed-sum")

    @Test fun strongSignalKeepsLogPValuesFinite() {
        check("lm-strong")
        val strong = load("lm-strong").anova().row("x1")!!
        assertTrue(strong.p < 1e-100 && strong.log10P < -100.0, "p=${strong.p} log10p=${strong.log10P}")
    }

    @Test fun collinearCovariateIsAliasedNotThrown() {
        val fixture = load("lm-collinear")
        val design = fixture.design()
        assertEquals(design.columnCount - 1, design.rank)
        assertEquals(listOf("x3"), design.aliasedColumnNames)
        assertEquals(listOf("x3"), design.aliasedTerms)
        val fit = LinearModelFit.of(design)
        assertTrue(fit.coefficients[design.columns.indexOfFirst { it.name == "x3" }].isNaN())
        assertEquals(fixture.rows - design.rank, fit.residualDf)
        check("lm-collinear")
    }

    @Test fun emptyFactorLevelAfterRowSelection() {
        val fixture = load("lm-empty-level")
        val design = fixture.design()
        assertTrue(design.rows < fixture.rows, "rows with keep = 0 are left out")
        // Levels C and D are gone from the fitted rows: g1 has two levels (A, B) so one column.
        assertEquals(1, design.termColumns(design.spec.terms.indexOfFirst { it.label == "g1" }).size)
        check("lm-empty-level")
    }

    @Test fun unbalancedReplicatesSequentialAndMarginalDiffer() {
        check("lm-unbalanced")
        val design = load("lm-unbalanced").design()
        val sequential = Anova.sequential(design).row("g1")!!.sumOfSquares
        val marginal = Anova.marginal(design).row("g1")!!.sumOfSquares
        assertTrue(abs(sequential - marginal) > 1.0, "sequential $sequential and marginal $marginal should differ")
    }

    @Test fun weightedLeastSquaresIsGlsWithDiagonalCovariance() = check("lm-weighted")

    @Test fun fixturesAreListedInTheManifest() {
        val outputs = OracleManifest.load().outputs.map { it.key }
        for (name in MODELS) {
            for (part in listOf("data.tsv", "model.txt", "coefficients.tsv", "anova.tsv")) {
                assertTrue("oracle/r/$name.$part" in outputs, "$name.$part is in the manifest")
            }
        }
    }

    private fun check(name: String) {
        val fixture = load(name)
        val design = fixture.design()
        val fit = fixture.fit(design)
        val manifest = OracleManifest.load()
        val tolerance = (
            manifest.tolerance(
                "oracle/r/$name.coefficients.tsv",
            ) as net.maizegenetics.tassel.testing.Tolerance.Relative
            ).value
        assertEquals(1e-10, tolerance)

        // Coefficients, by position.
        val coefficients = table(name, "coefficients")
        assertEquals(design.columnCount, coefficients.size, "$name: coefficient count")
        coefficients.forEachIndexed { index, row ->
            val label = "$name coefficient ${row.getValue("term")}"
            same(label + " estimate", row.getValue("estimate"), fit.coefficients[index], tolerance)
            same(label + " se", row.getValue("se"), fit.standardErrors[index], tolerance)
            val t = TTest.coefficient(fit, index)
            same(label + " t", row.getValue("t"), t.statistic, tolerance)
            same(label + " p", row.getValue("p"), t.p, tolerance)
            same(label + " log10p", row.getValue("log10p"), t.log10P, tolerance)
        }

        // Sequential table.
        val anova = fixture.anova(design)
        for (row in table(name, "anova")) {
            val term = row.getValue("term")
            if (term == "Residuals") {
                assertEquals(row.getValue("df").toInt(), fit.residualDf, "$name residual df")
                same("$name residual SS", row.getValue("ss"), fit.residualSumOfSquares, tolerance)
                assertEquals(fit.residualDf, anova.residualDf)
                continue
            }
            compareRow("$name sequential $term", row, anova.row(term), tolerance)
        }
        // A term R leaves out of its table (aliased) has no degrees of freedom here.
        val inR = table(name, "anova").map { it.getValue("term") }.toSet()
        for (row in anova.rows) if (row.term !in inR) assertEquals(0, row.df, "$name: ${row.term} is aliased")

        // Marginal table.
        if ("oracle/r/$name.drop1.tsv" in OracleManifest.load().outputs.map { it.key }) {
            val marginal = fixture.marginal(design)
            for (row in table(name, "drop1")) {
                val term = row.getValue("term")
                if (term == "<none>") {
                    same("$name drop1 RSS", row.getValue("rss"), fit.residualSumOfSquares, tolerance)
                } else {
                    compareRow("$name marginal $term", row, marginal.row(term), tolerance)
                }
            }
        }
    }

    private fun compareRow(label: String, row: Map<String, String>, actual: AnovaRow?, tolerance: Double) {
        assertTrue(actual != null, "$label: no row")
        assertEquals(row.getValue("df").toInt(), actual.df, "$label df")
        same("$label SS", row.getValue("ss"), actual.sumOfSquares, tolerance)
        same("$label F", row.getValue("f"), actual.f, tolerance)
        same("$label p", row.getValue("p"), actual.p, tolerance)
        same("$label log10p", row.getValue("log10p"), actual.log10P, tolerance)
    }

    private fun same(label: String, expected: String, actual: Double, tolerance: Double) {
        if (expected == "NA") {
            assertTrue(actual.isNaN(), "$label: R has NA, got $actual")
            return
        }
        val e = expected.toDouble()
        val scale = max(abs(e), 1e-300)
        assertTrue(abs(actual - e) <= tolerance * scale, "$label: R $e, got $actual")
    }

    private fun table(name: String, part: String): List<Map<String, String>> =
        rows(TestData.readOracleText("r/$name.$part.tsv"))

    private fun rows(text: String): List<Map<String, String>> {
        val lines = text.lines().filter { it.isNotEmpty() }
        val header = lines.first().split('\t')
        return lines.drop(1).map { line -> header.zip(line.split('\t')).toMap() }
    }

    private fun load(name: String): Fixture {
        val data = rows(TestData.readOracleText("r/$name.data.tsv"))
        val model = TestData.readOracleText("r/$name.model.txt").lines().filter { " = " in it }
            .associate { it.substringBefore(" = ") to it.substringAfter(" = ") }
        return Fixture(name, data, model)
    }

    private class Fixture(val name: String, val data: List<Map<String, String>>, val model: Map<String, String>) {
        val rows: Int get() = data.size
        private val formula = model.getValue("formula")
        private val response = formula.substringBefore("~").trim()
        private val variables = formula.substringAfter("~").split('+').map { it.trim() }
        private val contrast = Contrast.valueOf(model.getValue("contrast"))
        private val weights: DoubleArray? =
            if (model.getValue("weights") == "none") null else DoubleArray(rows) { data[it].getValue("w").toDouble() }

        private val phenotype: PhenotypeTable
        val spec: ModelSpec

        init {
            val attributes =
                mutableListOf<net.maizegenetics.tassel.phenotype.PhenotypeAttribute>(
                    PhenotypeTable.taxa(
                        List(rows) {
                            "T$it"
                        },
                    ),
                )
            val terms = mutableListOf<ModelTerm>(ModelTerm.Intercept)
            for (column in listOf(response) + variables) {
                val values = data.map { it.getValue(column) }
                if (values.all { it.toDoubleOrNull() != null }) {
                    val type = if (column == response) AttributeType.data else AttributeType.covariate
                    attributes += NumericAttribute(column, type, FloatArray(rows) { values[it].toDouble().toFloat() })
                    if (column != response) terms += ModelTerm.Covariate(column)
                } else {
                    attributes += FactorAttribute.of(column, values)
                    terms += ModelTerm.Factor(column, contrast)
                }
            }
            phenotype = PhenotypeTable(attributes)
            spec = ModelSpec(response, terms)
        }

        fun design(): DesignMatrix {
            val keep = (0 until rows).filter { data[it].getValue("keep") == "1" }.toIntArray()
            return DesignMatrix.build(spec, phenotype, keep)
        }

        private fun gls(design: DesignMatrix): DesignMatrix {
            val w = weights ?: return design
            val kept = design.observations
            return GeneralizedLeastSquares.whiten(design, Whitener.diagonal(DoubleArray(kept.size) { w[kept[it]] }))
        }

        fun fit(design: DesignMatrix) = LinearModelFit.of(gls(design))

        fun anova(design: DesignMatrix = design()) = Anova.sequential(gls(design))

        fun marginal(design: DesignMatrix) = Anova.marginal(gls(design))
    }

    private companion object {
        val MODELS = listOf(
            "lm-covariates", "lm-factors", "lm-mixed", "lm-mixed-sum", "lm-strong", "lm-collinear",
            "lm-empty-level", "lm-unbalanced", "lm-weighted",
        )
    }
}
