package net.maizegenetics.tassel.analysis

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.table.ResultTable
import kotlin.math.abs
import kotlin.math.sqrt
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/** [summarizePhenotype] against hand-computed values. */
class PhenotypeSummaryTest {
    private val nan = Float.NaN

    private fun numeric(name: String, vararg v: Float, type: AttributeType = AttributeType.data) =
        NumericAttribute(name, type, v)

    private fun table(taxa: List<String>, vararg attrs: net.maizegenetics.tassel.phenotype.PhenotypeAttribute) =
        PhenotypeTable(listOf(PhenotypeTable.taxa(taxa)) + attrs.toList(), "t")

    private fun ResultTable.row(key: String): List<Any?> = rows.first { it[0] == key }

    private fun near(expected: Double, actual: Any?) {
        val a = actual as Double
        assertTrue(abs(a - expected) <= 1e-12 * abs(expected).coerceAtLeast(1e-300), "expected $expected, got $a")
    }

    @Test
    fun replicatesWithMissingValuesAndFactor() = runTest {
        // A has 2 rows, B has 2 rows, C has 1. yield: 1, ?, 4, 7, 10 -> n=4 mean 5.5, median 5.5
        val t = table(
            listOf("A", "A", "B", "B", "C"),
            FactorAttribute.of("env", listOf("e2", "e1", null, "e1", "e1")),
            numeric("yield", 1f, nan, 4f, 7f, 10f),
        )
        val s = summarizePhenotype(t)
        assertEquals(
            listOf(
                "Attribute", "Type", "Observations", "Missing", "Proportion Missing",
                "Mean", "Std Dev", "Min", "Max", "Median", "Levels",
            ),
            s.attributes.columnNames,
        )
        assertEquals(listOf<Any?>("env", "factor", 4, 1, 0.2, null, null, null, null, null, 2), s.attributes.row("env"))
        val y = s.attributes.row("yield")
        assertEquals(listOf<Any?>("yield", "data", 4, 1, 0.2), y.take(5))
        near(5.5, y[5])
        near(sqrt((4.5 * 4.5 + 1.5 * 1.5 + 1.5 * 1.5 + 4.5 * 4.5) / 3.0), y[6])
        assertEquals(1.0, y[7])
        assertEquals(10.0, y[8])
        assertEquals(5.5, y[9])
        assertEquals(null, y[10])
        assertEquals(
            listOf(
                listOf<Any?>("env", "e1", 3),
                listOf<Any?>("env", "e2", 1),
                listOf<Any?>("env", "?", 1),
            ),
            s.factorLevels.rows,
        )
        assertEquals(
            listOf(listOf<Any?>("A", 2, 1), listOf<Any?>("B", 2, 1), listOf<Any?>("C", 1, 0)),
            s.taxa.rows,
        )
    }

    @Test
    fun taxaKeepFirstAppearanceOrder() = runTest {
        val t = table(listOf("Z", "A", "Z"), numeric("x", 1f, nan, nan))
        assertEquals(listOf("Z", "A"), summarizePhenotype(t).taxa.rows.map { it[0] })
        assertEquals(listOf(listOf<Any?>("Z", 2, 1), listOf<Any?>("A", 1, 1)), summarizePhenotype(t).taxa.rows)
    }

    @Test
    fun allMissingColumn() = runTest {
        val row = summarizePhenotype(table(listOf("A", "B"), numeric("x", nan, nan))).attributes.row("x")
        assertEquals(listOf<Any?>("x", "data", 0, 2, 1.0), row.take(5))
        row.subList(5, 10).forEach { assertTrue((it as Double).isNaN()) }
    }

    @Test
    fun singleObservationHasNanSd() = runTest {
        val row = summarizePhenotype(table(listOf("A"), numeric("x", 3f))).attributes.row("x")
        assertEquals(3.0, row[5])
        assertTrue((row[6] as Double).isNaN())
        assertEquals(listOf(3.0, 3.0, 3.0), row.subList(7, 10))
    }

    @Test
    fun medianOddAndEven() = runTest {
        val odd = summarizePhenotype(table(List(5) { "T$it" }, numeric("x", 9f, 1f, 5f, 3f, 7f))).attributes.row("x")
        assertEquals(5.0, odd[9])
        val even = summarizePhenotype(table(List(4) { "T$it" }, numeric("x", 8f, 1f, 2f, 4f))).attributes.row("x")
        assertEquals(3.0, even[9])
    }

    @Test
    fun covariateOnlyTable() = runTest {
        val s = summarizePhenotype(table(listOf("A", "B"), numeric("q", 0.5f, 1.5f, type = AttributeType.covariate)))
        assertEquals(1, s.attributes.rowCount)
        assertEquals("covariate", s.attributes[0, 1])
        assertEquals(0, s.factorLevels.rowCount)
        assertEquals(listOf(listOf<Any?>("A", 1, 0), listOf<Any?>("B", 1, 0)), s.taxa.rows)
    }

    @Test
    fun factorOnlyTableListsUnusedLevelsAndMissing() = runTest {
        val f = FactorAttribute("env", intArrayOf(0, 0, -1), listOf("a", "b"))
        val s = summarizePhenotype(table(listOf("A", "A", "B"), f))
        assertEquals(
            listOf(listOf<Any?>("env", "a", 2), listOf<Any?>("env", "b", 0), listOf<Any?>("env", "?", 1)),
            s.factorLevels.rows,
        )
        assertEquals(listOf<Any?>("env", "factor", 2, 1), s.attributes.rows[0].take(4))
        assertEquals(listOf(listOf<Any?>("A", 2, 0), listOf<Any?>("B", 1, 1)), s.taxa.rows)
    }

    @Test
    fun momentsMatchTwoPassReference() = runTest {
        val values = FloatArray(1000) { (it * 0.37f + (it % 7) * 1.1f) / 3f }
        val row = summarizePhenotype(table(List(1000) { "T$it" }, NumericAttribute("x", AttributeType.data, values)))
            .attributes.row("x")
        val d = values.map { it.toDouble() }
        val mean = d.sum() / d.size
        val sd = sqrt(d.sumOf { (it - mean) * (it - mean) } / (d.size - 1))
        near(mean, row[5])
        near(sd, row[6])
    }
}
