package net.maizegenetics.tassel.analysis

import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.operation
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.table.ResultTable
import net.maizegenetics.tassel.table.ResultTableBuilder
import kotlin.math.sqrt

/**
 * The tables [summarizePhenotype] produced.
 *
 * Each is a plain [ResultTable], so the CLI, the app, and the exporters show it without knowing it came from a
 * phenotype. TASSEL 5 has no phenotype summary, so nothing here has a golden.
 */
class PhenotypeSummary internal constructor(
    /**
     * One row per attribute other than the taxa attribute, titled `Attribute Summary`. Columns: `Attribute`, `Type`,
     * `Observations`, `Missing`, `Proportion Missing`, `Mean`, `Std Dev`, `Min`, `Max`, `Median`, `Levels`.
     * `Observations` counts the non-missing values, so `Observations + Missing` is the number of rows. The moment
     * columns are empty for a factor, and `Levels` is empty for a numeric attribute.
     */
    val attributes: ResultTable,
    /** One row per factor level, then one `?` row for missing values, per factor; titled `Factor Levels`. */
    val factorLevels: ResultTable,
    /** One row per taxon in first-appearance order, titled `Taxa Summary`. */
    val taxa: ResultTable,
) {
    /** The tables, in the order attributes, factor levels, taxa. */
    val tables: List<ResultTable> = listOf(attributes, factorLevels, taxa)
}

private val ATTRIBUTE_COLUMNS = listOf(
    "Attribute", "Type", "Observations", "Missing", "Proportion Missing",
    "Mean", "Std Dev", "Min", "Max", "Median", "Levels",
)

/**
 * Summarizes [table] as per-attribute counts and moments, factor level counts, and per-taxon counts.
 *
 * Moments are computed in `Double` from the stored `Float`s, skipping missing values: the mean first, then the sum of
 * squared deviations from it, divided by n - 1 as TASSEL 5's `TransformDataPlugin` does. A numeric attribute with no
 * values has NaN moments, minimum, maximum, and median; with one value its standard deviation is NaN. The median of an
 * even count is the mean of the two middle values. Proportions with no rows are NaN.
 *
 * Every declared factor level is listed, including one no observation uses (count 0).
 *
 * @see net.maizegenetics.analysis.numericaltransform.TransformDataPlugin
 * @see net.maizegenetics.analysis.data.GenotypeSummaryPlugin
 */
suspend fun summarizePhenotype(
    table: PhenotypeTable,
    progress: ProgressReporter = ProgressReporter.None,
): PhenotypeSummary {
    val columns = table.attributes.filter { it.type != AttributeType.taxa }
    val rows = table.numberOfObservations
    val total = columns.size + 2
    val attributes = ResultTableBuilder("Attribute Summary", ATTRIBUTE_COLUMNS)
    val levels = ResultTableBuilder("Factor Levels", listOf("Factor", "Level", "Count"))
    columns.forEachIndexed { index, attribute ->
        val missing = attribute.missingCount
        val common = listOf(
            attribute.name,
            attribute.type.name,
            rows - missing,
            missing,
            missing.toDouble() / rows.toDouble(),
        )
        when (attribute) {
            is NumericAttribute -> {
                val m = moments(attribute)
                attributes.row(common + listOf(m.mean, m.sd, m.min, m.max, m.median, null))
            }
            is FactorAttribute -> {
                attributes.row(common + listOf(null, null, null, null, null, attribute.levels.size))
                val counts = IntArray(attribute.levels.size)
                for (observation in 0 until rows) {
                    val code = attribute.code(observation)
                    if (code >= 0) counts[code]++
                }
                attribute.levels.forEachIndexed { code, level -> levels.row(attribute.name, level, counts[code]) }
                levels.row(attribute.name, FactorAttribute.MISSING, missing)
            }
            else -> error("Unexpected attribute ${attribute.type}")
        }
        progress.report(index + 1, total, "Summarizing attributes")
    }

    val firstAppearance = LinkedHashMap<String, IntArray>()
    for (observation in 0 until rows) {
        val counts = firstAppearance.getOrPut(table.taxon(observation)) { IntArray(2) }
        counts[0]++
        counts[1] += columns.count { it.isMissing(observation) }
    }
    val taxa = ResultTableBuilder("Taxa Summary", listOf("Taxa", "Observations", "Missing"))
    for ((taxon, counts) in firstAppearance) taxa.row(taxon, counts[0], counts[1])
    progress.report(total, total, "Summarizing taxa")
    return PhenotypeSummary(attributes.build(), levels.build(), taxa.build())
}

private class Moments(val mean: Double, val sd: Double, val min: Double, val max: Double, val median: Double)

private fun moments(attribute: NumericAttribute): Moments {
    val values = DoubleArray(attribute.size - attribute.missingCount)
    var n = 0
    for (observation in 0 until attribute.size) {
        if (!attribute.isMissing(observation)) values[n++] = attribute.value(observation)
    }
    if (n == 0) return Moments(Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN)
    var sum = 0.0
    for (v in values) sum += v
    val mean = sum / n
    var squares = 0.0
    for (v in values) squares += (v - mean) * (v - mean)
    val sd = if (n < 2) Double.NaN else sqrt(squares / (n - 1))
    values.sort()
    val median = if (n % 2 == 1) values[n / 2] else (values[n / 2 - 1] + values[n / 2]) / 2.0
    return Moments(mean, sd, values[0], values[n - 1], median)
}

/**
 * The spec for [summarizePhenotype]. Declared here; `TasselOperations` registers it (#166).
 *
 * @see net.maizegenetics.analysis.numericaltransform.TransformDataPlugin
 */
val phenotypeSummaryOperation = operation<PhenotypeTable, PhenotypeSummary>(
    id = "phenotype-summary",
    help = "Counts, missingness, and moments of each phenotype attribute, factor level counts, and counts by taxon",
) {
    run { input, progress -> summarizePhenotype(input, progress) }
}
