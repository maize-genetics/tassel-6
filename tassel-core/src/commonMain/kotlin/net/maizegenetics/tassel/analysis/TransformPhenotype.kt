package net.maizegenetics.tassel.analysis

import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.operation
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import kotlin.math.ln
import kotlin.math.pow
import kotlin.math.sqrt

/** The base of a logarithm transform. The names are TASSEL 5's, so `-log base_10` translates unchanged. */
@Suppress("ktlint:standard:enum-entry-name-case")
enum class LogBase { natural, base_2, base_10 }

/**
 * What [transformPhenotype] does to each chosen trait, in TASSEL 5's order: a power or a logarithm, then
 * standardization.
 *
 * @property traits the data attributes to transform; empty means every data attribute. Covariates and factors are
 *   never transformed, as in headless TASSEL 5.
 * @property factors factor attributes whose level combinations are the groups to standardize within; empty means one
 *   group per trait. Needs [standardize].
 * @property log a logarithm to take, or `null`. Ignored when [power] is given: power wins, as in TASSEL 5.
 * @property power an exponent to raise each value to, or `null`.
 * @property standardize subtracts the mean and divides by the n - 1 standard deviation, skipping missing values.
 */
data class PhenotypeTransform(
    val traits: List<String> = emptyList(),
    val factors: List<String> = emptyList(),
    val log: LogBase? = null,
    val power: Double? = null,
    val standardize: Boolean = false,
) {
    init {
        require(log != null || power != null || standardize) { "Nothing to do: give a log, a power, or standardize" }
        require(power == null || !power.isNaN()) { "power is not a number" }
        require(factors.isEmpty() || standardize) { "factors group the standardization, so standardize must be true" }
    }
}

/** The transformed table, and one warning per trait that lost values. */
class PhenotypeTransformResult(val table: PhenotypeTable, val warnings: List<String>)

/**
 * Transforms the data attributes of [table], returning a new table with every attribute in its place.
 *
 * Values are widened to `Double` from the stored floats, transformed, and narrowed to `Float`; a standardization
 * works from the narrowed power or logarithm. A missing value stays missing. Any other value whose result is NaN or
 * infinite, such as the logarithm of zero or of a negative number, becomes missing, and the result's warnings count
 * them per trait. A standardization has no result when its group has one value or a standard deviation of zero.
 * A row with a missing factor value belongs to no group, so it becomes missing when standardizing by factors.
 *
 * TASSEL 5 instead keeps a NaN it did not flag and a `-Infinity`, standardizes with a `Float` mean, and leaves a row
 * with a missing factor unstandardized (ADR-0023).
 *
 * @throws IllegalArgumentException for a trait that is not a data attribute, or a factor that is not a factor.
 * @see net.maizegenetics.analysis.numericaltransform.TransformDataPlugin
 */
suspend fun transformPhenotype(
    table: PhenotypeTable,
    transform: PhenotypeTransform,
    progress: ProgressReporter = ProgressReporter.None,
): PhenotypeTransformResult {
    val chosen = if (transform.traits.isEmpty()) {
        table.attributes.filter { it.type == AttributeType.data }.map { it.name }
    } else {
        transform.traits.distinct().onEach { name ->
            val attribute = table.attribute(name)
            require(attribute.type == AttributeType.data) {
                "'$name' is ${attribute.type}; only data attributes can be transformed"
            }
        }
    }
    val factors = transform.factors.distinct().map { name ->
        val attribute = table.attribute(name)
        require(attribute is FactorAttribute) { "'$name' is ${attribute.type}, not a factor" }
        attribute
    }
    val groups = if (transform.standardize) groupsOf(table.numberOfObservations, factors) else emptyList()

    val warnings = mutableListOf<String>()
    var done = 0
    val attributes = table.attributes.map { attribute ->
        if (attribute.name !in chosen) return@map attribute
        val result = transformTrait(attribute as NumericAttribute, transform, groups, factors.isNotEmpty())
        result.warning(attribute.name)?.let(warnings::add)
        progress.report(++done, chosen.size, "Transforming ${attribute.name}")
        result.attribute
    }
    return PhenotypeTransformResult(PhenotypeTable(attributes, "transformed_${table.name}"), warnings)
}

/** The observation indices of each group of identical factor levels; rows missing a factor are in none. */
private fun groupsOf(rows: Int, factors: List<FactorAttribute>): List<IntArray> {
    if (factors.isEmpty()) return listOf(IntArray(rows) { it })
    val groups = LinkedHashMap<List<Int>, MutableList<Int>>()
    for (observation in 0 until rows) {
        val codes = factors.map { it.code(observation) }
        if (codes.any { it < 0 }) continue
        groups.getOrPut(codes) { mutableListOf() }.add(observation)
    }
    return groups.values.map { it.toIntArray() }
}

private class TraitResult(val attribute: NumericAttribute, val lost: Int, val ungrouped: Int) {
    fun warning(name: String): String? {
        if (lost == 0 && ungrouped == 0) return null
        val parts = buildList {
            if (lost > 0) add("$lost values had no finite result")
            if (ungrouped > 0) add("$ungrouped have a missing factor value")
        }
        return "Trait '$name': ${parts.joinToString(" and ")}, so they are missing"
    }
}

private fun transformTrait(
    attribute: NumericAttribute,
    transform: PhenotypeTransform,
    groups: List<IntArray>,
    byFactors: Boolean,
): TraitResult {
    val size = attribute.size
    var values = DoubleArray(size) { attribute.value(it) }
    for (observation in 0 until size) {
        if (values[observation].isNaN()) continue
        val x = values[observation]
        values[observation] = when {
            transform.power != null -> narrow(x.pow(transform.power))
            transform.log == LogBase.natural -> narrow(ln(x))
            transform.log == LogBase.base_2 -> narrow(ln(x) / ln(2.0))
            transform.log == LogBase.base_10 -> narrow(ln(x) / ln(10.0))
            else -> x
        }
    }
    var ungrouped = 0
    if (transform.standardize) {
        val standardized = DoubleArray(size) { Double.NaN }
        val grouped = BooleanArray(size)
        for (group in groups) {
            val (mean, sd) = meanAndSd(values, group)
            for (observation in group) {
                grouped[observation] = true
                standardized[observation] = narrow((values[observation] - mean) / sd)
            }
        }
        if (byFactors) ungrouped = (0 until size).count { !grouped[it] && !values[it].isNaN() }
        values = standardized
    }
    val lost = (0 until size).count { !attribute.isMissing(it) && values[it].isNaN() } - ungrouped
    val result = NumericAttribute(attribute.name, attribute.type, FloatArray(size) { values[it].toFloat() })
    return TraitResult(result, lost, ungrouped)
}

/** [value] rounded to a `Float`, as the stored result is; NaN when that is not finite. */
private fun narrow(value: Double): Double = value.toFloat().let { if (it.isFinite()) it.toDouble() else Double.NaN }

/** The mean and n - 1 standard deviation of the non-missing [values] at [indices]; NaN with fewer than two. */
private fun meanAndSd(values: DoubleArray, indices: IntArray): Pair<Double, Double> {
    var n = 0
    var sum = 0.0
    for (i in indices) {
        if (!values[i].isNaN()) {
            sum += values[i]
            n++
        }
    }
    if (n == 0) return Double.NaN to Double.NaN
    val mean = sum / n
    var squares = 0.0
    for (i in indices) if (!values[i].isNaN()) squares += (values[i] - mean) * (values[i] - mean)
    return mean to if (n < 2) Double.NaN else sqrt(squares / (n - 1))
}

/**
 * The spec for [transformPhenotype]. Declared here; `TasselOperations` registers it (#166).
 *
 * The legacy names are those of TASSEL 5's `TransformDataPlugin`; `traits` and `factor` are comma separated.
 *
 * @see net.maizegenetics.analysis.numericaltransform.TransformDataPlugin
 */
val transformPhenotypeOperation = operation<PhenotypeTable, PhenotypeTransformResult>(
    id = "transform-phenotype",
    legacyPlugin = "TransformDataPlugin",
    help = "Takes a logarithm or a power of phenotype traits and standardizes them, optionally within factor levels",
) {
    val traits by stringList(legacy = "traits", help = "Data attributes to transform, comma separated; default all")
    val factors by stringList(legacy = "factor", help = "Factors to standardize within, comma separated")
    val log by enum<LogBase>(legacy = "log", help = "Logarithm to take: natural, base_2, or base_10")
    val power by double(legacy = "power", help = "Exponent to raise each value to; wins over log")
    val standardize by boolean(default = false, legacy = "standardize", help = "Subtract the mean, divide by the SD")

    run { input, progress ->
        transformPhenotype(
            input,
            PhenotypeTransform(
                traits = traits().orEmpty(),
                factors = factors().orEmpty(),
                log = log(),
                power = power(),
                standardize = standardize(),
            ),
            progress,
        )
    }
}
