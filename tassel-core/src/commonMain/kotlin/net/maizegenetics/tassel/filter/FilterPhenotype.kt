package net.maizegenetics.tassel.filter

import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.operation
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable

/**
 * Which observations and attributes a [filterPhenotype] call keeps, and which numeric attributes it retypes.
 *
 * The default `PhenotypeFilter()` changes nothing. Every filter is optional and they combine: an observation is kept
 * only when every observation filter keeps it. Observation filters name attributes of the input table, so they
 * work even when the same call drops those attributes. Attribute filters then choose the columns, and retyping comes
 * last, so [keepTypes] sees the input's types.
 *
 * Names are checked against the input table: an unknown taxon, attribute, or factor level is an error that lists the
 * valid ones. The taxa attribute is always kept.
 *
 * @property taxa taxa to keep or to drop, or `null` to start from every observation.
 * @property includeTaxa keeps the observations of [taxa] when true, drops them when false.
 * @property factor the factor attribute that [levels] refers to.
 * @property levels levels of [factor] to keep or to drop. An observation with a missing [factor] value is in no
 *   level, so it is dropped when [includeLevels] is true and kept when it is false.
 * @property includeLevels keeps the observations at [levels] when true, drops them when false.
 * @property attribute the numeric attribute that [min] and [max] bound.
 * @property min smallest value of [attribute] to keep, inclusive. A missing value fails the range.
 * @property max largest value of [attribute] to keep, inclusive. A missing value fails the range.
 * @property requireNonMissing attributes that must have a value; observations missing any of them are dropped.
 * @property requireAllNonMissing requires every attribute of the input table to have a value.
 * @property keepAttributes attributes to keep, or `null` to keep every one. They stay in table order.
 * @property dropAttributes attributes to drop, applied after [keepAttributes]. The taxa attribute cannot be dropped.
 * @property keepTypes attribute types to keep, or `null` for every type. [AttributeType.taxa] is kept regardless.
 * @property toCovariate numeric attributes to turn into [AttributeType.covariate].
 * @property toData numeric attributes to turn into [AttributeType.data].
 * @property excludeLast drops the input's last attribute that is not the taxa attribute.
 * @see net.maizegenetics.phenotype.PhenotypeBuilder
 * @see net.maizegenetics.analysis.filter.FilterTraitsPlugin
 */
data class PhenotypeFilter(
    val taxa: List<String>? = null,
    val includeTaxa: Boolean = true,
    val factor: String? = null,
    val levels: List<String>? = null,
    val includeLevels: Boolean = true,
    val attribute: String? = null,
    val min: Double? = null,
    val max: Double? = null,
    val requireNonMissing: List<String>? = null,
    val requireAllNonMissing: Boolean = false,
    val keepAttributes: List<String>? = null,
    val dropAttributes: List<String>? = null,
    val keepTypes: Set<AttributeType>? = null,
    val toCovariate: List<String>? = null,
    val toData: List<String>? = null,
    val excludeLast: Boolean = false,
) {
    init {
        require((factor == null) == (levels == null)) { "factor and levels must be given together" }
        require(attribute != null || (min == null && max == null)) { "min and max need an attribute" }
        require(attribute == null || min != null || max != null) { "attribute needs a min or a max" }
        require(min == null || !min.isNaN()) { "min is not a number" }
        require(max == null || !max.isNaN()) { "max is not a number" }
        require(min == null || max == null || min <= max) { "min $min is greater than max $max" }
        val both = toCovariate.orEmpty().intersect(toData.orEmpty().toSet())
        require(both.isEmpty()) { "Cannot turn $both into both covariate and data" }
    }
}

/**
 * Filters [table] by [filter], returning a new table. Returns [table] itself when nothing is removed or retyped.
 *
 * Observations keep table order, so replicates stay together and in place, and attributes keep table order. Filtering
 * every observation out is allowed: the result has the same attributes and no rows.
 *
 * TASSEL 5 has `keepTaxa`, `removeTaxa`, `keepAttributes`, `changeAttributeType`, and `removeMissingObservations` on
 * its phenotype builder, and an `excludeLast` switch on `FilterTraitsPlugin`. The factor-level and numeric-range
 * filters are new in TASSEL 6. `excludeLast` differs in one way: TASSEL 5 drops the last column even if it is the
 * taxa column; this drops the last column that is not.
 *
 * @throws IllegalArgumentException for an unknown taxon, attribute, factor, or level, for a non-numeric attribute where
 *   a numeric one is needed, or for naming the taxa attribute where it cannot be dropped or retyped.
 * @see net.maizegenetics.phenotype.PhenotypeBuilder
 * @see net.maizegenetics.phenotype.FilterPhenotype
 * @see net.maizegenetics.analysis.filter.FilterTraitsPlugin
 */
@Suppress("CyclomaticComplexMethod")
suspend fun filterPhenotype(
    table: PhenotypeTable,
    filter: PhenotypeFilter,
    progress: ProgressReporter = ProgressReporter.None,
): PhenotypeTable {
    val rows = keptObservations(table, filter)
    progress.report(1, 3, "Filtering observations")
    val columns = keptAttributes(table, filter)
    progress.report(2, 3, "Filtering attributes")
    val retyped = retypeTargets(table, filter)

    val sameRows = rows.size == table.numberOfObservations
    val sameColumns = columns.size == table.numberOfAttributes
    val retypes = retyped.filter { (name, type) -> name in columns && table.attribute(name).type != type }
    if (sameRows && sameColumns && retypes.isEmpty()) return table

    var result = if (sameRows) table else table.selectRows(rows)
    if (!sameColumns) result = result.selectAttributes(columns)
    if (retypes.isNotEmpty()) {
        result = PhenotypeTable(
            result.attributes.map { attribute ->
                val type = retypes[attribute.name]
                if (type == null || attribute !is NumericAttribute) {
                    attribute
                } else {
                    NumericAttribute(attribute.name, type, attribute.toFloatArray())
                }
            },
            result.name,
        )
    }
    progress.report(3, 3, "Filtering phenotype")
    return result
}

/**
 * The spec for [filterPhenotype]. Declared here; `TasselOperations` registers it (#166).
 *
 * TASSEL 5's plugin has one parameter, `excludeLastTrait`; every other parameter is new. Lists are comma separated.
 * `keepTypes` takes `factor`, `covariate`, `data`, or `taxa`.
 *
 * @see net.maizegenetics.analysis.filter.FilterTraitsPlugin
 */
val filterPhenotypeOperation = operation<PhenotypeTable, PhenotypeTable>(
    id = "filter-phenotype",
    legacyPlugin = "FilterTraitsPlugin",
    help = "Keeps a subset of a phenotype's observations and attributes, and can retype numeric attributes",
) {
    val taxa by stringList(help = "Taxa to keep or to drop, comma separated; by default every observation is kept")
    val includeTaxa by boolean(default = true, help = "Keep the named taxa when true, drop them when false")
    val factor by string(help = "Factor attribute that levels refers to")
    val levels by stringList(help = "Levels of the factor to keep or to drop, comma separated")
    val includeLevels by boolean(default = true, help = "Keep the named levels when true, drop them when false")
    val attribute by string(help = "Numeric attribute that min and max bound")
    val min by double(help = "Smallest value of the attribute to keep, inclusive; missing values fail the range")
    val max by double(help = "Largest value of the attribute to keep, inclusive; missing values fail the range")
    val requireNonMissing by stringList(help = "Attributes that must have a value, comma separated")
    val requireAllNonMissing by boolean(default = false, help = "Require every attribute to have a value")
    val keepAttributes by stringList(help = "Attributes to keep, comma separated; the taxa attribute is always kept")
    val dropAttributes by stringList(help = "Attributes to drop, comma separated")
    val keepTypes by stringList(help = "Attribute types to keep: factor, covariate, data")
    val toCovariate by stringList(help = "Numeric attributes to turn into covariates, comma separated")
    val toData by stringList(help = "Numeric attributes to turn into data, comma separated")
    val excludeLast by boolean(
        default = false,
        legacy = "excludeLastTrait",
        help = "Drop the last attribute",
    )

    run { input, progress ->
        filterPhenotype(
            input,
            PhenotypeFilter(
                taxa = taxa(),
                includeTaxa = includeTaxa(),
                factor = factor(),
                levels = levels(),
                includeLevels = includeLevels(),
                attribute = attribute(),
                min = min(),
                max = max(),
                requireNonMissing = requireNonMissing(),
                requireAllNonMissing = requireAllNonMissing(),
                keepAttributes = keepAttributes(),
                dropAttributes = dropAttributes(),
                keepTypes = keepTypes()?.map(::attributeType)?.toSet(),
                toCovariate = toCovariate(),
                toData = toData(),
                excludeLast = excludeLast(),
            ),
            progress,
        )
    }
}

private fun attributeType(text: String): AttributeType =
    requireNotNull(AttributeType.entries.firstOrNull { it.name == text }) {
        "Unknown attribute type '$text'; the types are ${AttributeType.entries.map { it.name }}"
    }

/** Observation indices that pass every observation filter, in table order. */
private fun keptObservations(table: PhenotypeTable, filter: PhenotypeFilter): IntArray {
    val keep = BooleanArray(table.numberOfObservations) { true }

    filter.taxa?.let { names ->
        val known = table.taxa.names.toSet()
        val unknown = names.filter { it !in known }
        require(unknown.isEmpty()) { "Unknown taxa $unknown; the taxa are ${table.taxa.names}" }
        val named = names.toSet()
        for (observation in keep.indices) {
            if ((table.taxon(observation) in named) != filter.includeTaxa) keep[observation] = false
        }
    }

    if (filter.factor != null && filter.levels != null) {
        val attribute = table.attribute(filter.factor)
        require(attribute is FactorAttribute) { "'${filter.factor}' is ${attribute.type}, not a factor" }
        val unknown = filter.levels.filter { it !in attribute.levels }
        require(unknown.isEmpty()) {
            "Unknown levels $unknown of factor '${filter.factor}'; the levels are ${attribute.levels}"
        }
        val named = filter.levels.toSet()
        for (observation in keep.indices) {
            // A missing value is at no level, so it is dropped by "keep these" and kept by "drop these".
            val inNamed = attribute.level(observation)?.let { it in named } ?: false
            if (inNamed != filter.includeLevels) keep[observation] = false
        }
    }

    if (filter.attribute != null) {
        val attribute = numeric(table, filter.attribute, "bounded by min and max")
        // Values are stored as Float, so bounds are compared as Float: min = 0.1 keeps a stored 0.1f.
        val min = filter.min?.toFloat()
        val max = filter.max?.toFloat()
        for (observation in keep.indices) {
            val value = attribute.floatValue(observation)
            val inRange = !value.isNaN() && (min == null || value >= min) && (max == null || value <= max)
            if (!inRange) keep[observation] = false
        }
    }

    val required = when {
        filter.requireAllNonMissing -> table.attributes
        else -> filter.requireNonMissing.orEmpty().map { table.attribute(it) }
    }
    for (attribute in required) {
        for (observation in keep.indices) if (attribute.isMissing(observation)) keep[observation] = false
    }

    return keep.indices.filter { keep[it] }.toIntArray()
}

/** The names of the attributes that stay, in table order, taxa attribute included. */
private fun keptAttributes(table: PhenotypeTable, filter: PhenotypeFilter): List<String> {
    val taxaName = table.taxaAttribute.name
    var kept: List<PhenotypeAttribute> = table.attributes

    filter.keepAttributes?.let { names ->
        names.forEach { table.attribute(it) }
        kept = kept.filter { it.name in names || it.name == taxaName }
    }
    filter.dropAttributes?.let { names ->
        names.forEach { table.attribute(it) }
        require(taxaName !in names) { "The taxa attribute '$taxaName' cannot be dropped" }
        kept = kept.filter { it.name !in names }
    }
    filter.keepTypes?.let { types ->
        kept = kept.filter { it.type == AttributeType.taxa || it.type in types }
    }
    if (filter.excludeLast) {
        val last = table.attributes.lastOrNull { it.type != AttributeType.taxa }
        requireNotNull(last) { "Nothing to exclude: the only attribute is the taxa attribute '$taxaName'" }
        kept = kept.filter { it !== last }
    }
    return kept.map { it.name }
}

/** The new type of each attribute named by [PhenotypeFilter.toCovariate] or [PhenotypeFilter.toData]. */
private fun retypeTargets(table: PhenotypeTable, filter: PhenotypeFilter): Map<String, AttributeType> {
    val targets = LinkedHashMap<String, AttributeType>()
    filter.toCovariate.orEmpty().forEach { name ->
        numeric(table, name, "retyped")
        targets[name] = AttributeType.covariate
    }
    filter.toData.orEmpty().forEach { name ->
        numeric(table, name, "retyped")
        targets[name] = AttributeType.data
    }
    return targets
}

private fun numeric(table: PhenotypeTable, name: String, what: String): NumericAttribute {
    val attribute = table.attribute(name)
    require(attribute is NumericAttribute) {
        "'$name' is ${attribute.type}; only data and covariate attributes can be $what"
    }
    return attribute
}
