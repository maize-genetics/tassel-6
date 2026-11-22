package net.maizegenetics.tassel.join

import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.operation
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.phenotype.TaxaAttribute

/** Which rows a join keeps. */
@Suppress("ktlint:standard:enum-entry-name-case")
enum class JoinType {
    /** Only matched row pairs (an inner join). */
    intersect,

    /** Matched pairs, unmatched left rows, and unmatched right rows (a full outer join). */
    union,
}

/** What a join does when a key repeats on both sides. */
@Suppress("ktlint:standard:enum-entry-name-case")
enum class JoinMultiple {
    /** Fail, naming the keys. */
    error,

    /** Emit every left-right pair. */
    cartesian,
}

/** What a join does with a non-key column name that both sides have. */
@Suppress("ktlint:standard:enum-entry-name-case")
enum class JoinCollision {
    /** Fail, listing the names. */
    error,

    /** Rename each side's column to `name.1`, `name.2`, ..., numbered by input position. */
    suffix,

    /** Keep one column holding the non-missing value; fail when both are present and differ. */
    coalesce,
}

/** The choices of [joinPhenotypes]. The defaults are the safe ones: a join the data does not allow fails. */
data class JoinOptions(
    val type: JoinType = JoinType.intersect,
    val multiple: JoinMultiple = JoinMultiple.error,
    val collision: JoinCollision = JoinCollision.error,
)

/**
 * Joins [tables] by taxon and shared factors, as ADR-0024 defines. Two inputs join as left and right; more than two
 * are a left fold in list order.
 *
 * The key is the taxon plus every factor the two sides share by name, matched on level label. A missing factor value
 * matches nothing. Rows come out left first, each followed by its right matches, then (for a union) the unmatched
 * right rows; columns are the left's, then the right's new ones. Nothing is sorted. Cells a side cannot supply are
 * missing.
 *
 * @throws IllegalArgumentException for fewer than two tables, a key that repeats on both sides under
 *   [JoinMultiple.error], a shared non-key name under [JoinCollision.error], a coalesce of unequal values or of
 *   columns of different types, or a suffixed name that is already taken.
 * @see net.maizegenetics.phenotype.PhenotypeBuilder.mergeTwoPhenotypes
 */
suspend fun joinPhenotypes(
    tables: List<PhenotypeTable>,
    options: JoinOptions = JoinOptions(),
    progress: ProgressReporter = ProgressReporter.None,
): PhenotypeTable {
    require(tables.size >= 2) { "A join needs at least two phenotype tables, got ${tables.size}" }
    val inputs = resolveCollisions(tables, options.collision)
    var result = inputs[0]
    for (next in 1 until inputs.size) {
        result = joinPair(result, inputs[next], options)
        progress.report(next, inputs.size - 1, "Joined ${tables[next].name}")
    }
    return PhenotypeTable(result.attributes, "joined")
}

/**
 * Checks, over the original inputs, the names held by two or more tables that are not factors in all of them, and
 * applies the collision rule: fail, rename, or (when the holders have one type) leave for the pairwise joins to
 * coalesce. Doing it up front makes `suffix` independent of how the fold groups the inputs.
 */
private fun resolveCollisions(tables: List<PhenotypeTable>, rule: JoinCollision): List<PhenotypeTable> {
    val holders = LinkedHashMap<String, MutableList<Pair<Int, PhenotypeAttribute>>>()
    tables.forEachIndexed { position, table ->
        table.attributes.filter { it.type != AttributeType.taxa }.forEach {
            holders.getOrPut(it.name) { mutableListOf() }.add(position to it)
        }
    }
    val colliding = holders.filter { (_, held) -> held.size > 1 && held.any { it.second.type != AttributeType.factor } }
    if (colliding.isEmpty()) return tables
    when (rule) {
        JoinCollision.error -> throw IllegalArgumentException(
            "Columns ${colliding.keys.quoted()} are in more than one table and are not shared factors; " +
                "set collision to suffix or coalesce",
        )

        JoinCollision.coalesce -> colliding.forEach { (name, held) ->
            val types = held.map { it.second.type }.distinct()
            require(types.size == 1) {
                "Column '$name' is ${types.joinToString(" and ")} in different tables, which coalesce cannot merge; " +
                    "use collision suffix"
            }
        }

        JoinCollision.suffix -> {
            val taken = holders.keys.toMutableSet().also { names -> tables.forEach { names += it.taxaAttribute.name } }
            val renamed = HashMap<Pair<Int, String>, String>()
            colliding.forEach { (name, held) ->
                held.forEach { (position, _) ->
                    val newName = "$name.${position + 1}"
                    require(taken.add(newName)) { "Cannot rename '$name' to '$newName': that name is already in use" }
                    renamed[position to name] = newName
                }
            }
            return tables.mapIndexed { position, table ->
                PhenotypeTable(
                    table.attributes.map { attribute ->
                        renamed[position to attribute.name]?.let { attribute.renamed(it) } ?: attribute
                    },
                    table.name,
                )
            }
        }
    }
    return tables
}

private fun PhenotypeAttribute.renamed(name: String): PhenotypeAttribute = when (this) {
    is TaxaAttribute -> TaxaAttribute(name, names)
    is NumericAttribute -> NumericAttribute(name, type, toFloatArray())
    is FactorAttribute -> FactorAttribute.of(name, (0 until size).map { level(it) })
}

/** One join of [left] with [right]; shared non-key names reaching here are for `coalesce` only. */
private fun joinPair(left: PhenotypeTable, right: PhenotypeTable, options: JoinOptions): PhenotypeTable {
    val leftTaxa = left.taxaAttribute.name
    val rightTaxa = right.taxaAttribute.name
    require(leftTaxa == rightTaxa || (right.indexOf(leftTaxa) < 0 && left.indexOf(rightTaxa) < 0)) {
        "A taxa column is named '$leftTaxa' or '$rightTaxa' and is another kind of column in the other table"
    }
    val keyNames = left.attributes
        .filter {
            it is FactorAttribute &&
                right.indexOf(it.name).let { c -> c >= 0 && right.attributes[c] is FactorAttribute }
        }
        .map { it.name }
    val leftKeys = keysOf(left, keyNames.map { left.attribute(it) as FactorAttribute })
    val rightKeys = keysOf(right, keyNames.map { right.attribute(it) as FactorAttribute })

    val rightRows = LinkedHashMap<List<String>, MutableList<Int>>()
    rightKeys.forEachIndexed { row, key -> if (key != null) rightRows.getOrPut(key) { mutableListOf() }.add(row) }
    val leftCounts = HashMap<List<String>, Int>()
    leftKeys.forEach { if (it != null) leftCounts[it] = (leftCounts[it] ?: 0) + 1 }

    if (options.multiple == JoinMultiple.error) {
        val ambiguous = leftCounts.keys.filter { leftCounts.getValue(it) > 1 && (rightRows[it]?.size ?: 0) > 1 }
        require(ambiguous.isEmpty()) {
            val ordered = leftKeys.filterNotNull().distinct().filter { it in ambiguous }
            val shown = ordered.take(MAX_KEYS_SHOWN).joinToString("; ") { describe(it, keyNames) }
            val more = if (ordered.size > MAX_KEYS_SHOWN) "; and ${ordered.size - MAX_KEYS_SHOWN} more" else ""
            "Keys repeat in both tables, so the join is ambiguous: $shown$more; set multiple to cartesian to pair them all"
        }
    }

    val leftOf = IntBuffer()
    val rightOf = IntBuffer()
    val rightMatched = BooleanArray(right.numberOfObservations)
    for (row in 0 until left.numberOfObservations) {
        val matches = leftKeys[row]?.let { rightRows[it] }
        if (matches == null) {
            if (options.type == JoinType.union) {
                leftOf.add(row)
                rightOf.add(NONE)
            }
            continue
        }
        for (match in matches) {
            leftOf.add(row)
            rightOf.add(match)
            rightMatched[match] = true
        }
    }
    if (options.type == JoinType.union) {
        for (row in 0 until right.numberOfObservations) {
            if (!rightMatched[row]) {
                leftOf.add(NONE)
                rightOf.add(row)
            }
        }
    }
    val rows = leftOf.size
    val l = leftOf.toArray()
    val r = rightOf.toArray()

    val attributes = ArrayList<PhenotypeAttribute>()
    attributes += TaxaAttribute(leftTaxa, List(rows) { if (l[it] != NONE) left.taxon(l[it]) else right.taxon(r[it]) })
    for (attribute in left.attributes) {
        if (attribute.type == AttributeType.taxa) continue
        val other = right.indexOf(attribute.name).let { if (it < 0) null else right.attributes[it] }
        attributes += when {
            attribute is FactorAttribute && other is FactorAttribute -> FactorAttribute.of(
                attribute.name,
                List(rows) { if (l[it] != NONE) attribute.level(l[it]) else other.level(r[it]) },
            )

            attribute is FactorAttribute -> FactorAttribute.of(
                attribute.name,
                List(rows) { if (l[it] != NONE) attribute.level(l[it]) else null },
            )

            attribute is NumericAttribute && other != null -> coalesce(
                attribute,
                other as NumericAttribute,
                l,
                r,
                left,
                keyNames,
            )
            attribute is NumericAttribute -> NumericAttribute(
                attribute.name,
                attribute.type,
                FloatArray(rows) { if (l[it] != NONE) attribute.floatValue(l[it]) else Float.NaN },
            )

            else -> error("unreachable: ${attribute.name}")
        }
    }
    for (attribute in right.attributes) {
        if (attribute.type == AttributeType.taxa || left.indexOf(attribute.name) >= 0) continue
        attributes += when (attribute) {
            is FactorAttribute -> FactorAttribute.of(
                attribute.name,
                List(rows) { if (r[it] != NONE) attribute.level(r[it]) else null },
            )

            is NumericAttribute -> NumericAttribute(
                attribute.name,
                attribute.type,
                FloatArray(rows) { if (r[it] != NONE) attribute.floatValue(r[it]) else Float.NaN },
            )

            is TaxaAttribute -> error("unreachable: ${attribute.name}")
        }
    }
    return PhenotypeTable(attributes, "joined")
}

/** The column [from] holds in a pair of rows: the left's value, or the right's when the left's is missing. */
private fun coalesce(
    from: NumericAttribute,
    other: NumericAttribute,
    l: IntArray,
    r: IntArray,
    left: PhenotypeTable,
    keyNames: List<String>,
): NumericAttribute {
    val values = FloatArray(l.size) { row ->
        val a = if (l[row] != NONE) from.floatValue(l[row]) else Float.NaN
        val b = if (r[row] != NONE) other.floatValue(r[row]) else Float.NaN
        when {
            a.isNaN() -> b
            b.isNaN() -> a
            a == b -> a
            else -> {
                val key = left.taxon(l[row]) + keyNames.joinToString("") { name ->
                    ", $name=" + (left.attribute(name) as FactorAttribute).format(l[row])
                }
                throw IllegalArgumentException("Column '${from.name}' has $a and $b for $key; coalesce cannot choose")
            }
        }
    }
    return NumericAttribute(from.name, from.type, values)
}

/** Each row's key (taxon, then the level of each factor), or `null` when a factor value is missing. */
private fun keysOf(table: PhenotypeTable, factors: List<FactorAttribute>): List<List<String>?> =
    List(table.numberOfObservations) { row ->
        val levels = factors.map { it.level(row) }
        if (levels.any { it == null }) null else listOf(table.taxon(row)) + levels.filterNotNull()
    }

private fun describe(key: List<String>, factorNames: List<String>): String = if (factorNames.isEmpty()) {
    key[0]
} else {
    key[0] + " (" +
        factorNames.indices.joinToString(", ") { "${factorNames[it]}=${key[it + 1]}" } +
        ")"
}

private fun Collection<String>.quoted(): String = joinToString(", ") { "'$it'" }

private const val NONE = -1
private const val MAX_KEYS_SHOWN = 10

/** A growable list of ints, so a large join does not box every index. */
private class IntBuffer {
    private var data = IntArray(16)
    var size = 0
        private set

    fun add(value: Int) {
        if (size == data.size) data = data.copyOf(size * 2)
        data[size++] = value
    }

    fun toArray(): IntArray = data.copyOf(size)
}

/**
 * The spec for [joinPhenotypes]. Declared here; `TasselOperations` registers it (#166).
 *
 * The input is a [DataSet] of at least two phenotype tables and nothing else. The legacy plugin is TASSEL 5's
 * `IntersectionAlignmentPlugin`, which is the default `type`; the translator maps `UnionAlignmentPlugin` to this
 * operation with `type = union`. The three parameters are TASSEL 6's and have no legacy names.
 *
 * @see net.maizegenetics.analysis.data.IntersectionAlignmentPlugin
 * @see net.maizegenetics.analysis.data.UnionAlignmentPlugin
 */
val joinPhenotypesOperation = operation<DataSet, PhenotypeTable>(
    id = "join-phenotypes",
    legacyPlugin = "IntersectionAlignmentPlugin",
    help = "Joins two or more phenotype tables by taxon and shared factors, keeping matched rows or all rows",
) {
    val type by enum(JoinType.intersect, help = "intersect keeps matched rows; union keeps every row")
    val multiple by enum(
        JoinMultiple.error,
        help = "When a key repeats in both tables: error, or cartesian to pair all",
    )
    val collision by enum(
        JoinCollision.error,
        help = "A column both tables have: error, suffix (name.1, name.2), or coalesce",
    )

    run { input, progress ->
        val tables = input.phenotypes()
        require(tables.size >= 2 && tables.size == input.size) {
            "join-phenotypes needs at least two phenotype tables and nothing else; " +
                "got ${tables.size} phenotype tables among ${input.size} items"
        }
        joinPhenotypes(tables, JoinOptions(type(), multiple(), collision()), progress)
    }
}
