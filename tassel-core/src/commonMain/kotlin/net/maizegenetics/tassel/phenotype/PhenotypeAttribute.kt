package net.maizegenetics.tassel.phenotype

import net.maizegenetics.tassel.table.floatToString

/**
 * What an attribute means. The constant names are TASSEL 5's, because they are the words written
 * on the second line of a `<Phenotype>` file and read back by its loader.
 *
 * @see net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE
 */
@Suppress("ktlint:standard:enum-entry-name-case")
enum class AttributeType {
    /** The observation's taxon. A phenotype has exactly one. */
    taxa,

    /** A categorical attribute, such as an environment or a replicate. */
    factor,

    /** A numeric attribute used as an explanatory variable, such as a population-structure Q value. */
    covariate,

    /** A numeric trait. */
    data,
}

/**
 * One column of a [PhenotypeTable]: a name, an [AttributeType], and one value per observation.
 *
 * Attributes are immutable. Each kind stores its values in the form TASSEL 5 does, with one
 * missing-value representation (ADR-0023):
 * - [TaxaAttribute]: taxon names; never missing.
 * - [NumericAttribute]: `Float`s; missing is NaN.
 * - [FactorAttribute]: level codes into a sorted level list; missing is code -1.
 */
sealed class PhenotypeAttribute(
    /** The column name. Unique within a table. */
    val name: String,
) {
    init {
        require(name.isNotEmpty()) { "An attribute needs a name" }
    }

    /** What the attribute means. */
    abstract val type: AttributeType

    /** Number of observations. */
    abstract val size: Int

    /** Whether the value of [observation] is missing. */
    abstract fun isMissing(observation: Int): Boolean

    /** Number of missing values. */
    val missingCount: Int by lazy { (0 until size).count { isMissing(it) } }

    /** The value of [observation] as a `<Phenotype>` file writes it. */
    abstract fun format(observation: Int): String

    /** A new attribute holding the values of [observations], in that order; indices may repeat. */
    abstract fun select(observations: IntArray): PhenotypeAttribute

    override fun toString(): String = "${this::class.simpleName}($name, $size observations)"
}

/**
 * The taxon of each observation. A taxon may repeat, which is how replicates and environments are
 * held: one observation per row. Unlike TASSEL 5, which always names this column `Taxa`, the name
 * is whatever the file or caller gives it.
 *
 * @see net.maizegenetics.phenotype.TaxaAttribute
 */
class TaxaAttribute(name: String, taxa: List<String>) : PhenotypeAttribute(name) {
    private val taxa: List<String> = taxa.toList()

    init {
        this.taxa.forEachIndexed { observation, taxon ->
            require(taxon.isNotEmpty()) { "Taxa attribute '$name' has an empty taxon at observation $observation" }
        }
    }

    override val type: AttributeType get() = AttributeType.taxa
    override val size: Int get() = taxa.size

    /** The taxon name of [observation]. */
    fun taxon(observation: Int): String = taxa[observation]

    /** Every observation's taxon, in observation order. */
    val names: List<String> get() = taxa

    override fun isMissing(observation: Int): Boolean = false

    override fun format(observation: Int): String = taxa[observation]

    override fun select(observations: IntArray): TaxaAttribute = TaxaAttribute(name, observations.map { taxa[it] })

    override fun equals(other: Any?): Boolean = other is TaxaAttribute && other.name == name && other.taxa == taxa

    override fun hashCode(): Int = 31 * name.hashCode() + taxa.hashCode()
}

/**
 * A numeric attribute, [AttributeType.data] or [AttributeType.covariate], stored as `Float` the
 * way TASSEL 5 stores it. Analyses read values as `Double` through [value].
 *
 * NaN is the only missing marker. Infinite values are stored as NaN too, so there is a single
 * representation of "no usable value"; TASSEL 5 can hold an `Infinity` that is not flagged missing.
 *
 * @see net.maizegenetics.phenotype.NumericAttribute
 */
class NumericAttribute(name: String, override val type: AttributeType, values: FloatArray) : PhenotypeAttribute(name) {
    private val values: FloatArray =
        FloatArray(values.size) { values[it].let { v -> if (v.isFinite()) v else Float.NaN } }

    init {
        require(type == AttributeType.data || type == AttributeType.covariate) {
            "Numeric attribute '$name' must be data or covariate, not $type"
        }
    }

    /** A covariate attribute of [values], which is how TASSEL 5 stores principal components. */
    constructor(name: String, values: FloatArray) : this(name, AttributeType.covariate, values)

    override val size: Int get() = values.size

    /** The value of [observation], NaN when missing. */
    fun value(observation: Int): Double = values[observation].toDouble()

    /** The stored `Float` of [observation], NaN when missing. */
    fun floatValue(observation: Int): Float = values[observation]

    /** A copy of the stored values. */
    fun toFloatArray(): FloatArray = values.copyOf()

    override fun isMissing(observation: Int): Boolean = values[observation].isNaN()

    /** Java's `Float.toString` of the value, which writes a missing value as `NaN`. */
    override fun format(observation: Int): String = floatToString(values[observation])

    override fun select(observations: IntArray): NumericAttribute =
        NumericAttribute(name, type, FloatArray(observations.size) { values[observations[it]] })

    override fun equals(other: Any?): Boolean =
        other is NumericAttribute && other.name == name && other.type == type && other.values.contentEquals(values)

    override fun hashCode(): Int = (31 * name.hashCode() + type.hashCode()) * 31 + values.contentHashCode()
}

/**
 * A categorical attribute: one level per observation, held as a code into [levels].
 *
 * Levels are unique and sorted by `String.compareTo`, as TASSEL 5 sorts them, so `"10"` comes
 * before `"2"`. A missing value has code -1 and is written as [MISSING]. TASSEL 5 has no missing
 * factor value: there, `?` is an ordinary level.
 *
 * @see net.maizegenetics.phenotype.CategoricalAttribute
 */
class FactorAttribute(name: String, codes: IntArray, levels: List<String>) : PhenotypeAttribute(name) {
    private val codes: IntArray = codes.copyOf()

    /** The levels, sorted and unique. A level need not be used by any observation. */
    val levels: List<String> = levels.toList()

    init {
        this.levels.forEachIndexed { index, level ->
            require(level.isNotEmpty()) { "Factor '$name' has an empty level" }
            require(level != MISSING) { "Factor '$name' cannot have the level '$MISSING', which means missing" }
            require(index == 0 || this.levels[index - 1] < level) {
                "Factor '$name' levels must be unique and sorted; '${this.levels[index - 1]}' comes before '$level'"
            }
        }
        this.codes.forEachIndexed { observation, code ->
            require(code >= -1 && code < this.levels.size) {
                "Factor '$name' has code $code at observation $observation, outside -1 until ${this.levels.size}"
            }
        }
    }

    override val type: AttributeType get() = AttributeType.factor
    override val size: Int get() = codes.size

    /** The level code of [observation]: an index into [levels], or -1 when missing. */
    fun code(observation: Int): Int = codes[observation]

    /** The level of [observation], or `null` when missing. */
    fun level(observation: Int): String? = codes[observation].let { if (it < 0) null else levels[it] }

    override fun isMissing(observation: Int): Boolean = codes[observation] < 0

    /** The level, or [MISSING]. */
    override fun format(observation: Int): String = level(observation) ?: MISSING

    /** The selected observations, with only the levels they use. */
    override fun select(observations: IntArray): FactorAttribute = of(name, observations.map { level(it) })

    override fun equals(other: Any?): Boolean = other is FactorAttribute &&
        other.name == name &&
        other.levels == levels &&
        other.codes.contentEquals(codes)

    override fun hashCode(): Int = (31 * name.hashCode() + levels.hashCode()) * 31 + codes.contentHashCode()

    companion object {
        /** How a missing factor value is written: TASSEL 5's spelling, so TASSEL 5 can read the file. */
        const val MISSING: String = "?"

        /** A factor whose levels are the distinct non-null [labels], sorted; `null` is missing. */
        fun of(name: String, labels: List<String?>): FactorAttribute {
            val levels = labels.filterNotNull().distinct().sorted()
            val index = HashMap<String, Int>(levels.size * 2)
            levels.forEachIndexed { code, level -> index[level] = code }
            return FactorAttribute(name, IntArray(labels.size) { labels[it]?.let(index::getValue) ?: -1 }, levels)
        }
    }
}
