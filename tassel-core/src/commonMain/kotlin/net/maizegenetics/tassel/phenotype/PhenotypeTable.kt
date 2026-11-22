package net.maizegenetics.tassel.phenotype

import net.maizegenetics.tassel.taxa.TaxaList

/**
 * A phenotype: one row per observation and one typed, named [PhenotypeAttribute] per column. A
 * taxon may have several observations (replicates, environments), told apart by factor attributes.
 *
 * A table is immutable. It holds exactly one [TaxaAttribute], and attribute names are unique.
 * [selectRows] and [selectAttributes] return new tables. ADR-0023 records the model and how it
 * differs from TASSEL 5's.
 *
 * This is not a [net.maizegenetics.tassel.table.ResultTable]: numeric attributes hold `Float`s, as
 * TASSEL 5's do, and [writeDelimited] writes them with Java's `Float.toString` rather than
 * rounding them.
 *
 * @see net.maizegenetics.phenotype.Phenotype
 * @see net.maizegenetics.phenotype.CorePhenotype
 */
class PhenotypeTable(
    attributes: List<PhenotypeAttribute>,
    /** A display name, such as the file it was read from. Not written to the file. */
    val name: String = DEFAULT_NAME,
) {
    /** The attributes, in column order. */
    val attributes: List<PhenotypeAttribute> = attributes.toList()

    private val indexByName: Map<String, Int>

    init {
        require(this.attributes.isNotEmpty()) { "A phenotype needs at least one attribute" }
        val sizes = this.attributes.map { it.size }.distinct()
        require(sizes.size == 1) { "Every attribute needs the same number of observations, found $sizes" }
        val taxaColumns = this.attributes.count { it.type == AttributeType.taxa }
        require(taxaColumns == 1) { "A phenotype needs exactly one taxa attribute, found $taxaColumns" }
        indexByName = HashMap<String, Int>(this.attributes.size * 2).also { index ->
            this.attributes.forEachIndexed { column, attribute ->
                require(index.put(attribute.name, column) == null) {
                    "Attribute name '${attribute.name}' appears more than once"
                }
            }
        }
    }

    /** Number of rows. */
    val numberOfObservations: Int get() = attributes[0].size

    /** Number of columns, the taxa attribute included. */
    val numberOfAttributes: Int get() = attributes.size

    /** The one taxa attribute. */
    val taxaAttribute: TaxaAttribute = this.attributes.first { it.type == AttributeType.taxa } as TaxaAttribute

    /** The taxon of [observation]. */
    fun taxon(observation: Int): String = taxaAttribute.taxon(observation)

    /** The distinct taxa, in the order they first appear. */
    val taxa: TaxaList by lazy { TaxaList.of(taxaAttribute.names.distinct()) }

    /** Whether some taxon has more than one observation. */
    val hasReplicatedTaxa: Boolean get() = taxa.size < numberOfObservations

    /** The column index of the attribute called [name], or -1. */
    fun indexOf(name: String): Int = indexByName[name] ?: -1

    /** The attribute called [name]. */
    fun attribute(name: String): PhenotypeAttribute {
        val column = indexOf(name)
        require(column >= 0) { "No attribute named '$name'; the attributes are ${attributes.map { it.name }}" }
        return attributes[column]
    }

    /** The attributes of [type], in column order. */
    fun attributes(type: AttributeType): List<PhenotypeAttribute> = attributes.filter { it.type == type }

    /** A table of the observations at [observations], in that order; an index may repeat. */
    fun selectRows(observations: IntArray): PhenotypeTable {
        observations.forEach { observation ->
            require(observation in 0 until numberOfObservations) {
                "Observation $observation is outside 0 until $numberOfObservations"
            }
        }
        return PhenotypeTable(attributes.map { it.select(observations) }, name)
    }

    /**
     * A table of the attributes called [names], in that order. The taxa attribute is kept even when
     * it is not named, as the first column.
     */
    fun selectAttributes(names: List<String>): PhenotypeTable {
        require(names.size == names.distinct().size) { "Attribute names repeat in $names" }
        val chosen = names.map { attribute(it) }
        val withTaxa = if (taxaAttribute in chosen) chosen else listOf(taxaAttribute) + chosen
        return PhenotypeTable(withTaxa, name)
    }

    /**
     * Appends the phenotype as TASSEL 5's `PhenotypeUtils.write` writes it: a `<Phenotype>` line,
     * a line of attribute types, a line of attribute names, then one line per observation, each
     * cell as [PhenotypeAttribute.format] spells it. Lines end with `\n`, including the last.
     *
     * The format splits on whitespace, and reads `NA`, `NaN`, `.`, and `?` as missing, so a name,
     * taxon, or factor level that would not read back as itself is rejected.
     *
     * @see net.maizegenetics.phenotype.PhenotypeUtils.write
     * @see net.maizegenetics.util.TableReportUtils.saveDelimitedTableReport
     */
    fun writeDelimited(out: Appendable, delimiter: String = "\t") {
        checkWritable()
        out.append(HEADER).append('\n')
        appendLine(out, delimiter) { it.type.name }
        appendLine(out, delimiter) { it.name }
        for (observation in 0 until numberOfObservations) {
            appendLine(out, delimiter) { it.format(observation) }
        }
    }

    /** [writeDelimited] into a string. */
    fun toDelimitedString(delimiter: String = "\t"): String = buildString { writeDelimited(this, delimiter) }

    override fun equals(other: Any?): Boolean =
        other is PhenotypeTable && other.name == name && other.attributes == attributes

    override fun hashCode(): Int = 31 * name.hashCode() + attributes.hashCode()

    override fun toString(): String =
        "PhenotypeTable($name, $numberOfObservations observations, $numberOfAttributes attributes)"

    private inline fun appendLine(out: Appendable, delimiter: String, cell: (PhenotypeAttribute) -> String) {
        attributes.forEachIndexed { column, attribute ->
            if (column != 0) out.append(delimiter)
            out.append(cell(attribute))
        }
        out.append('\n')
    }

    private fun checkWritable() {
        attributes.forEach { attribute ->
            requireWritable(attribute.name, "Attribute name")
            when (attribute) {
                is TaxaAttribute -> attribute.names.forEach { requireWritable(it, "Taxon") }
                is FactorAttribute -> attribute.levels.forEach { level ->
                    requireWritable(level, "Level of factor '${attribute.name}'")
                    require(!PhenotypeFiles.isDefaultMissingToken(level)) {
                        "Level '$level' of factor '${attribute.name}' would be read back as missing"
                    }
                }
                is NumericAttribute -> Unit
            }
        }
    }

    private fun requireWritable(text: String, what: String) {
        val clean = text.none { it.isFieldSeparator() }
        require(clean) { "$what '$text' contains whitespace, which the file format splits on" }
    }

    companion object {
        /** The first line of the file, which is what TASSEL 5's loader recognizes the format by. */
        const val HEADER: String = "<Phenotype>"

        /** TASSEL 5's name for the taxa attribute, which TASSEL 6 uses when nothing else names it. */
        const val TAXA: String = "Taxa"

        /** The name of a table nobody named, as in TASSEL 5. */
        const val DEFAULT_NAME: String = "Phenotype"

        /** A taxa attribute named [TAXA] holding [names]. */
        fun taxa(names: List<String>): TaxaAttribute = TaxaAttribute(TAXA, names)
    }
}

/** The characters Java's `\s` matches, which is how TASSEL 5 splits phenotype lines. */
internal fun Char.isFieldSeparator(): Boolean = this == ' ' || this in '\t'..'\r'
