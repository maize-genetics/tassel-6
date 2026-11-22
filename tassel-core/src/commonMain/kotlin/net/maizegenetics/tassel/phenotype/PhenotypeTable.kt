package net.maizegenetics.tassel.table

/**
 * A TASSEL 5 `Phenotype`: one row per observation, and one named, typed attribute per column.
 *
 * This is not a [ResultTable]. TASSEL writes a `Phenotype` through `PhenotypeUtils.write`, which
 * puts a `<Phenotype>` line and a line of attribute types ahead of the column names, and it holds
 * numeric attributes as `Float`. [ResultTable] rejects a `Float` cell on purpose, so that a float
 * cannot slip past [DoubleFormat]; a `Phenotype` is the place floats legitimately go, and they are
 * written with [floatToString] rather than rounded.
 *
 * @see net.maizegenetics.phenotype.Phenotype
 * @see net.maizegenetics.phenotype.PhenotypeUtils.write
 */
class PhenotypeTable(
    /** The attributes, in column order. The first is normally the [AttributeType.taxa] one. */
    val attributes: List<Attribute>,
) {
    init {
        require(attributes.isNotEmpty()) { "A phenotype needs at least one attribute" }
        val sizes = attributes.map { it.values.size }.distinct()
        require(sizes.size == 1) {
            "Every attribute needs the same number of observations, found $sizes"
        }
    }

    /** Number of rows. */
    val numberOfObservations: Int get() = attributes[0].values.size

    /** Number of columns. */
    val numberOfAttributes: Int get() = attributes.size

    /**
     * Appends the phenotype as TASSEL 5's `PhenotypeUtils.write` writes it: a `<Phenotype>` line,
     * a line of attribute types, a line of attribute names, then one line per observation. Lines
     * end with `\n`, including the last. Cells follow `saveDelimitedTableReport`: `null` is
     * written as nothing, a `Double` through [DoubleFormat], and everything else as its string,
     * which for a `Float` is Java's `Float.toString` spelling.
     *
     * @see net.maizegenetics.util.TableReportUtils.saveDelimitedTableReport
     */
    fun writeDelimited(out: Appendable, delimiter: String = "\t") {
        out.append(HEADER).append('\n')
        appendLine(out, delimiter) { it.type.name }
        appendLine(out, delimiter) { it.name }
        for (observation in 0 until numberOfObservations) {
            attributes.forEachIndexed { column, attribute ->
                if (column != 0) out.append(delimiter)
                out.append(formatCell(attribute.values[observation]))
            }
            out.append('\n')
        }
    }

    /** [writeDelimited] into a string. */
    fun toDelimitedString(delimiter: String = "\t"): String = buildString { writeDelimited(this, delimiter) }

    override fun toString(): String =
        "PhenotypeTable($numberOfObservations observations, $numberOfAttributes attributes)"

    private inline fun appendLine(out: Appendable, delimiter: String, cell: (Attribute) -> String) {
        attributes.forEachIndexed { column, attribute ->
            if (column != 0) out.append(delimiter)
            out.append(cell(attribute))
        }
        out.append('\n')
    }

    /**
     * What an attribute means. The constant names are TASSEL 5's, because they are the words
     * written on the second line of the file and read back by its loader.
     *
     * @see net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE
     */
    @Suppress("ktlint:standard:enum-entry-name-case")
    enum class AttributeType { taxa, factor, covariate, data }

    /**
     * One column. A cell is a `String`, `Float`, `Double`, `Int`, `Long`, or `null`; the
     * constructor rejects anything else, as [ResultTable] does.
     */
    class Attribute(val name: String, val type: AttributeType, val values: List<Any?>) {
        init {
            values.forEach { cell ->
                val usable =
                    cell == null || cell is String || cell is Float || cell is Double || cell is Int || cell is Long
                require(usable) {
                    "Attribute '$name' holds a ${cell!!::class.simpleName} cell; " +
                        "use a String, Float, Double, Int, Long, or null"
                }
            }
        }

        /** A covariate column of [values], which is how TASSEL 5 stores principal components. */
        constructor(name: String, values: FloatArray) :
            this(name, AttributeType.covariate, values.toList())
    }

    companion object {
        /** The first line of the file, which is what TASSEL 5's loader recognizes the format by. */
        const val HEADER: String = "<Phenotype>"

        /** TASSEL 5's name for the taxa attribute. */
        const val TAXA: String = "Taxa"

        /** A taxa attribute named [TAXA] holding [names]. */
        fun taxa(names: List<String>): Attribute = Attribute(TAXA, AttributeType.taxa, names)

        /** One cell as TASSEL 5 writes it; see [writeDelimited]. */
        fun formatCell(cell: Any?): String = when (cell) {
            null -> ""
            is Double -> DoubleFormat.format(cell)
            is Float -> floatToString(cell)
            else -> cell.toString()
        }
    }
}
