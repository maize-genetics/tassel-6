package net.maizegenetics.tassel.distance

import net.maizegenetics.tassel.table.floatToString
import net.maizegenetics.tassel.taxa.TaxaList

/**
 * A symmetric matrix of distances, kinship coefficients, or any other pairwise measure over [taxa],
 * with free-text annotations describing how it was made.
 *
 * Values are [Float], as TASSEL 5 stores them. That is the precision ceiling of every TASSEL
 * distance matrix and of the text files TASSEL writes, so a port that rounded elsewhere would not
 * reproduce them. It is also why the kinship golden is compared at 8e-6 rather than something
 * tighter (ADR-0009).
 *
 * @see net.maizegenetics.taxa.distance.DistanceMatrix
 * @see net.maizegenetics.taxa.distance.DistanceMatrixBuilder
 */
class DistanceMatrix(
    /** The taxa, in row and column order. */
    val taxa: TaxaList,
    /**
     * The lower triangle including the diagonal, row by row: the value for `row` and `column`, with
     * `column <= row`, is at `row * (row + 1) / 2 + column`.
     */
    lowerTriangle: FloatArray,
    annotations: Map<String, String> = emptyMap(),
) {
    private val values = lowerTriangle.copyOf()

    /** Annotations by key, in the order TASSEL 5 writes them, which is sorted by key. */
    val annotations: Map<String, String> =
        annotations.entries.sortedBy { it.key }.associate { it.key to it.value }

    init {
        require(values.size == triangleSize(taxa.size)) {
            "A ${taxa.size}-taxon matrix needs ${triangleSize(taxa.size)} values, not ${values.size}"
        }
    }

    /** Number of taxa, which is both the row and the column count. */
    val numberOfTaxa: Int get() = taxa.size

    /** The value for [row] and [column], in either order. */
    fun distance(row: Int, column: Int): Float {
        if (row !in 0 until numberOfTaxa) throw IndexOutOfBoundsException("Row $row of $numberOfTaxa")
        if (column !in 0 until numberOfTaxa) throw IndexOutOfBoundsException("Column $column of $numberOfTaxa")
        val lower = maxOf(row, column)
        val upper = minOf(row, column)
        return values[lower * (lower + 1) / 2 + upper]
    }

    /**
     * Writes TASSEL 5's square-matrix text: one `##key=value` line per annotation, the taxa count,
     * then one line per taxon holding its name and its whole row. Lines end with `\n`, including the
     * last, and values are spelled as Java's `Float.toString` spells them.
     *
     * @see net.maizegenetics.taxa.distance.WriteDistanceMatrix.saveDelimitedDistanceMatrix
     */
    fun writeDelimited(out: Appendable, delimiter: String = "\t") {
        annotations.forEach { (key, value) ->
            out.append("##").append(key).append('=').append(value).append('\n')
        }
        out.append(numberOfTaxa.toString()).append('\n')
        for (row in 0 until numberOfTaxa) {
            out.append(taxa[row].name)
            for (column in 0 until numberOfTaxa) {
                out.append(delimiter).append(floatToString(distance(row, column)))
            }
            out.append('\n')
        }
    }

    /** [writeDelimited] into a string. */
    fun toDelimitedString(delimiter: String = "\t"): String = buildString { writeDelimited(this, delimiter) }

    /**
     * The matrix over [names], in the order given, as a new matrix with the same annotations.
     *
     * A name asked for twice is kept at its first position and listed in [DistanceMatrixSubset.repeated].
     * A name this matrix lacks is an error, unless [missing] is [MissingTaxa.Drop]; then it is left
     * out and listed in [DistanceMatrixSubset.missing]. Names match exactly, case included.
     *
     * @throws IllegalArgumentException for a missing name under [MissingTaxa.Error]; the message
     *   gives the count and the first few names.
     */
    fun subset(names: List<String>, missing: MissingTaxa = MissingTaxa.Error): DistanceMatrixSubset {
        val chosen = ArrayList<Int>(names.size)
        val kept = HashSet<String>(names.size)
        val absent = ArrayList<String>()
        val repeated = ArrayList<String>()
        for (name in names) {
            val index = taxa.indexOf(name)
            when {
                index < 0 -> if (name !in absent) absent += name
                kept.add(name) -> chosen += index
                else -> if (name !in repeated) repeated += name
            }
        }
        if (absent.isNotEmpty() && missing == MissingTaxa.Error) {
            val shown = absent.take(MAX_NAMED_MISSING).joinToString { "'$it'" }
            val more = if (absent.size > MAX_NAMED_MISSING) ", ..." else ""
            throw IllegalArgumentException("${absent.size} taxa are not in the matrix: $shown$more")
        }
        val size = chosen.size
        val subsetValues = FloatArray(triangleSize(size))
        var next = 0
        for (row in 0 until size) {
            for (column in 0..row) {
                subsetValues[next++] = distance(chosen[row], chosen[column])
            }
        }
        val subsetTaxa = TaxaList(chosen.map { taxa[it] })
        return DistanceMatrixSubset(DistanceMatrix(subsetTaxa, subsetValues, annotations), absent, repeated)
    }

    companion object {
        private const val MAX_NAMED_MISSING = 5

        /** Annotation key naming the kind of matrix, for example `Centered_IBS`. */
        const val MATRIX_TYPE: String = "Matrix_Type"

        /** Annotation key for Centered IBS's sum of allele-frequency terms, which scales the matrix. */
        const val CENTERED_IBS_SUMPK: String = "Centered_IBS.SumPk"

        /** Number of values in the lower triangle, including the diagonal, for [numberOfTaxa] taxa. */
        fun triangleSize(numberOfTaxa: Int): Int = numberOfTaxa * (numberOfTaxa + 1) / 2
    }
}

/** What [DistanceMatrix.subset] does with a requested name the matrix lacks. */
enum class MissingTaxa {
    /** Fail, naming the taxa. */
    Error,

    /** Leave them out and report them. */
    Drop,
}

/**
 * A matrix restricted to some taxa. [missing] holds the requested names the source lacked, which
 * are empty unless [MissingTaxa.Drop] was passed, and [repeated] the names requested more than once.
 */
class DistanceMatrixSubset(val matrix: DistanceMatrix, val missing: List<String>, val repeated: List<String>)
