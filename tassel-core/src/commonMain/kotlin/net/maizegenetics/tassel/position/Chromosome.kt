package net.maizegenetics.tassel.position

/**
 * A chromosome (or contig), named as TASSEL 5 names it.
 *
 * [of] normalizes the name: it trims it, upper-cases it, drops a leading `CHROMOSOME` or `CHR`, and keeps only the text
 * before the first space, so `chr1`, `Chromosome1`, and `1 extra` are all chromosome `1`.
 *
 * Ordering follows TASSEL 5: two chromosomes with integer names compare as integers; otherwise names compare as
 * text, with a leading number zero-padded to five digits so that `2A` sorts before `10A`. Like TASSEL 5's, this
 * ordering is not transitive across a mix of numeric and non-numeric names.
 *
 * @see net.maizegenetics.dna.map.Chromosome
 */
class Chromosome private constructor(
    /** The normalized name, such as `1` or `SCAFFOLD_7`. */
    val name: String,
) : Comparable<Chromosome> {
    /** The name as an integer, or null if it is not one. */
    val number: Int? = name.toIntOrNull()?.takeIf { it != Int.MAX_VALUE }

    private val compareKey: String = run {
        val digits = name.takeWhile { it in '0'..'9' }.length
        if (digits > 0) "0".repeat(maxOf(0, NUMBER_WIDTH - digits)) + name else name
    }

    override fun compareTo(other: Chromosome): Int {
        if (this === other) return 0
        val number = number
        val otherNumber = other.number
        if (number != null && otherNumber != null) return number.compareTo(otherNumber)
        return compareKey.compareTo(other.compareKey)
    }

    override fun equals(other: Any?): Boolean = other is Chromosome && name == other.name

    override fun hashCode(): Int = name.hashCode()

    override fun toString(): String = name

    companion object {
        private const val NUMBER_WIDTH = 5

        /**
         * The chromosome called [name], normalized as described on [Chromosome].
         *
         * @throws IllegalArgumentException if [name] is empty
         */
        fun of(name: String): Chromosome {
            require(name.isNotEmpty()) { "Chromosome name can't be empty" }
            return Chromosome(normalize(name))
        }

        /** The chromosome with an integer name. */
        fun of(number: Int): Chromosome = of(number.toString())

        private fun normalize(name: String): String {
            // Java's String.trim: strips every character up to and including space, not only whitespace.
            var result = name.trim { it <= ' ' }.uppercase()
            result = result.removePrefix("CHROMOSOME")
            result = result.removePrefix("CHR")
            val space = result.indexOf(' ')
            return if (space > 0) result.substring(0, space) else result
        }
    }
}
