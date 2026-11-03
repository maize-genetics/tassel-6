package net.maizegenetics.tassel.position

/**
 * Strand of a site. Declared in TASSEL 5's byte order (unknown, minus, plus), which is the order positions sort in.
 *
 * @see net.maizegenetics.dna.map.Position
 */
enum class Strand(
    /** The symbol TASSEL 5 writes: `+`, `-`, or `N`. */
    val symbol: String,
) {
    UNKNOWN("N"),
    MINUS("-"),
    PLUS("+"),
    ;

    companion object {
        /** The strand written as [symbol] (`+`, `-`, or `N`). */
        fun of(symbol: String): Strand =
            entries.firstOrNull { it.symbol == symbol } ?: throw IllegalArgumentException("Unknown strand: $symbol")
    }
}

/**
 * A site: a physical position on a chromosome, with an optional SNP name.
 *
 * A site without a name is called `S<chromosome>_<position>`, as in TASSEL 5; passing that default as [snpId] is the
 * same as passing null. Positions sort by chromosome, then position, then strand, then name, and two positions are
 * equal when all four match.
 *
 * @see net.maizegenetics.dna.map.Position
 * @see net.maizegenetics.dna.map.GeneralPosition
 */
class Position(
    val chromosome: Chromosome,
    /** Physical position on [chromosome]. */
    val position: Int,
    snpId: String? = null,
    val strand: Strand = Strand.PLUS,
    /** Variants listed for the site in the input file, such as `["A", "G"]` from a HapMap `A/G`. */
    val knownVariants: List<String> = emptyList(),
) : Comparable<Position> {
    /** The SNP name given for this site, or null if it has only the default name. */
    val snpId: String? = snpId?.takeUnless { it == defaultName(chromosome, position) }

    /** The site name: [snpId], or `S<chromosome>_<position>` if there is none. */
    val name: String get() = snpId ?: defaultName(chromosome, position)

    override fun compareTo(other: Position): Int {
        chromosome.compareTo(other.chromosome).let { if (it != 0) return it }
        position.compareTo(other.position).let { if (it != 0) return it }
        strand.compareTo(other.strand).let { if (it != 0) return it }
        if (snpId == null && other.snpId == null) return 0
        return name.compareTo(other.name)
    }

    override fun equals(other: Any?): Boolean = other is Position &&
        chromosome == other.chromosome &&
        position == other.position &&
        strand == other.strand &&
        snpId == other.snpId

    override fun hashCode(): Int = ((chromosome.hashCode() * 31 + position) * 31 + strand.hashCode()) * 31 +
        snpId.hashCode()

    override fun toString(): String = "$name ($chromosome:$position)"

    private companion object {
        fun defaultName(chromosome: Chromosome, position: Int): String = "S${chromosome.name}_$position"
    }
}
