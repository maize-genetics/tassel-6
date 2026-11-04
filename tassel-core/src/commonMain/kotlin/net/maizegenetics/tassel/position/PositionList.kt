package net.maizegenetics.tassel.position

/**
 * The ordered sites of a genotype table. Sites must be sorted (see [Position.compareTo]), and each chromosome's sites
 * must be contiguous. Several sites may share a physical position.
 *
 * @see net.maizegenetics.dna.map.PositionList
 */
class PositionList(positions: List<Position>) : AbstractList<Position>() {
    private val positions: List<Position> = positions.toList()
    private val physical: IntArray = IntArray(this.positions.size) { this.positions[it].position }
    private val siteRanges: Map<Chromosome, IntRange>

    /** Chromosomes in site order. */
    val chromosomes: List<Chromosome>

    init {
        val ranges = LinkedHashMap<Chromosome, IntRange>()
        var start = 0
        for (site in 1..this.positions.size) {
            if (site < this.positions.size) {
                val previous = this.positions[site - 1]
                val current = this.positions[site]
                require(previous <= current) { "Positions are not in order: $previous before $current" }
                if (current.chromosome == previous.chromosome) continue
            }
            val chromosome = this.positions[start].chromosome
            require(ranges.put(chromosome, start until site) == null) {
                "Sites of chromosome $chromosome are not contiguous"
            }
            start = site
        }
        siteRanges = ranges
        chromosomes = ranges.keys.toList()
    }

    override val size: Int get() = positions.size

    override fun get(index: Int): Position = positions[index]

    /** Sites on [chromosome], or an empty range if it has none. */
    fun siteRange(chromosome: Chromosome): IntRange = siteRanges[chromosome] ?: IntRange.EMPTY

    /**
     * The first site at [position] on [chromosome], with TASSEL 5's return values: if no site is at [position], the
     * result is `-(insertion site) - 1`, where the insertion site is the index the position would have; if
     * [chromosome] has no sites, the result is [Int.MIN_VALUE].
     *
     * @see net.maizegenetics.dna.map.PositionList.siteOfPhysicalPosition
     */
    fun siteOfPhysicalPosition(chromosome: Chromosome, position: Int): Int {
        val range = siteRanges[chromosome] ?: return Int.MIN_VALUE
        var low = range.first
        var high = range.last + 1
        while (low < high) {
            val middle = (low + high) ushr 1
            if (physical[middle] < position) low = middle + 1 else high = middle
        }
        return if (low <= range.last && physical[low] == position) low else -low - 1
    }
}
