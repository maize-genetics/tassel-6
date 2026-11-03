package net.maizegenetics.tassel.taxa

/**
 * One sample (line, individual) in a genotype table, identified by name.
 *
 * @see net.maizegenetics.taxa.Taxon
 */
data class Taxon(val name: String) : Comparable<Taxon> {
    override fun compareTo(other: Taxon): Int = name.compareTo(other.name)

    override fun toString(): String = name
}

/**
 * The ordered taxa of a genotype table. Names must be unique; [indexOf] by name is constant time.
 *
 * @see net.maizegenetics.taxa.TaxaList
 */
class TaxaList(taxa: List<Taxon>) : AbstractList<Taxon>() {
    private val taxa: List<Taxon> = taxa.toList()
    private val indexByName: Map<String, Int> = HashMap<String, Int>(taxa.size).also { index ->
        this.taxa.forEachIndexed { i, taxon ->
            require(index.put(taxon.name, i) == null) { "Taxon [${taxon.name}] appears more than once" }
        }
    }

    override val size: Int get() = taxa.size

    override fun get(index: Int): Taxon = taxa[index]

    /** Index of the taxon called [name], or -1. */
    fun indexOf(name: String): Int = indexByName[name] ?: -1

    override fun indexOf(element: Taxon): Int = indexOf(element.name)

    override fun lastIndexOf(element: Taxon): Int = indexOf(element.name)

    override fun contains(element: Taxon): Boolean = element.name in indexByName

    /** Taxa names, in order. */
    val names: List<String> get() = taxa.map { it.name }

    companion object {
        /** A list of taxa with the given names, in order. */
        fun of(names: Iterable<String>): TaxaList = TaxaList(names.map(::Taxon))

        /** A list of taxa with the given names, in order. */
        fun of(vararg names: String): TaxaList = of(names.asIterable())
    }
}
