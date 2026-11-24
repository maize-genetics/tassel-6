package net.maizegenetics.tassel.dataset

import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.phenotype.PhenotypeTable

/**
 * The data an operation with several inputs works on: an ordered list of core values, such as the phenotype tables
 * selected for a join. ADR-0024 records why this is a plain value and not a change to `Operation` or `Engine`.
 *
 * Items keep their order, which is the order a front end sent them in and the order joins read them as left and
 * right. A `DataSet` inside [items] is flattened into the one being built. The accessors do not check counts; each
 * operation checks its own arity and says what it needs.
 *
 * @see net.maizegenetics.plugindef.DataSet
 */
class DataSet(items: List<Any>) {
    /** The items, in order, with nested data sets flattened. */
    val items: List<Any> = items.flatMap { if (it is DataSet) it.items else listOf(it) }

    /** The phenotype tables among the items, in order. */
    fun phenotypes(): List<PhenotypeTable> = items.filterIsInstance<PhenotypeTable>()

    /** The genotype tables among the items, in order. */
    fun genotypes(): List<GenotypeTable> = items.filterIsInstance<GenotypeTable>()

    /** Number of items. */
    val size: Int get() = items.size

    override fun equals(other: Any?): Boolean = other is DataSet && other.items == items

    override fun hashCode(): Int = items.hashCode()

    override fun toString(): String = "DataSet(${items.joinToString { it::class.simpleName ?: "?" }})"

    companion object {
        /** A data set of [items], in the order given. */
        fun of(vararg items: Any): DataSet = DataSet(items.toList())
    }
}
