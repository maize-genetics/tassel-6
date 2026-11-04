package net.maizegenetics.tassel.filter

import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.operation
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.taxa.TaxaList

/**
 * Which taxa a [filterTaxa] call keeps.
 *
 * Defaults are TASSEL 5's, so `TaxaFilter()` keeps the table unchanged. [taxaNames] is applied first, then the
 * statistics thresholds, which are inclusive at both ends.
 *
 * @property taxaNames taxa to keep or to drop, or `null` to start from every taxon. Names the table does not have
 *   are ignored.
 * @property includeTaxa keeps [taxaNames] when true, drops them when false. Has no effect without [taxaNames].
 * @property minNotMissing minimum proportion of sites with a called genotype.
 * @property minHeterozygous minimum proportion of heterozygous sites.
 * @property maxHeterozygous maximum proportion of heterozygous sites.
 * @see net.maizegenetics.dna.snp.FilterTaxa
 */
data class TaxaFilter(
    val taxaNames: List<String>? = null,
    val includeTaxa: Boolean = true,
    val minNotMissing: Double = 0.0,
    val minHeterozygous: Double = 0.0,
    val maxHeterozygous: Double = 1.0,
) {
    init {
        require(minNotMissing in 0.0..1.0) { "minNotMissing must be between 0.0 and 1.0: $minNotMissing" }
        require(minHeterozygous in 0.0..1.0) { "minHeterozygous must be between 0.0 and 1.0: $minHeterozygous" }
        require(maxHeterozygous in 0.0..1.0) { "maxHeterozygous must be between 0.0 and 1.0: $maxHeterozygous" }
        require(taxaNames == null || taxaNames.size == taxaNames.toSet().size) {
            "taxaNames repeats a name: $taxaNames"
        }
    }

    /** True when a taxon's statistics decide whether it is kept, which is TASSEL 5's test for computing them. */
    internal val usesStats: Boolean
        get() = minNotMissing != 0.0 || minHeterozygous != 0.0 || maxHeterozygous != 1.0
}

/**
 * Keeps the taxa of [table] that [filter] selects, with all sites.
 *
 * Taxa named by [TaxaFilter.taxaNames] and kept come out in the order they are named, as TASSEL 5 does; taxa selected
 * any other way stay in table order. Returns [table] itself when [filter] would keep every taxon in place.
 *
 * @see net.maizegenetics.analysis.filter.FilterByTaxa
 * @see net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin
 */
fun filterTaxa(
    table: GenotypeTable,
    filter: TaxaFilter,
    progress: ProgressReporter = ProgressReporter.None,
): GenotypeTable {
    val named = namedTaxa(table.taxa, filter)
    val selected = if (filter.usesStats) keepByStats(table, filter, named, progress) else named
    if (selected.size == table.numberOfTaxa && selected.withIndex().all { (index, taxon) -> index == taxon }) {
        return table
    }

    val numberOfSites = table.numberOfSites
    val genotypes = ByteArray(numberOfSites * selected.size)
    var offset = 0
    for (block in table.siteBlocks(taxa = selected)) {
        for (site in block.sites) {
            block.copySiteGenotypes(site, genotypes, offset)
            offset += selected.size
        }
        progress.report(block.sites.last + 1, numberOfSites, "Filtering taxa")
    }
    val taxa = TaxaList(selected.map { table.taxa[it] })
    return InMemoryGenotypeTable.wrap(taxa, table.positions, genotypes)
}

/**
 * The spec for [filterTaxa], registered in [net.maizegenetics.tassel.ops.TasselOperations].
 *
 * TASSEL 5's `taxaList` is missing because the operation DSL has no string-list parameter type, so naming taxa means
 * calling [filterTaxa] with a [TaxaFilter] directly. [includeTaxa] is declared anyway, so the legacy name resolves;
 * like TASSEL 5, it does nothing without a taxon list.
 *
 * @see net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin
 */
val filterTaxaOperation = operation<GenotypeTable, GenotypeTable>(
    id = "filter-taxa",
    legacyPlugin = "FilterTaxaBuilderPlugin",
    help = "Keeps a subset of a genotype table's taxa",
) {
    val includeTaxa by boolean(
        default = true,
        legacy = "includeTaxa",
        help = "Keep the named taxa when true, drop them when false",
    )
    val minNotMissing by double(
        default = 0.0,
        range = 0.0..1.0,
        legacy = "minNotMissing",
        help = "Minimum proportion of sites with a called genotype",
    )
    val minHeterozygous by double(
        default = 0.0,
        range = 0.0..1.0,
        legacy = "minHeterozygous",
        help = "Minimum proportion of heterozygous sites",
    )
    val maxHeterozygous by double(
        default = 1.0,
        range = 0.0..1.0,
        legacy = "maxHeterozygous",
        help = "Maximum proportion of heterozygous sites",
    )

    run { input, progress ->
        filterTaxa(
            input,
            TaxaFilter(
                includeTaxa = includeTaxa(),
                minNotMissing = minNotMissing(),
                minHeterozygous = minHeterozygous(),
                maxHeterozygous = maxHeterozygous(),
            ),
            progress,
        )
    }
}

/** Taxon indices left by [TaxaFilter.taxaNames], in TASSEL 5's order. */
private fun namedTaxa(taxa: TaxaList, filter: TaxaFilter): IntArray {
    val names = filter.taxaNames ?: return IntArray(taxa.size) { it }
    return if (filter.includeTaxa) {
        // TASSEL 5 keeps the list's order here, and drops names the table does not have.
        names.map { taxa.indexOf(it) }.filter { it >= 0 }.toIntArray()
    } else {
        val dropped = names.toSet()
        taxa.indices.filter { taxa[it].name !in dropped }.toIntArray()
    }
}

/**
 * The entries of [selected] whose taxon passes [TaxaFilter]'s statistics thresholds, in the same order.
 *
 * The counts are accumulated over site blocks (ADR-0012) rather than by reading each taxon's whole row. They are the
 * two [net.maizegenetics.tassel.genotype.AlleleStats] counts TASSEL 5's taxon filter uses, so
 * `FilterTaxaTest.statisticsMatchAlleleStats` checks them against `GenotypeTable.taxonStats`.
 */
private fun keepByStats(
    table: GenotypeTable,
    filter: TaxaFilter,
    selected: IntArray,
    progress: ProgressReporter,
): IntArray {
    val numberOfSites = table.numberOfSites
    val unknown = IntArray(selected.size)
    val heterozygous = IntArray(selected.size)
    val buffer = ByteArray(selected.size)
    for (block in table.siteBlocks(taxa = selected)) {
        for (site in block.sites) {
            block.copySiteGenotypes(site, buffer)
            for (taxon in buffer.indices) {
                val genotype = buffer[taxon]
                if (genotype == Genotype.UNKNOWN) {
                    unknown[taxon]++
                    continue
                }
                val first = Genotype.firstAllele(genotype)
                val second = Genotype.secondAllele(genotype)
                if (first != second && first != Allele.UNKNOWN && second != Allele.UNKNOWN) heterozygous[taxon]++
            }
        }
        progress.report(block.sites.last + 1, numberOfSites, "Counting taxa")
    }

    val sites = numberOfSites.toDouble()
    return selected.filterIndexed { index, _ ->
        val notMissing = (numberOfSites - unknown[index]).toDouble() / sites
        if (notMissing < filter.minNotMissing) return@filterIndexed false
        val proportion = heterozygous[index].toDouble() / sites
        proportion >= filter.minHeterozygous && proportion <= filter.maxHeterozygous
    }.toIntArray()
}
