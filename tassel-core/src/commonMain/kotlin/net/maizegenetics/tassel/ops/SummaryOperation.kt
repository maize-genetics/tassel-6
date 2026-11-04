package net.maizegenetics.tassel.ops

import net.maizegenetics.tassel.analysis.GenotypeSummary
import net.maizegenetics.tassel.analysis.genotypeSummary
import net.maizegenetics.tassel.genotype.GenotypeTable

/**
 * The genotype summary: an overview table, an allele table, a site table, and a taxa table (ADR-0002).
 *
 * @see net.maizegenetics.analysis.data.GenotypeSummaryPlugin
 */
val summaryOperation = operation<GenotypeTable, GenotypeSummary>(
    id = "summary",
    legacyPlugin = "GenotypeSummaryPlugin",
    help = "Counts and proportions over a genotype table, by site and by taxon",
) {
    val overview by boolean(default = true, legacy = "overview", help = "Get Overview Report")
    val siteSummary by boolean(default = true, legacy = "siteSummary", help = "Get Site Summary")
    val taxaSummary by boolean(default = true, legacy = "taxaSummary", help = "Get Taxa Summary")
    run { input, progress ->
        genotypeSummary(input, overview(), siteSummary(), taxaSummary(), progress)
    }
}
