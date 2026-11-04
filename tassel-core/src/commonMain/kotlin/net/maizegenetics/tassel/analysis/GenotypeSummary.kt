package net.maizegenetics.tassel.analysis

import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.AlleleStats
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.table.ResultTable
import net.maizegenetics.tassel.table.ResultTableBuilder
import net.maizegenetics.tassel.table.resultTable

/**
 * The tables [genotypeSummary] produced, in the order TASSEL 5 returns them. A table is null when its report was
 * turned off; [overview] and [alleles] are both produced or both absent.
 */
class GenotypeSummary internal constructor(
    /** Counts and proportions over the whole table, titled `Overall Summary`. */
    val overview: ResultTable?,
    /** Genotype and major/minor allele pair counts, titled `Allele Summary`. */
    val alleles: ResultTable?,
    /** One row per site, titled `Site Summary`. */
    val sites: ResultTable?,
    /** One row per taxon, titled `Taxa Summary`. */
    val taxa: ResultTable?,
) {
    /** The tables that were produced, in TASSEL 5's order: overview, alleles, sites, taxa. */
    val tables: List<ResultTable> = listOfNotNull(overview, alleles, sites, taxa)
}

/**
 * Summarizes [genotypes] as the four tables TASSEL 5's genotype summary produces.
 *
 * Each report can be turned off, exactly as the TASSEL 5 plugin parameters do. The [overview] table reuses the pass
 * the site summary already made, and TASSEL 5 counts missing gametes differently in its fallback pass, so turning
 * [siteSummary] off can change the overview on data that holds rare or undefined allele codes. That quirk is
 * reproduced; the tutorial dataset has no such codes, so the goldens are the same either way.
 *
 * @see net.maizegenetics.analysis.data.GenotypeSummaryPlugin
 */
fun genotypeSummary(
    genotypes: GenotypeTable,
    overview: Boolean = true,
    siteSummary: Boolean = true,
    taxaSummary: Boolean = true,
    progress: ProgressReporter = ProgressReporter.None,
): GenotypeSummary {
    val totals = GenotypeSummaryTotals()
    val sites = if (siteSummary) siteTable(genotypes, totals, progress) else null
    val taxa = if (taxaSummary) taxaTable(genotypes, progress) else null
    var overviewTable: ResultTable? = null
    var alleleTable: ResultTable? = null
    if (overview) {
        val counts = overviewCounts(genotypes, totals, accumulate = !siteSummary, progress = progress)
        overviewTable = overviewTable(genotypes, totals, counts)
        alleleTable = alleleTable(genotypes, counts)
    }
    return GenotypeSummary(overviewTable, alleleTable, sites, taxa)
}

/**
 * Running totals the overview table needs, filled in by whichever pass over the sites happens to run.
 *
 * TASSEL 5 keeps these as plugin fields and fills them in either from the site summary or, when the site summary is
 * off, from a second loop that defines one of them differently: [gametesMissing] then counts allele codes outside A,
 * C, G, T, insert, and gap (rare and undefined) as present rather than missing.
 */
private class GenotypeSummaryTotals {
    /** Gametes with no called allele, summed over sites. */
    var gametesMissing: Long = 0

    /** Genotypes whose two allele codes differ, including one-sided missing calls, summed over sites. */
    var heterozygous: Long = 0

    /** Mean of the per-site minor allele frequency. */
    var averageMinorAlleleFrequency: Double = 0.0
}

/** The counts over every genotype that the overview and allele tables share. */
private class OverviewCounts(
    /** IUPAC genotype code to the number of times it occurs, most frequent first. */
    val genotypes: List<Pair<String, Long>>,
    /** `major:minor` allele pair to the number of sites with that pair, most frequent first. */
    val majorMinorPairs: List<Pair<String, Long>>,
) {
    /** Genotypes with both alleles missing. */
    val missing: Long = genotypes.firstOrNull { it.first == UNKNOWN_GENOTYPE }?.second ?: 0L
}

private val SITE_COLUMNS_BEFORE_ALLELES = listOf(
    "Site Number", "Site Name", "Chromosome", "Physical Position", "Number of Taxa", "Ref", "Alt",
    "Major Allele", "Major Allele Gametes", "Major Allele Proportion", "Major Allele Frequency",
    "Minor Allele", "Minor Allele Gametes", "Minor Allele Proportion", "Minor Allele Frequency",
)

private val SITE_COLUMNS_AFTER_ALLELES = listOf(
    "Gametes Missing",
    "Proportion Missing",
    "Number Heterozygous",
    "Proportion Heterozygous",
    "Inbreeding Coefficient",
    "Inbreeding Coefficient Scaled by Missing",
)

private val TAXA_COLUMNS = listOf(
    "Taxa", "Taxa Name", "Number of Sites", "Gametes Missing", "Proportion Missing",
    "Number Heterozygous", "Proportion Heterozygous", "Inbreeding Coefficient",
    "Inbreeding Coefficient Scaled by Missing",
)

/** TASSEL 5 writes a placeholder where the inbreeding coefficients will go; see `GenotypeSummaryPlugin`. */
private const val TBD = "TBD"
private const val INBREEDING_COEFFICIENT = "Inbreeding Coefficient"
private const val INBREEDING_COEFFICIENT_SCALED = "ICSBM"

/** The allele columns of an absent allele, as TASSEL 5's `NA`, `ZERO_INT`, and `ZERO_DOUBLE` constants. */
private const val NA = "NA"
private const val ZERO_INT: Int = 0
private const val ZERO_DOUBLE: Double = 0.0

/** TASSEL 5 `GenotypeTable.UNKNOWN_ALLELE_STR`, which is also the IUPAC code of a fully missing genotype. */
private const val UNKNOWN_GENOTYPE = "N"

/**
 * Maximum alleles a site can report, TASSEL 5's `NucleotideGenotypeCallTable.maxNumAlleles()`. A nucleotide table does
 * not retain rare alleles, so this is also the number of states [AlleleStats] counts.
 */
private const val MAX_ALLELES = Allele.NUCLEOTIDE_STATES

/**
 * One row per site, and the running totals the overview needs.
 *
 * @see net.maizegenetics.analysis.data.GenotypeSummaryPlugin.getSiteSummary
 */
private fun siteTable(
    genotypes: GenotypeTable,
    totals: GenotypeSummaryTotals,
    progress: ProgressReporter,
): ResultTable {
    val columnNames = buildList {
        addAll(SITE_COLUMNS_BEFORE_ALLELES)
        for (allele in 2 until MAX_ALLELES) {
            val heading = "Allele ${allele + 1}"
            add(heading)
            add("$heading Gametes")
            add("$heading Proportion")
            add("$heading Frequency")
        }
        addAll(SITE_COLUMNS_AFTER_ALLELES)
    }
    val numberOfSites = genotypes.numberOfSites
    val numberOfTaxa = genotypes.numberOfTaxa
    val totalGametes = numberOfTaxa * 2
    val builder = ResultTableBuilder("Site Summary", columnNames)
    for (block in genotypes.siteBlocks()) {
        for (site in block.sites) {
            val stats = block.siteStats(site)
            val position = genotypes.positions[site]
            val cells = ArrayList<Any?>(columnNames.size)
            cells += site
            cells += position.name
            cells += position.chromosome.name
            cells += position.position
            cells += numberOfTaxa
            // A Position carries no reference or alternate allele, so both read as unknown, as they do for a TASSEL 5
            // genotype table built from HapMap.
            cells += Nucleotide.alleleString(Allele.UNKNOWN)
            cells += Nucleotide.alleleString(Allele.UNKNOWN)
            val alleleGametes = stats.alleleGameteCount
            stats.alleles.forEachIndexed { index, allele ->
                val count = stats.alleleCounts[index]
                val frequency = count.toDouble() / alleleGametes.toDouble()
                cells += Nucleotide.alleleString(allele)
                cells += count
                cells += count.toDouble() / totalGametes.toDouble()
                cells += frequency
                if (index == 1) totals.averageMinorAlleleFrequency += frequency
            }
            repeat(MAX_ALLELES - stats.alleles.size) {
                cells += NA
                cells += ZERO_INT
                cells += ZERO_DOUBLE
                cells += ZERO_DOUBLE
            }
            val gametesMissing = totalGametes - alleleGametes
            totals.gametesMissing += gametesMissing.toLong()
            cells += gametesMissing
            cells += gametesMissing.toDouble() / totalGametes.toDouble()
            val heterozygous = heterozygousCount(stats)
            totals.heterozygous += heterozygous.toLong()
            cells += heterozygous
            cells += heterozygous.toDouble() / stats.nonMissingCount.toDouble()
            cells += TBD
            cells += TBD
            builder.row(cells)
        }
        progress.report(block.sites.last + 1, numberOfSites, "Site summary")
    }
    totals.averageMinorAlleleFrequency /= numberOfSites.toDouble()
    return builder.build()
}

/**
 * One row per taxon.
 *
 * @see net.maizegenetics.analysis.data.GenotypeSummaryPlugin.getTaxaSummary
 */
private fun taxaTable(genotypes: GenotypeTable, progress: ProgressReporter): ResultTable {
    val numberOfSites = genotypes.numberOfSites
    val numberOfTaxa = genotypes.numberOfTaxa
    val totalGametes = numberOfSites * 2
    val gametesNotMissing = IntArray(numberOfTaxa)
    val sitesNotMissing = IntArray(numberOfTaxa)
    val heterozygous = IntArray(numberOfTaxa)
    val siteGenotypes = ByteArray(numberOfTaxa)
    for (block in genotypes.siteBlocks()) {
        for (site in block.sites) {
            block.copySiteGenotypes(site, siteGenotypes)
            for (taxon in 0 until numberOfTaxa) {
                val genotype = siteGenotypes[taxon]
                val first = Genotype.firstAllele(genotype)
                val second = Genotype.secondAllele(genotype)
                if (first != Allele.UNKNOWN || second != Allele.UNKNOWN) {
                    sitesNotMissing[taxon]++
                    if (first != Allele.UNKNOWN) gametesNotMissing[taxon]++
                    if (second != Allele.UNKNOWN) gametesNotMissing[taxon]++
                }
                if (first != second) heterozygous[taxon]++
            }
        }
        progress.report(block.sites.last + 1, numberOfSites, "Taxa summary")
    }
    val builder = ResultTableBuilder("Taxa Summary", TAXA_COLUMNS)
    for (taxon in 0 until numberOfTaxa) {
        val gametesMissing = totalGametes - gametesNotMissing[taxon]
        builder.row(
            taxon,
            genotypes.taxa[taxon].name,
            numberOfSites,
            gametesMissing,
            gametesMissing.toDouble() / totalGametes.toDouble(),
            heterozygous[taxon],
            heterozygous[taxon].toDouble() / sitesNotMissing[taxon].toDouble(),
            INBREEDING_COEFFICIENT,
            INBREEDING_COEFFICIENT_SCALED,
        )
    }
    return builder.build()
}

/**
 * Counts every genotype by its IUPAC code and every site by its major/minor allele pair. When [accumulate] is true,
 * also fills in [totals] the way TASSEL 5 does when the site summary is off.
 *
 * @see net.maizegenetics.dna.snp.genotypecall.AbstractGenotypeCallTable.genoCounts
 * @see net.maizegenetics.dna.snp.genotypecall.AbstractGenotypeCallTable.majorMinorCounts
 */
private fun overviewCounts(
    genotypes: GenotypeTable,
    totals: GenotypeSummaryTotals,
    accumulate: Boolean,
    progress: ProgressReporter,
): OverviewCounts {
    val numberOfSites = genotypes.numberOfSites
    val numberOfTaxa = genotypes.numberOfTaxa
    val totalGametes = numberOfTaxa * 2
    val byGenotype = LongArray(GENOTYPE_CODES)
    val byAllelePair = LongArray(ALLELE_CODES * ALLELE_CODES)
    val siteGenotypes = ByteArray(numberOfTaxa)
    for (block in genotypes.siteBlocks()) {
        for (site in block.sites) {
            block.copySiteGenotypes(site, siteGenotypes)
            for (genotype in siteGenotypes) byGenotype[genotype.toInt() and 0xFF]++
            val stats = AlleleStats.of(siteGenotypes)
            when (stats.alleles.size) {
                0 -> Unit
                1 -> byAllelePair[pairIndex(stats.majorAllele, stats.majorAllele)]++
                else -> byAllelePair[pairIndex(stats.majorAllele, stats.minorAllele)]++
            }
            if (accumulate) {
                totals.gametesMissing += (totalGametes - stats.nonMissingGameteCount).toLong()
                totals.heterozygous += heterozygousCount(stats).toLong()
                totals.averageMinorAlleleFrequency += stats.minorAlleleFrequency
            }
        }
        progress.report(block.sites.last + 1, numberOfSites, "Overall summary")
    }
    if (accumulate) totals.averageMinorAlleleFrequency /= numberOfSites.toDouble()

    val genotypeCounts = LinkedHashMap<String, Long>()
    for (genotype in 0 until GENOTYPE_CODES) {
        if (byGenotype[genotype] == 0L) continue
        val code = Nucleotide.iupac(genotype.toByte()).toString()
        genotypeCounts[code] = (genotypeCounts[code] ?: 0L) + byGenotype[genotype]
    }
    val pairCounts = LinkedHashMap<String, Long>()
    for (pair in byAllelePair.indices) {
        if (byAllelePair[pair] == 0L) continue
        val first = Nucleotide.alleleString((pair / ALLELE_CODES).toByte())
        val second = Nucleotide.alleleString((pair % ALLELE_CODES).toByte())
        pairCounts["$first:$second"] = byAllelePair[pair]
    }
    // TASSEL 5 bubble-sorts both lists by count, descending, which is stable. The allele-pair list starts in
    // ascending (major, minor) code order there as it does here, so ties come out the same way. The genotype list
    // starts in Java HashMap order there, which no pure-Kotlin code can reproduce; ties here keep ascending genotype
    // code order instead. No golden has a tie.
    return OverviewCounts(
        genotypes = genotypeCounts.toList().sortedByDescending { it.second },
        majorMinorPairs = pairCounts.toList().sortedByDescending { it.second },
    )
}

/**
 * The whole-table counts and proportions.
 *
 * @see net.maizegenetics.analysis.data.GenotypeSummaryPlugin.getOverallSummary
 */
private fun overviewTable(
    genotypes: GenotypeTable,
    totals: GenotypeSummaryTotals,
    counts: OverviewCounts,
): ResultTable {
    val numberOfSites = genotypes.numberOfSites.toLong()
    val numberOfTaxa = genotypes.numberOfTaxa.toLong()
    val totalDiploids = numberOfSites * numberOfTaxa
    val diploidsNotMissing = totalDiploids - counts.missing
    val totalGametes = totalDiploids * 2L
    val gametesNotMissing = totalGametes - totals.gametesMissing
    return resultTable("Overall Summary", listOf("Stat Type", "Value")) {
        row("Number of Taxa", numberOfTaxa.toDouble())
        row("Number of Sites", numberOfSites.toDouble())
        row("Sites x Taxa", totalDiploids.toDouble())
        row("Number Not Missing", diploidsNotMissing.toDouble())
        row("Proportion Not Missing", diploidsNotMissing.toDouble() / totalDiploids.toDouble())
        row("Number Missing", counts.missing.toDouble())
        row("Proportion Missing", counts.missing.toDouble() / totalDiploids.toDouble())
        row("Number Gametes", totalGametes.toDouble())
        row("Gametes Not Missing", gametesNotMissing.toDouble())
        row("Proportion Gametes Not Missing", gametesNotMissing.toDouble() / totalGametes.toDouble())
        row("Gametes Missing", totals.gametesMissing.toDouble())
        row("Proportion Gametes Missing", totals.gametesMissing.toDouble() / totalGametes.toDouble())
        row("Number Heterozygous", totals.heterozygous.toDouble())
        row("Proportion Heterozygous", totals.heterozygous.toDouble() / totalDiploids.toDouble())
        row("Average Minor Allele Frequency", totals.averageMinorAlleleFrequency)
    }
}

/**
 * Genotype counts over the whole table, then the major/minor allele pair counts over the sites. The pair rows have no
 * frequency, which TASSEL 5 leaves as a null cell and writes as nothing.
 *
 * @see net.maizegenetics.analysis.data.GenotypeSummaryPlugin.getOverallSummary
 */
private fun alleleTable(genotypes: GenotypeTable, counts: OverviewCounts): ResultTable {
    val numberOfSites = genotypes.numberOfSites.toLong()
    val totalDiploids = numberOfSites * genotypes.numberOfTaxa.toLong()
    val diploidsNotMissing = totalDiploids - counts.missing
    val builder = ResultTableBuilder("Allele Summary", listOf("Alleles", "Number", "Proportion", "Frequency"))
    for ((code, count) in counts.genotypes) {
        builder.row(
            code,
            count,
            count.toDouble() / totalDiploids.toDouble(),
            count.toDouble() / diploidsNotMissing.toDouble(),
        )
    }
    for ((pair, count) in counts.majorMinorPairs) {
        builder.row(pair, count, count.toDouble() / numberOfSites.toDouble(), null)
    }
    return builder.build()
}

/**
 * Genotypes whose two allele codes differ, as TASSEL 5's `GenotypeTable.heterozygousCount(site)` counts them: a
 * one-sided missing call such as `A/N` is counted, which [AlleleStats.heterozygousCount] does not do.
 */
private fun heterozygousCount(stats: AlleleStats): Int = stats.heterozygousCount + stats.partlyUnknownCount

private const val GENOTYPE_CODES = 256
private const val ALLELE_CODES = 16

private fun pairIndex(first: Byte, second: Byte): Int =
    (first.toInt() and 0xF) * ALLELE_CODES + (second.toInt() and 0xF)
