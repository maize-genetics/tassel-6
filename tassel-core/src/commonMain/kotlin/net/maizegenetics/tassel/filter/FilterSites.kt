package net.maizegenetics.tassel.filter

import net.maizegenetics.tassel.genotype.AlleleStats
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.operation
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList

/** Which coordinates [SiteFilter] uses to pick a contiguous run of sites. */
enum class SiteRangeFilterType {
    /** Keep every site. */
    NONE,

    /** Keep sites by table index, `startSite` through `endSite`. */
    SITES,

    /** Keep sites by chromosome and physical position, `startChromosome`/`startPosition` through the end pair. */
    POSITIONS,
}

/**
 * Which sites a [filterSites] call keeps, and whether it rewrites genotypes.
 *
 * Defaults are TASSEL 5's, so `SiteFilter()` keeps the table unchanged. The filters are applied in TASSEL 5's order:
 * the range first, then the statistics thresholds (all inclusive at both ends), then [removeMinorSnpStates].
 *
 * @property rangeFilterType which coordinates [startSite]/[endSite] or the chromosome and position pairs use.
 * @property startSite first site index to keep, used when [rangeFilterType] is [SiteRangeFilterType.SITES].
 * @property endSite last site index to keep, inclusive.
 * @property startChromosome chromosome of the first site to keep, used with [SiteRangeFilterType.POSITIONS].
 * @property startPosition physical position of the first site to keep; `null` means the first site of
 *   [startChromosome]. A position with no site starts at the next site on that chromosome.
 * @property endChromosome chromosome of the last site to keep.
 * @property endPosition physical position of the last site to keep, inclusive; `null` means the last site of
 *   [endChromosome]. A position with no site ends at the previous site on that chromosome.
 * @property minCount minimum number of taxa with a called genotype at the site. TASSEL 5 compares
 *   `2 * minCount` against the site's non-missing gamete count, so one missing allele counts as half a taxon.
 * @property minMaf minimum minor allele frequency.
 * @property maxMaf maximum minor allele frequency.
 * @property minHeterozygous minimum proportion of heterozygous taxa.
 * @property maxHeterozygous maximum proportion of heterozygous taxa.
 * @property removeMinorSnpStates sets every genotype with an allele other than the site's major or minor allele to
 *   unknown. This changes genotypes without dropping sites.
 * @see net.maizegenetics.dna.snp.FilterSite
 */
data class SiteFilter(
    val rangeFilterType: SiteRangeFilterType = SiteRangeFilterType.NONE,
    val startSite: Int = 0,
    val endSite: Int = 0,
    val startChromosome: String? = null,
    val startPosition: Int? = null,
    val endChromosome: String? = null,
    val endPosition: Int? = null,
    val minCount: Int = 0,
    val minMaf: Double = 0.0,
    val maxMaf: Double = 1.0,
    val minHeterozygous: Double = 0.0,
    val maxHeterozygous: Double = 1.0,
    val removeMinorSnpStates: Boolean = false,
) {
    init {
        require(minCount >= 0) { "minCount must not be negative: $minCount" }
        require(minMaf in 0.0..1.0) { "minMaf must be between 0.0 and 1.0: $minMaf" }
        require(maxMaf in 0.0..1.0) { "maxMaf must be between 0.0 and 1.0: $maxMaf" }
        require(maxMaf >= minMaf) { "maxMaf $maxMaf is less than minMaf $minMaf" }
        require(minHeterozygous in 0.0..1.0) { "minHeterozygous must be between 0.0 and 1.0: $minHeterozygous" }
        require(maxHeterozygous in 0.0..1.0) { "maxHeterozygous must be between 0.0 and 1.0: $maxHeterozygous" }
        when (rangeFilterType) {
            SiteRangeFilterType.SITES -> require(startSite in 0..endSite) {
                "startSite $startSite must be between 0 and endSite $endSite"
            }

            SiteRangeFilterType.POSITIONS -> require(startChromosome != null && endChromosome != null) {
                "startChromosome and endChromosome are required when filtering by positions"
            }

            SiteRangeFilterType.NONE -> Unit
        }
    }

    /** True when a site's statistics decide whether it is kept, which is TASSEL 5's test for computing them. */
    internal val usesStats: Boolean
        get() = minCount != 0 || minMaf != 0.0 || maxMaf != 1.0 || minHeterozygous != 0.0 || maxHeterozygous != 1.0

    /** True when a site with these [stats] passes every statistics threshold. */
    internal fun keeps(stats: AlleleStats): Boolean {
        // TASSEL 5 compares gametes, not taxa: siteMinCount 270 keeps a site with 540 non-missing allele copies.
        if (stats.alleleGameteCount < minCount * 2) return false
        val maf = stats.minorAlleleFrequency
        if (maf < minMaf || maf > maxMaf) return false
        val heterozygous = stats.proportionHeterozygous
        return heterozygous >= minHeterozygous && heterozygous <= maxHeterozygous
    }
}

/**
 * Keeps the sites of [table] that [filter] selects, in table order, with all taxa.
 *
 * Returns [table] itself when [filter] would keep every site and change no genotype.
 *
 * @throws IllegalArgumentException if a site range falls outside the table or names a chromosome it does not have.
 * @see net.maizegenetics.analysis.filter.FilterBySites
 * @see net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin
 */
fun filterSites(
    table: GenotypeTable,
    filter: SiteFilter,
    progress: ProgressReporter = ProgressReporter.None,
): GenotypeTable {
    val range = siteRange(table, filter)
    val rewrites = filter.usesStats || filter.removeMinorSnpStates
    if (!rewrites && range.first == 0 && range.last == table.numberOfSites - 1) return table

    val numberOfTaxa = table.numberOfTaxa
    val keptPositions = ArrayList<Position>()
    val keptGenotypes = ArrayList<ByteArray>()
    for (block in table.siteBlocks()) {
        if (block.sites.last < range.first) continue
        if (block.sites.first > range.last) break
        for (site in block.sites) {
            if (site < range.first || site > range.last) continue
            val genotypes = block.siteGenotypes(site)
            if (rewrites) {
                val stats = AlleleStats.of(genotypes)
                if (filter.usesStats && !filter.keeps(stats)) continue
                // TASSEL 5 masks after subsetting sites, but a site's statistics do not depend on the other sites.
                if (filter.removeMinorSnpStates) maskMinorSnpStates(genotypes, stats)
            }
            keptPositions += table.positions[site]
            keptGenotypes += genotypes
        }
        progress.report(block.sites.last + 1, table.numberOfSites, "Filtering sites")
    }

    val combined = ByteArray(keptGenotypes.size * numberOfTaxa)
    keptGenotypes.forEachIndexed { site, genotypes -> genotypes.copyInto(combined, site * numberOfTaxa) }
    return InMemoryGenotypeTable.wrap(table.taxa, PositionList(keptPositions), combined)
}

/**
 * The spec for [filterSites], registered in [net.maizegenetics.tassel.ops.TasselOperations].
 *
 * @see net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin
 */
val filterSitesOperation = operation<GenotypeTable, GenotypeTable>(
    id = "filter-sites",
    legacyPlugin = "FilterSiteBuilderPlugin",
    help = "Keeps a subset of a genotype table's sites",
) {
    val rangeFilterType by enum(
        default = SiteRangeFilterType.NONE,
        legacy = "siteRangeFilterType",
        help = "Keep a run of sites by table index or by chromosome and position",
    )
    val startSite by int(default = 0, range = 0..Int.MAX_VALUE, legacy = "startSite", help = "First site to keep")
    val endSite by int(default = 0, range = 0..Int.MAX_VALUE, legacy = "endSite", help = "Last site to keep")
    val startChr by string(legacy = "startChr", help = "Chromosome of the first site to keep")
    val endChr by string(legacy = "endChr", help = "Chromosome of the last site to keep")
    val startPosition by int(
        range = 0..Int.MAX_VALUE,
        legacy = "startPos",
        help = "Physical position of the first site to keep",
    )
    val endPosition by int(
        range = 0..Int.MAX_VALUE,
        legacy = "endPos",
        help = "Physical position of the last site to keep",
    )
    val minCount by int(
        default = 0,
        range = 0..Int.MAX_VALUE,
        legacy = "siteMinCount",
        help = "Minimum number of taxa with a called genotype",
    )
    val minMaf by double(
        default = 0.0,
        range = 0.0..1.0,
        legacy = "siteMinAlleleFreq",
        help = "Minimum minor allele frequency",
    )
    val maxMaf by double(
        default = 1.0,
        range = 0.0..1.0,
        legacy = "siteMaxAlleleFreq",
        help = "Maximum minor allele frequency",
    )
    val minHeterozygous by double(
        default = 0.0,
        range = 0.0..1.0,
        legacy = "minHeterozygous",
        help = "Minimum proportion of heterozygous taxa",
    )
    val maxHeterozygous by double(
        default = 1.0,
        range = 0.0..1.0,
        legacy = "maxHeterozygous",
        help = "Maximum proportion of heterozygous taxa",
    )
    val removeMinorSnpStates by boolean(
        legacy = "removeMinorSNPStates",
        help = "Set genotypes with an allele other than the site's major or minor allele to unknown",
    )

    run { input, progress ->
        // TASSEL 5 FilterSiteBuilderPlugin.preProcessParameters infers the range type from the supplied coordinates.
        val inferred = when {
            rangeFilterType() != SiteRangeFilterType.NONE -> rangeFilterType()
            startChr() != null || endChr() != null || startPosition() != null || endPosition() != null ->
                SiteRangeFilterType.POSITIONS
            startSite() != 0 || endSite() != 0 -> SiteRangeFilterType.SITES
            else -> SiteRangeFilterType.NONE
        }
        filterSites(
            input,
            SiteFilter(
                rangeFilterType = inferred,
                startSite = startSite(),
                endSite = endSite(),
                startChromosome = startChr(),
                startPosition = startPosition(),
                endChromosome = endChr(),
                endPosition = endPosition(),
                minCount = minCount(),
                minMaf = minMaf(),
                maxMaf = maxMaf(),
                minHeterozygous = minHeterozygous(),
                maxHeterozygous = maxHeterozygous(),
                removeMinorSnpStates = removeMinorSnpStates(),
            ),
            progress,
        )
    }
}

/** The sites [filter] selects, as an inclusive index range. */
private fun siteRange(table: GenotypeTable, filter: SiteFilter): IntRange {
    val numberOfSites = table.numberOfSites
    val positions = table.positions
    return when (filter.rangeFilterType) {
        SiteRangeFilterType.NONE -> 0..numberOfSites - 1

        SiteRangeFilterType.SITES -> {
            require(filter.endSite < numberOfSites) {
                "endSite ${filter.endSite} is outside a table of $numberOfSites sites"
            }
            filter.startSite..filter.endSite
        }

        SiteRangeFilterType.POSITIONS -> {
            val startChromosome = Chromosome.of(filter.startChromosome!!)
            val endChromosome = Chromosome.of(filter.endChromosome!!)
            val start = if (filter.startPosition == null) {
                val sites = positions.siteRange(startChromosome)
                require(!sites.isEmpty()) { "The table has no sites on chromosome $startChromosome" }
                sites.first
            } else {
                val site = positions.siteOfPhysicalPosition(startChromosome, filter.startPosition)
                // No site at that position: TASSEL 5 starts at the site the position would be inserted before.
                if (site >= 0) site else minOf(-(site + 1), numberOfSites - 1)
            }
            val end = if (filter.endPosition == null) {
                val sites = positions.siteRange(endChromosome)
                require(!sites.isEmpty()) { "The table has no sites on chromosome $endChromosome" }
                sites.last
            } else {
                val site = positions.siteOfPhysicalPosition(endChromosome, filter.endPosition)
                // No site at that position: TASSEL 5 ends at the site before the insertion point.
                if (site >= 0) site else minOf(-site - 2, numberOfSites - 1)
            }
            require(start in 0..end && end < numberOfSites) {
                "Positions ${filter.startChromosome}:${filter.startPosition} to " +
                    "${filter.endChromosome}:${filter.endPosition} select sites $start to $end, which is not a " +
                    "range inside a table of $numberOfSites sites"
            }
            start..end
        }
    }
}

/**
 * Sets every genotype of the site in [genotypes] to unknown unless both of its alleles are the site's major or minor
 * allele.
 *
 * The comparison keeps TASSEL 5's `MaskMatrixBuilder.getInstanceRemoveMinorSNPs`, including the sign extension in its
 * `genotype >>> 4`: a first allele of 8 or above (rare or unknown) never matches, so a genotype such as `N/A` is
 * masked even at a site whose minor allele is unknown.
 */
private fun maskMinorSnpStates(genotypes: ByteArray, stats: AlleleStats) {
    val major = stats.majorAllele.toInt()
    val minor = stats.minorAllele.toInt()
    for (index in genotypes.indices) {
        val genotype = genotypes[index].toInt()
        val second = genotype and 0xF
        val masked = if (second != major && second != minor) {
            true
        } else {
            val first = genotype ushr 4
            first != major && first != minor
        }
        if (masked) genotypes[index] = Genotype.UNKNOWN
    }
}
