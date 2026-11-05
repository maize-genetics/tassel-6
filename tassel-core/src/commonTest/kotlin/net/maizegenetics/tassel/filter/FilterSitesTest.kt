package net.maizegenetics.tassel.filter

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readString
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.ops.TasselOperations
import net.maizegenetics.tassel.testing.SuspendLazy
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.assertMatchesGolden
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertSame
import kotlin.test.assertTrue

/**
 * `filter-sites` against the six TASSEL 5 `FilterSiteBuilderPlugin` parameter sets in
 * `testdata/golden/filter-sites/` (ADR-0009, `testdata/golden/MANIFEST`). Each golden is a HapMap export compared
 * exactly, so every test filters the tutorial file, writes it with [HapMap.write], and compares the text.
 */
class FilterSitesTest {
    @Test
    fun keepsSitesWithAtLeastMinCountTaxa() = runTest {
        assertGolden("min-count-270", SiteFilter(minCount = 270), sites = 2033)
    }

    @Test
    fun keepsSitesInsideTheMinorAlleleFrequencyRange() = runTest {
        assertGolden("maf-0.05-0.4", SiteFilter(minMaf = 0.05, maxMaf = 0.4), sites = 2046)
    }

    @Test
    fun keepsARangeOfSiteIndices() = runTest {
        val filter = SiteFilter(rangeFilterType = SiteRangeFilterType.SITES, startSite = 100, endSite = 599)
        assertGolden("site-range-100-599", filter, sites = 500)
    }

    @Test
    fun keepsARangeOfChromosomePositions() = runTest {
        val filter = SiteFilter(
            rangeFilterType = SiteRangeFilterType.POSITIONS,
            startChromosome = "5",
            startPosition = 1000000,
            endChromosome = "6",
            endPosition = 50000000,
        )
        assertGolden("position-range-5-6", filter, sites = 379)
    }

    @Test
    fun removingMinorSnpStatesKeepsEverySite() = runTest {
        assertGolden("remove-minor-snp-states", SiteFilter(removeMinorSnpStates = true), sites = 3093)
    }

    @Test
    fun combinesCountFrequencyAndHeterozygosity() = runTest {
        val filter = SiteFilter(minCount = 250, minMaf = 0.1, maxHeterozygous = 0.02)
        assertGolden("combined", filter, sites = 1953)
    }

    @Test
    fun removingMinorSnpStatesOnlyChangesGenotypes() = runTest {
        val filtered = filterSites(tutorial(), SiteFilter(removeMinorSnpStates = true))
        assertEquals(tutorial().positions.map { it.name }, filtered.positions.map { it.name })
        assertEquals(tutorial().taxa.names, filtered.taxa.names)
        var changed = 0
        for (site in 0 until tutorial().numberOfSites) {
            val before = tutorial().siteGenotypes(site)
            val after = filtered.siteGenotypes(site)
            before.indices.forEach { taxon ->
                if (before[taxon] != after[taxon]) {
                    assertEquals(Genotype.UNKNOWN, after[taxon], "site $site taxon $taxon")
                    changed++
                }
            }
        }
        assertTrue(changed > 0, "removeMinorSnpStates changed no genotype")
    }

    @Test
    fun aFilterThatKeepsEverythingReturnsTheSameTable() = runTest {
        assertSame(tutorial(), filterSites(tutorial(), SiteFilter()))
    }

    @Test
    fun aSiteRangeOutsideTheTableFails() = runTest {
        val filter = SiteFilter(rangeFilterType = SiteRangeFilterType.SITES, startSite = 3000, endSite = 3093)
        assertFailsWith<IllegalArgumentException> { filterSites(tutorial(), filter) }
    }

    @Test
    fun anUnknownChromosomeFails() = runTest {
        val filter = SiteFilter(
            rangeFilterType = SiteRangeFilterType.POSITIONS,
            startChromosome = "11",
            endChromosome = "11",
        )
        assertFailsWith<IllegalArgumentException> { filterSites(tutorial(), filter) }
    }

    @Test
    fun invalidThresholdsFail() = runTest {
        assertFailsWith<IllegalArgumentException> { SiteFilter(minCount = -1) }
        assertFailsWith<IllegalArgumentException> { SiteFilter(minMaf = 0.4, maxMaf = 0.1) }
        assertFailsWith<IllegalArgumentException> { SiteFilter(minMaf = 1.5) }
        assertFailsWith<IllegalArgumentException> { SiteFilter(maxHeterozygous = -0.1) }
        assertFailsWith<IllegalArgumentException> {
            SiteFilter(rangeFilterType = SiteRangeFilterType.SITES, startSite = 10, endSite = 5)
        }
        assertFailsWith<IllegalArgumentException> {
            SiteFilter(rangeFilterType = SiteRangeFilterType.POSITIONS, startChromosome = "1")
        }
    }

    @Test
    fun specCarriesTassel5Names() = runTest {
        assertEquals("filter-sites", filterSitesOperation.id)
        assertEquals("FilterSiteBuilderPlugin", filterSitesOperation.legacyPlugin)
        assertSame(filterSitesOperation, TasselOperations.registry["filter-sites"])
        assertSame(filterSitesOperation, TasselOperations.registry.byLegacyPlugin("FilterSiteBuilderPlugin"))
        val legacy = filterSitesOperation.parameters.associate { checkNotNull(it.legacy) to it.name }
        assertEquals(
            mapOf(
                "siteRangeFilterType" to "rangeFilterType",
                "startSite" to "startSite",
                "endSite" to "endSite",
                "startChr" to "startChr",
                "endChr" to "endChr",
                "startPos" to "startPosition",
                "endPos" to "endPosition",
                "siteMinCount" to "minCount",
                "siteMinAlleleFreq" to "minMaf",
                "siteMaxAlleleFreq" to "maxMaf",
                "minHeterozygous" to "minHeterozygous",
                "maxHeterozygous" to "maxHeterozygous",
                "removeMinorSNPStates" to "removeMinorSnpStates",
            ),
            legacy,
        )
    }

    @Test
    fun specRunsTheCombinedFilterFromCommandLineText() = runTest {
        val values = mapOf("minCount" to "250", "minMaf" to "0.1", "maxHeterozygous" to "0.02")
        assertEquals(1953, filterSitesOperation.run(tutorial(), values).numberOfSites)
    }

    @Test
    fun specInfersTheSiteRangeTypeFromTheSuppliedSites() = runTest {
        val filtered = filterSitesOperation.run(tutorial(), mapOf("startSite" to "100", "endSite" to "599"))
        assertEquals(500, filtered.numberOfSites)
    }

    /**
     * The `position-range-5-6` golden through the spec, which needs the string parameters ADR-0018 added. The values
     * are the command-line text of the MANIFEST command, and the result is compared to the golden, not just counted.
     */
    @Test
    fun specFiltersByChromosomePositionAgainstTheGolden() = runTest {
        val values = mapOf(
            "siteRangeFilterType" to "POSITIONS",
            "startChr" to "5",
            "startPos" to "1000000",
            "endChr" to "6",
            "endPos" to "50000000",
        ).mapKeys { (legacy, _) -> checkNotNull(filterSitesOperation.parameterByLegacyName(legacy)).name }

        val filtered = filterSitesOperation.run(tutorial(), values)
        assertEquals(379, filtered.numberOfSites)
        assertMatchesGolden(
            "golden/filter-sites/position-range-5-6.hmp.txt",
            Buffer().also { HapMap.write(filtered, it) }.readString(),
        )
    }

    @Test
    fun specInfersThePositionRangeTypeFromTheSuppliedChromosomes() = runTest {
        val values = mapOf(
            "startChr" to "5",
            "startPosition" to "1000000",
            "endChr" to "6",
            "endPosition" to "50000000",
        )
        assertEquals(379, filterSitesOperation.run(tutorial(), values).numberOfSites)
    }

    @Test
    fun specNeedsBothChromosomesToFilterByPosition() = runTest {
        assertFailsWith<IllegalArgumentException> {
            filterSitesOperation.run(tutorial(), mapOf("startPosition" to "1000000"))
        }
        assertFailsWith<IllegalArgumentException> {
            filterSitesOperation.run(tutorial(), mapOf("startChr" to "5"))
        }
    }

    private suspend fun assertGolden(name: String, filter: SiteFilter, sites: Int) {
        val filtered = filterSites(tutorial(), filter)
        assertEquals(sites, filtered.numberOfSites, "number of sites")
        assertEquals(281, filtered.numberOfTaxa, "number of taxa")
        val written = Buffer().also { HapMap.write(filtered, it) }.readString()
        assertMatchesGolden("golden/filter-sites/$name.hmp.txt", written)
    }

    private companion object {
        /** The tutorial file is read once and shared: parsing it per test is slow on Node and Native. */
        val tutorial = SuspendLazy {
            HapMap.read(Buffer().also { it.write(TestData.readBytes("tutorial/mdp_genotype.hmp.txt")) })
        }
    }
}
