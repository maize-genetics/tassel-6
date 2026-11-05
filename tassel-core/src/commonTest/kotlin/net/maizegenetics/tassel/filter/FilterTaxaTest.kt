package net.maizegenetics.tassel.filter

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readString
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.ops.TasselOperations
import net.maizegenetics.tassel.testing.SuspendLazy
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.assertMatchesGolden
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertSame

/**
 * `filter-taxa` against the four TASSEL 5 `FilterTaxaBuilderPlugin` parameter sets in
 * `testdata/golden/filter-taxa/` (ADR-0009, `testdata/golden/MANIFEST`). Each golden is a HapMap export compared
 * exactly, so every test filters the tutorial file, writes it with [HapMap.write], and compares the text.
 */
class FilterTaxaTest {
    @Test
    fun keepsTaxaWithEnoughCalledSites() = runTest {
        assertGolden("min-not-missing-0.95", TaxaFilter(minNotMissing = 0.95), taxa = 242)
    }

    @Test
    fun keepsTaxaInsideTheHeterozygosityRange() = runTest {
        val filter = TaxaFilter(minHeterozygous = 0.005, maxHeterozygous = 0.02)
        assertGolden("het-0.005-0.02", filter, taxa = 268)
    }

    @Test
    fun keepsTheNamedTaxaInTheOrderTheyWereNamed() = runTest {
        val filter = TaxaFilter(taxaNames = NAMES, includeTaxa = true)
        val filtered = assertGolden("include-list", filter, taxa = 4)
        assertEquals(NAMES, filtered.taxa.names)
    }

    @Test
    fun dropsTheNamedTaxaAndKeepsTableOrder() = runTest {
        val filter = TaxaFilter(taxaNames = NAMES, includeTaxa = false)
        val filtered = assertGolden("exclude-list", filter, taxa = 277)
        assertEquals(tutorial().taxa.names.filter { it !in NAMES }, filtered.taxa.names)
    }

    @Test
    fun namesTheTableDoesNotHaveAreIgnored() = runTest {
        val filter = TaxaFilter(taxaNames = listOf("W22", "not-a-taxon", "B73"), includeTaxa = true)
        assertEquals(listOf("W22", "B73"), filterTaxa(tutorial(), filter).taxa.names)
    }

    @Test
    fun aFilterThatKeepsEverythingReturnsTheSameTable() = runTest {
        assertSame(tutorial(), filterTaxa(tutorial(), TaxaFilter()))
        assertSame(tutorial(), filterTaxa(tutorial(), TaxaFilter(taxaNames = tutorial().taxa.names)))
    }

    /**
     * The taxon counters in `FilterTaxa.kt` accumulate over site blocks instead of reading each taxon's whole row, so
     * this pins them to the `AlleleStats` values TASSEL 5's `Stats` returns.
     */
    @Test
    fun statisticsMatchAlleleStats() = runTest {
        val keptByStats = filterTaxa(tutorial(), TaxaFilter(minNotMissing = 0.95)).taxa.names.toSet()
        val keptByAlleleStats = tutorial().taxa.indices
            .filter { tutorial().taxonStats(it).proportionNotMissing >= 0.95 }
            .map { tutorial().taxa[it].name }
            .toSet()
        assertEquals(keptByAlleleStats, keptByStats)

        val keptByHet = filterTaxa(tutorial(), TaxaFilter(minHeterozygous = 0.005, maxHeterozygous = 0.02))
            .taxa.names.toSet()
        val keptByAlleleStatsHet = tutorial().taxa.indices
            .filter { tutorial().taxonStats(it).proportionHeterozygous in 0.005..0.02 }
            .map { tutorial().taxa[it].name }
            .toSet()
        assertEquals(keptByAlleleStatsHet, keptByHet)
    }

    @Test
    fun invalidThresholdsFail() = runTest {
        assertFailsWith<IllegalArgumentException> { TaxaFilter(minNotMissing = 1.5) }
        assertFailsWith<IllegalArgumentException> { TaxaFilter(minHeterozygous = -0.1) }
        assertFailsWith<IllegalArgumentException> { TaxaFilter(maxHeterozygous = 2.0) }
        assertFailsWith<IllegalArgumentException> { TaxaFilter(taxaNames = listOf("B73", "B73")) }
    }

    @Test
    fun specCarriesTassel5Names() = runTest {
        assertEquals("filter-taxa", filterTaxaOperation.id)
        assertEquals("FilterTaxaBuilderPlugin", filterTaxaOperation.legacyPlugin)
        assertSame(filterTaxaOperation, TasselOperations.registry["filter-taxa"])
        assertSame(filterTaxaOperation, TasselOperations.registry.byLegacyPlugin("FilterTaxaBuilderPlugin"))
        assertEquals(
            mapOf(
                "taxaList" to "taxaList",
                "includeTaxa" to "includeTaxa",
                "minNotMissing" to "minNotMissing",
                "minHeterozygous" to "minHeterozygous",
                "maxHeterozygous" to "maxHeterozygous",
            ),
            filterTaxaOperation.parameters.associate { checkNotNull(it.legacy) to it.name },
        )
        // maxHeterozygous means the proportion of heterozygous sites here and of heterozygous taxa in filter-sites.
        val maxHeterozygous = checkNotNull(filterTaxaOperation.parameterByLegacyName("maxHeterozygous"))
        assertEquals("maxHeterozygous", maxHeterozygous.name)
        assertEquals(1.0, maxHeterozygous.default as Double)
    }

    @Test
    fun specRunsTheHeterozygosityFilterFromCommandLineText() = runTest {
        val values = mapOf("minHeterozygous" to "0.005", "maxHeterozygous" to "0.02")
        assertEquals(268, filterTaxaOperation.run(tutorial(), values).numberOfTaxa)
    }

    /**
     * The `include-list` and `exclude-list` goldens through the spec, which needs the string-list parameter ADR-0018
     * added. `taxaList` is the one comma-separated value TASSEL 5 takes on the command line.
     */
    @Test
    fun specSelectsTaxaByNameAgainstTheGoldens() = runTest {
        val text = NAMES.joinToString(",")
        assertEquals("W22,B73,CML247,33-16", text)

        val included = filterTaxaOperation.run(tutorial(), mapOf("taxaList" to text, "includeTaxa" to "true"))
        assertEquals(NAMES, included.taxa.names)
        assertMatchesGolden(
            "golden/filter-taxa/include-list.hmp.txt",
            Buffer().also { HapMap.write(included, it) }.readString(),
        )

        val excluded = filterTaxaOperation.run(tutorial(), mapOf("taxaList" to text, "includeTaxa" to "false"))
        assertEquals(277, excluded.numberOfTaxa)
        assertMatchesGolden(
            "golden/filter-taxa/exclude-list.hmp.txt",
            Buffer().also { HapMap.write(excluded, it) }.readString(),
        )
    }

    @Test
    fun specTakesATaxonListAsAlreadyTypedValuesToo() = runTest {
        assertEquals(NAMES, filterTaxaOperation.run(tutorial(), mapOf("taxaList" to NAMES)).taxa.names)
    }

    private suspend fun assertGolden(name: String, filter: TaxaFilter, taxa: Int): GenotypeTable {
        val filtered = filterTaxa(tutorial(), filter)
        assertEquals(taxa, filtered.numberOfTaxa, "number of taxa")
        assertEquals(3093, filtered.numberOfSites, "number of sites")
        val written = Buffer().also { HapMap.write(filtered, it) }.readString()
        assertMatchesGolden("golden/filter-taxa/$name.hmp.txt", written)
        return filtered
    }

    private companion object {
        /** The order TASSEL 5 was given on the command line, which the `include-list` golden keeps. */
        val NAMES = listOf("W22", "B73", "CML247", "33-16")

        /** The tutorial file is read once and shared: parsing it per test is slow on Node and Native. */
        val tutorial = SuspendLazy {
            HapMap.read(Buffer().also { it.write(TestData.readBytes("tutorial/mdp_genotype.hmp.txt")) })
        }
    }
}
