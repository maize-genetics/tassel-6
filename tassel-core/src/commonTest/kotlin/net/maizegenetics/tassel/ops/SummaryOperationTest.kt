package net.maizegenetics.tassel.ops

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.analysis.GenotypeSummary
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNull
import kotlin.test.assertTrue

/** The `summary` spec: its TASSEL 5 names, its defaults, and that running it returns the four tables. */
class SummaryOperationTest {
    @Test
    fun carriesTheTassel5PluginAndParameterNames() {
        assertEquals("GenotypeSummaryPlugin", summaryOperation.legacyPlugin)
        assertEquals(GenotypeTable::class, summaryOperation.inputType)
        assertEquals(GenotypeSummary::class, summaryOperation.outputType)
        assertEquals(listOf("overview", "siteSummary", "taxaSummary"), summaryOperation.parameters.map { it.name })
        summaryOperation.parameters.forEach {
            assertEquals(it.name, it.legacy, "TASSEL 5 name of '${it.name}'")
            assertEquals(true, it.default, "default of '${it.name}'")
            assertTrue(it.help.isNotEmpty(), "help of '${it.name}'")
        }
        assertEquals(summaryOperation, TasselOperations.registry["summary"])
    }

    @Test
    fun runsWithDefaultsAndWithEveryReportOff() = runTest {
        val all = summaryOperation.run(table)
        assertEquals(
            listOf("Overall Summary", "Allele Summary", "Site Summary", "Taxa Summary"),
            all.tables.map { it.title },
        )
        val off = mapOf("overview" to "false", "siteSummary" to "false", "taxaSummary" to "false")
        val none = summaryOperation.run(table, off)
        assertTrue(none.tables.isEmpty())
        assertNull(none.sites)
    }

    private companion object {
        val table: GenotypeTable = InMemoryGenotypeTable.of(
            TaxaList.of("T0", "T1"),
            PositionList(listOf(Position(Chromosome.of("1"), 100, "s1"))),
            byteArrayOf(Nucleotide.genotype("AA"), Nucleotide.genotype("AC")),
        )
    }
}
