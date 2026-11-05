package net.maizegenetics.tassel.app.grid

import net.maizegenetics.tassel.app.testing.cellGenotype
import net.maizegenetics.tassel.app.testing.testTable
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * What a frame of the grid actually reads out of a [GenotypeTable]. These are the tests that pin down the claim that
 * a frame is cheap: the reads are counted, and they come out equal to the number of visible cells wherever the grid
 * is scrolled and however big the table is.
 */
class ReadWindowTest {
    @Test
    fun readsEveryVisibleCellExactlyOnceAndNothingElse() {
        val table = RecordingGenotypeTable(testTable(numberOfTaxa = 50, numberOfSites = 400))
        val window = GridWindow(firstTaxon = 4, taxonCount = 3, firstSite = 100, siteCount = 5)

        table.readWindow(window)

        val expected = window.sites.flatMap { site -> window.taxa.map { taxon -> taxon to site } }
        assertEquals(expected.size, table.reads.size)
        assertEquals(expected.toSet(), table.reads.toSet())
        assertEquals(window.cellCount, table.reads.size)
    }

    @Test
    fun aWindowFarFromSiteZeroReadsNothingBeforeIt() {
        // siteBlocks would have to walk every block from site 0 to reach this offset, and allocate and fill each one.
        val table = RecordingGenotypeTable(testTable(numberOfTaxa = 50, numberOfSites = 50_000))
        val window = GridWindow(firstTaxon = 0, taxonCount = 40, firstSite = 49_000, siteCount = 60)

        table.readWindow(window)

        assertEquals(window.cellCount, table.reads.size)
        assertTrue(table.reads.all { (_, site) -> site >= 49_000 }, "read a site outside the window")
    }

    @Test
    fun theCostOfAFrameDoesNotGrowWithTheTable() {
        val window = GridWindow(firstTaxon = 0, taxonCount = 30, firstSite = 0, siteCount = 60)
        val small = RecordingGenotypeTable(testTable(numberOfTaxa = 50, numberOfSites = 100))
        val large = RecordingGenotypeTable(testTable(numberOfTaxa = 50, numberOfSites = 50_000))

        small.readWindow(window)
        large.readWindow(window)

        assertEquals(small.reads.size, large.reads.size)
        assertEquals(window.cellCount, large.reads.size)
    }

    @Test
    fun theWindowHoldsTheCallsOfTheCellsItCovers() {
        val table = testTable(numberOfTaxa = 20, numberOfSites = 60)
        val window = GridWindow(firstTaxon = 5, taxonCount = 4, firstSite = 30, siteCount = 7)

        val genotypes = table.readWindow(window)

        for (taxon in window.taxa) {
            for (site in window.sites) {
                assertEquals(cellGenotype(taxon, site), genotypes.atTable(taxon, site), "cell ($taxon, $site)")
                assertEquals(
                    table.genotype(taxon, site),
                    genotypes[taxon - window.firstTaxon, site - window.firstSite],
                )
            }
        }
    }

    @Test
    fun readingOutsideTheWindowFails() {
        val genotypes = testTable(numberOfTaxa = 20, numberOfSites = 60)
            .readWindow(GridWindow(firstTaxon = 5, taxonCount = 4, firstSite = 30, siteCount = 7))

        assertFailsWith<IndexOutOfBoundsException> { genotypes[4, 0] }
        assertFailsWith<IndexOutOfBoundsException> { genotypes[0, 7] }
        assertFailsWith<IndexOutOfBoundsException> { genotypes[-1, 0] }
    }

    @Test
    fun aWindowOutsideTheTableFails() {
        val table = testTable(numberOfTaxa = 10, numberOfSites = 10)
        assertFailsWith<IllegalArgumentException> {
            table.readWindow(GridWindow(firstTaxon = 8, taxonCount = 5, firstSite = 0, siteCount = 1))
        }
        assertFailsWith<IllegalArgumentException> {
            table.readWindow(GridWindow(firstTaxon = 0, taxonCount = 1, firstSite = 9, siteCount = 2))
        }
    }

    @Test
    fun anEmptyWindowReadsNothing() {
        val table = RecordingGenotypeTable(testTable(numberOfTaxa = 10, numberOfSites = 10))

        val genotypes = table.readWindow(GridWindow.Empty)

        assertEquals(0, table.reads.size)
        assertEquals(GridWindow.Empty, genotypes.window)
    }

    @Test
    fun equalCallsDownAColumnBecomeOneRun() {
        // Every taxon is A/A here, so the column the grid paints is one rectangle rather than four.
        val call = cellGenotype(0, 0)
        val uniform = InMemoryGenotypeTable.of(
            taxa = TaxaList.of(List(4) { "taxon$it" }),
            positions = PositionList(List(1) { Position(Chromosome.of("1"), it + 1) }),
        ) { _, _ -> call }
        val window = GridWindow(firstTaxon = 0, taxonCount = 4, firstSite = 0, siteCount = 1)

        val runs = uniform.readWindow(window).columnRuns(column = 0)

        assertEquals(1, runs.size)
        assertEquals(ColumnRun(column = 0, firstRow = 0, rowCount = 4, genotype = call), runs.single())
    }

    @Test
    fun differingCallsDownAColumnStaySeparateRuns() {
        val table = testTable(numberOfTaxa = 5, numberOfSites = 1)
        val window = GridWindow(firstTaxon = 0, taxonCount = 5, firstSite = 0, siteCount = 1)

        val runs = table.readWindow(window).columnRuns(column = 0)

        // testTable rotates through A, C, G, T down a column, so taxa 0 and 4 share a call but are not adjacent.
        assertEquals(listOf(1, 1, 1, 1, 1), runs.map { it.rowCount })
        assertEquals(List(5) { cellGenotype(it, 0) }, runs.map { it.genotype })
        assertEquals(5, runs.sumOf { it.rowCount })
    }

    /** A [GenotypeTable] that remembers every cell read through it, in order. */
    private class RecordingGenotypeTable(private val table: GenotypeTable) : GenotypeTable {
        val reads = mutableListOf<Pair<Int, Int>>()

        override val taxa: TaxaList get() = table.taxa
        override val positions: PositionList get() = table.positions

        override fun genotype(taxon: Int, site: Int): Byte {
            reads.add(taxon to site)
            return table.genotype(taxon, site)
        }
    }
}
