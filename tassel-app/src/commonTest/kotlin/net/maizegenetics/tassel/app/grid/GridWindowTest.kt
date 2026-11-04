package net.maizegenetics.tassel.app.grid

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/** The mapping from a scroll offset and a viewport size to the block of table cells the grid shows. */
class GridWindowTest {
    @Test
    fun unscrolledWindowStartsAtTheFirstCell() {
        val window = window(taxonOffset = 0f, siteOffset = 0f)
        assertEquals(GridWindow(firstTaxon = 0, taxonCount = 10, firstSite = 0, siteCount = 10), window)
        assertEquals(0..9, window.taxa)
        assertEquals(0..9, window.sites)
        assertEquals(100, window.cellCount)
    }

    @Test
    fun scrollingByWholeCellsMovesTheWindowWithoutGrowingIt() {
        val window = window(taxonOffset = 14f * 3, siteOffset = 10f * 20)
        assertEquals(GridWindow(firstTaxon = 3, taxonCount = 10, firstSite = 20, siteCount = 10), window)
    }

    @Test
    fun partlyVisibleCellsAtBothEdgesAreIncluded() {
        // Half a cell of scroll leaves a sliver of the cell above and a sliver of the one below the viewport.
        val window = window(taxonOffset = 7f, siteOffset = 5f)
        assertEquals(GridWindow(firstTaxon = 0, taxonCount = 11, firstSite = 0, siteCount = 11), window)
    }

    @Test
    fun theWindowStopsAtTheEndOfTheTable() {
        val window = window(taxonOffset = 14f * 95, siteOffset = 10f * 995)
        assertEquals(GridWindow(firstTaxon = 95, taxonCount = 5, firstSite = 995, siteCount = 5), window)
        assertTrue(window.taxa.last < TAXA)
        assertTrue(window.sites.last < SITES)
    }

    @Test
    fun overshootingTheEndClampsToTheLastCell() {
        // A fling can hand us an offset past the end; the window has to stay inside the table rather than throw.
        val window = window(taxonOffset = 1e6f, siteOffset = 1e6f)
        assertEquals(GridWindow(firstTaxon = TAXA - 1, taxonCount = 1, firstSite = SITES - 1, siteCount = 1), window)
    }

    @Test
    fun aWindowSizeIsSetByTheViewportAndNotByTheTable() {
        val small = gridWindow(TAXA, SITES, 0f, 0f, WIDTH, HEIGHT, CELL_WIDTH, CELL_HEIGHT)
        val huge = gridWindow(1_000_000, 1_000_000, 0f, 0f, WIDTH, HEIGHT, CELL_WIDTH, CELL_HEIGHT)
        assertEquals(small.cellCount, huge.cellCount)
    }

    @Test
    fun anEmptyTableOrAnUnmeasuredViewportGivesAnEmptyWindow() {
        assertEquals(GridWindow.Empty, gridWindow(0, SITES, 0f, 0f, WIDTH, HEIGHT, CELL_WIDTH, CELL_HEIGHT))
        assertEquals(GridWindow.Empty, gridWindow(TAXA, 0, 0f, 0f, WIDTH, HEIGHT, CELL_WIDTH, CELL_HEIGHT))
        assertEquals(GridWindow.Empty, gridWindow(TAXA, SITES, 0f, 0f, 0f, HEIGHT, CELL_WIDTH, CELL_HEIGHT))
        assertEquals(GridWindow.Empty, gridWindow(TAXA, SITES, 0f, 0f, WIDTH, -1f, CELL_WIDTH, CELL_HEIGHT))
        assertEquals(0, GridWindow.Empty.cellCount)
    }

    @Test
    fun cellsMustHaveAPositiveSize() {
        assertFailsWith<IllegalArgumentException> {
            gridWindow(TAXA, SITES, 0f, 0f, WIDTH, HEIGHT, cellWidth = 0f, cellHeight = CELL_HEIGHT)
        }
    }

    private fun window(taxonOffset: Float, siteOffset: Float): GridWindow =
        gridWindow(TAXA, SITES, taxonOffset, siteOffset, WIDTH, HEIGHT, CELL_WIDTH, CELL_HEIGHT)

    private companion object {
        const val TAXA = 100
        const val SITES = 1000
        const val CELL_WIDTH = 10f
        const val CELL_HEIGHT = 14f

        /** Exactly ten cells wide and ten tall, so an unscrolled window is 10 by 10. */
        const val WIDTH = 10 * CELL_WIDTH
        const val HEIGHT = 10 * CELL_HEIGHT
    }
}
