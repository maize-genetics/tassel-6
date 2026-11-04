package net.maizegenetics.tassel.app.data

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.app.grid.GridWindow
import net.maizegenetics.tassel.app.grid.readWindow
import net.maizegenetics.tassel.app.theme.LightAlleleColors
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * The desktop half of [GenotypeFiles], as far as it goes without a window: the tutorial dataset really is on the
 * classpath where `bundleTutorialData` puts it, and it really parses.
 *
 * [DesktopGenotypeFiles.choose] opens a native `FileDialog` and is not covered here; the maintainer exercises it.
 */
class DesktopGenotypeFilesTest {
    @Test
    fun theBundledTutorialDatasetLoads() = runTest {
        val table = DesktopGenotypeFiles().tutorial().readGenotypes()

        assertEquals(281, table.numberOfTaxa)
        assertEquals(3093, table.numberOfSites)
        assertEquals(TUTORIAL_FILE_NAME, DesktopGenotypeFiles().tutorial().name)
    }

    @Test
    fun aViewportOfTheTutorialDatasetReadsAndColors() = runTest {
        // The end of the loading path the app walks: bundled bytes, a window, a color per cell. 1,200 cells out of
        // 281 x 3,093 = 869,133, which is why the grid can draw a frame of this dataset at all.
        val table = DesktopGenotypeFiles().tutorial().readGenotypes()
        val window = GridWindow(firstTaxon = 0, taxonCount = 40, firstSite = 3_000, siteCount = 30)

        val genotypes = table.readWindow(window)
        val colors = window.taxa.flatMap { taxon ->
            window.sites.map { site -> LightAlleleColors.of(genotypes.atTable(taxon, site)) }
        }

        assertEquals(window.cellCount, colors.size)
        assertEquals(1_200, window.cellCount)
        assertTrue(colors.toSet().size > 1, "a window near the end of the dataset should not be one flat color")
    }
}
