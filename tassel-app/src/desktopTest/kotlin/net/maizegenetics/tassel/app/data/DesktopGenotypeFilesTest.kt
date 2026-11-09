package net.maizegenetics.tassel.app.data

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.app.grid.GridWindow
import net.maizegenetics.tassel.app.grid.readWindow
import net.maizegenetics.tassel.app.theme.LightAlleleColors
import kotlin.test.Test
import kotlin.test.assertContentEquals
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

    @Test
    fun aBgzippedVcfReadsTheSameTableAsThePlainVcf() = runTest {
        val root = System.getenv("TASSEL_TESTDATA_DIR")?.trimEnd('/')
            ?: error("TASSEL_TESTDATA_DIR is not set; run tests through Gradle so it points at testdata/")
        val plain = pathGenotypeFile("$root/tutorial/mdp_genotype_w_depth.vcf").readGenotypes()
        val gzipped = pathGenotypeFile("$root/htslib/tutorial/mdp_genotype_w_depth.vcf.gz").readGenotypes()

        assertEquals(281, plain.numberOfTaxa)
        assertEquals(3_093, plain.numberOfSites)
        assertEquals(plain.taxa.names, gzipped.taxa.names)
        assertEquals(plain.numberOfSites, gzipped.numberOfSites)
        for (site in 0 until plain.numberOfSites) {
            assertContentEquals(plain.siteGenotypes(site), gzipped.siteGenotypes(site), "site $site")
        }
    }
}
