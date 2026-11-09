package net.maizegenetics.tassel.app.data

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.runBlocking
import kotlinx.coroutines.test.runTest
import kotlinx.coroutines.withContext
import kotlinx.coroutines.withTimeout
import net.maizegenetics.tassel.app.grid.GridWindow
import net.maizegenetics.tassel.app.grid.readWindow
import net.maizegenetics.tassel.app.theme.LightAlleleColors
import net.maizegenetics.tassel.io.IndexedVcf
import net.maizegenetics.tassel.io.openRandomAccessSource
import net.maizegenetics.tassel.io.parseVcfRegion
import net.maizegenetics.tassel.ops.ProgressReporter
import java.nio.file.Files
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertNull
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
    fun siblingIndexPrefersCsiThenTbi() {
        val dir = Files.createTempDirectory("tassel-index")
        try {
            val vcf = dir.resolve("sample.vcf.gz")
            Files.write(vcf, byteArrayOf(1))
            assertNull(siblingIndex(vcf.toString()))

            val tbi = dir.resolve("sample.vcf.gz.tbi")
            Files.write(tbi, byteArrayOf(1))
            assertEquals(tbi.toString(), siblingIndex(vcf.toString()))

            val csi = dir.resolve("sample.vcf.gz.csi")
            Files.write(csi, byteArrayOf(1))
            assertEquals(csi.toString(), siblingIndex(vcf.toString()))
        } finally {
            java.io.File(dir.toString()).deleteRecursively()
        }
    }

    @Test
    fun aWheatRegionPastTwoToThe29MatchesIndexedVcf() = runTest {
        val root = System.getenv("TASSEL_TESTDATA_DIR")?.trimEnd('/')
            ?: error("TASSEL_TESTDATA_DIR is not set; run tests through Gradle so it points at testdata/")
        val dataPath = "$root/htslib/wheat/wheat.vcf.gz"
        val indexPath = "$dataPath.csi"
        val region = parseVcfRegion("chr3B:599000000-601000000")
        val expected = openRandomAccessSource(dataPath).use { data ->
            openRandomAccessSource(indexPath).use { index ->
                IndexedVcf.open(data, index).use { indexed -> indexed.read(region).table }
            }
        }

        val actual = openRegionSource(pathGenotypeFile(dataPath), pathGenotypeFile(indexPath)).use { session ->
            assertEquals(851_934_019L, session.contigs.first { it.name == "chr3B" }.length)
            assertNull(session.contigs.find { it.name == "chr6D" })
            session.read(region, ProgressReporter.None)
        }

        assertTrue(actual.numberOfSites > 0)
        assertEquals(expected.taxa.names, actual.taxa.names)
        assertEquals(expected.numberOfSites, actual.numberOfSites)
        for (site in 0 until expected.numberOfSites) {
            assertEquals(expected.positions[site], actual.positions[site], "site $site")
            assertEquals(expected.positions[site].knownVariants, actual.positions[site].knownVariants, "site $site")
            assertContentEquals(expected.siteGenotypes(site), actual.siteGenotypes(site), "site $site")
        }
    }

    @Test
    fun aRegionReadFinishesOnTheBackgroundDispatcher() = runBlocking {
        // The BGZF readahead window is a child of the caller's job. withContext waits for children, so the reader
        // has to be closed before that call returns. This is the shape AppViewModel uses.
        val root = System.getenv("TASSEL_TESTDATA_DIR")?.trimEnd('/')
            ?: error("TASSEL_TESTDATA_DIR is not set; run tests through Gradle so it points at testdata/")
        val dataPath = "$root/htslib/wheat/wheat.vcf.gz"
        val session = openRegionSource(pathGenotypeFile(dataPath), pathGenotypeFile("$dataPath.csi"))
        val table = withTimeout(15_000) {
            withContext(Dispatchers.Default) {
                try {
                    session.read(parseVcfRegion("chr3B:599000000-601000000"), ProgressReporter.None)
                } finally {
                    session.close()
                }
            }
        }
        assertTrue(table.numberOfSites > 0)
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
