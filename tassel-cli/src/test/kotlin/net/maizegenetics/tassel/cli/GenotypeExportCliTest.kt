package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlinx.coroutines.runBlocking
import kotlinx.io.buffered
import kotlinx.io.files.SystemFileSystem
import kotlinx.io.readString
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.Bcf
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.io.decompressing
import java.nio.file.Path
import kotlin.io.path.exists
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * `--output` on the tutorial HapMap, against TASSEL 5's VCF export of the same file (ADR-0022).
 *
 * `filter-sites` with its defaults returns the input table, so the suffix is the only thing this exercises.
 * `.vcf.gz` is BGZF, so the check is the inflated text, not TASSEL 5's `GZIPOutputStream` bytes.
 */
class GenotypeExportCliTest {
    private val tutorial = testDataPath("tutorial/mdp_genotype.hmp.txt")

    @Test
    fun writesTheTassel5VcfExportOfTheTutorialHapMap() {
        assertEquals(testDataText(GOLDEN), outputText(filterSites("mdp_genotype.vcf")))
    }

    @Test
    fun aBgzippedVcfInflatesToTheTassel5VcfExport() {
        assertEquals(testDataText(GOLDEN), normalizeLineEndings(inflate(filterSites("mdp_genotype.vcf.gz"))))
    }

    @Test
    fun aBcfReadsBackEqualToTheInputTable() {
        val written = readBcf(filterSites("mdp_genotype.bcf"))
        assertSameGenotypeTable(readHapMap(tutorial), written)
    }

    @Test
    fun helpStatesTheSuffixRule() {
        val root = tasselCommand().test("--help")
        assertEquals(0, root.statusCode, root.output)
        assertContains(root.output, ".hmp.txt")
        assertContains(root.output, ".vcf,")
        assertContains(root.output, ".vcf.gz")
        assertContains(root.output, ".bcf")

        val sites = tasselCommand().test("filter-sites --help")
        assertEquals(0, sites.statusCode, sites.output)
        assertContains(sites.output, "filtered.hmp.txt")
        assertContains(sites.output, "filtered.vcf,")
        assertContains(sites.output, "filtered.vcf.gz")
        assertContains(sites.output, "filtered.bcf")
    }

    private fun filterSites(fileName: String): Path {
        val directory = tempDirectory("tassel-genotype-export")
        val output = directory.resolve(fileName)
        val result = tasselCommand().test("--quiet filter-sites --input $tutorial --output $output")
        assertEquals(0, result.statusCode, result.output)
        assertTrue(output.exists(), "filter-sites wrote no $fileName")
        return output
    }

    private fun readHapMap(path: String): GenotypeTable = runBlocking {
        SystemFileSystem.source(kotlinx.io.files.Path(path)).buffered().use { HapMap.read(it) }
    }

    private fun readBcf(path: Path): GenotypeTable = runBlocking {
        SystemFileSystem.source(kotlinx.io.files.Path(path.toString())).buffered().use { Bcf.read(it).table }
    }

    private fun inflate(path: Path): String = runBlocking {
        SystemFileSystem.source(kotlinx.io.files.Path(path.toString())).buffered().use { raw ->
            val opened = raw.decompressing()
            try {
                opened.source.readString()
            } finally {
                if (opened.source !== raw) opened.source.close()
            }
        }
    }

    private companion object {
        const val GOLDEN = "golden/vcf-export/mdp_genotype.vcf"
    }
}
