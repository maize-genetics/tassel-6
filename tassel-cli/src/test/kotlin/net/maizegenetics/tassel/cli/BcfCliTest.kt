package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.tassel.io.MISSING_BGZF_EOF
import java.io.ByteArrayOutputStream
import java.io.PrintStream
import java.nio.file.Files
import java.nio.file.Path
import kotlin.io.path.exists
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

/**
 * `--input` and `--region` on a BCF, through the same generated `filter-sites` command as
 * [IndexedVcfCliTest]. Defaults keep every site, so the HapMap it writes is the read.
 */
class BcfCliTest {
    private val tutorialBcf = testDataPath("htslib/bcf/tutorial/mdp_genotype_w_depth.bcf")
    private val wheatBcf = testDataPath("htslib/bcf/wheat/wheat.bcf")
    private val wheatGz = testDataPath("htslib/wheat/wheat.vcf.gz")
    private val tutorialTbi = testDataPath("htslib/tutorial/mdp_genotype_w_depth.vcf.gz.tbi")

    @Test
    fun filterSitesWritesTheTassel5HapMapExportOfTheTutorialBcf() {
        assertTassel5HapMapExport(tutorialBcf)
    }

    @Test
    fun aBgzippedBcfWithAnUnusualNameIsStillBcf() {
        val directory = tempDirectory("tassel-bcf-alias")
        val alias = directory.resolve("notes.bin")
        Files.copy(Path.of(tutorialBcf), alias)
        assertTassel5HapMapExport(alias.toString())
    }

    @Test
    fun wheatRegionMatchesTheSameRegionOfTheVcfGz() {
        val window = "chr3B:599000000-601000000"
        assertEquals(filterSitesHapMap(wheatGz, window), filterSitesHapMap(wheatBcf, window), window)
        assertEquals(filterSitesHapMap(wheatGz, "chr3B"), filterSitesHapMap(wheatBcf, "chr3B"), "chr3B")
    }

    @Test
    fun aRawBcfMatchesTheDecodedVcf() {
        val raw = testDataPath("htslib/bcf/small/small.raw.bcf")
        val view = testDataPath("htslib/bcf/small/small.view.vcf")
        assertEquals(filterSitesHapMap(view), filterSitesHapMap(raw))
    }

    @Test
    fun exitCodesMatchTheIndexedVcfCli() {
        assertUserError("filter-sites --input no/such/file.bcf", "no such file")
        assertUserError("filter-sites --input no/such/file.bcf --region 1:1-1000", "no such file")

        val planted = plantedWithTbi()
        val missing = tasselCommand().test(
            "filter-sites --input $planted --region 1:157104-157104",
        )
        assertUserError(
            missing.statusCode,
            missing.output,
            "no index found for $planted",
            "bcftools index -c $planted",
        )
        assertFalse(missing.output.contains("tabix"), missing.output)

        assertUserError(
            "filter-sites --input $tutorialBcf --region 99:1-1000",
            "99",
            "not in the index",
        )
        assertUserError(
            "filter-sites --input $tutorialBcf --region 1:10068726-157104",
            "1:10068726-157104",
            "inverted",
        )
        assertUserError(
            "filter-sites --input $tutorialBcf --region 1:157104-157104 --index $tutorialTbi",
            "cannot read index",
            "TBI index already contains contig names",
        )
        val raw = testDataPath("htslib/bcf/small/small.raw.bcf")
        assertUserError(
            "filter-sites --input $raw --region 1:1-100",
            raw,
            "not bgzip compressed",
        )
    }

    @Test
    fun aBcfPathThatIsHapMapDoesNotReadAsHapMap() {
        val directory = tempDirectory("tassel-bcf-hapmap")
        val bcf = directory.resolve("sample.bcf")
        Files.copy(Path.of(testDataPath("tutorial/mdp_genotype.hmp.txt")), bcf)
        val output = directory.resolve("out")
        val result = tasselCommand().test("filter-sites --input $bcf --output $output")
        assertUserError(result.statusCode, result.output, "Invalid BCF2 magic string")
        assertTrue(!directory.resolve("out.hmp.txt").exists(), result.output)
    }

    @Test
    fun aMissingEofBlockIsAWarning() {
        val directory = tempDirectory("tassel-bcf-no-eof")
        val input = testDataPath("htslib/bcf/corrupt/no-eof.bcf")
        val stderr = captureStderr {
            val result = tasselCommand().test(
                "--quiet filter-sites --input $input --output ${directory.resolve("out")}",
            )
            assertEquals(0, result.statusCode, result.output)
        }
        assertContains(stderr, MISSING_BGZF_EOF)
    }

    @Test
    fun helpMentionsBcf() {
        val root = tasselCommand().test("--help")
        assertEquals(0, root.statusCode, root.output)
        assertContains(root.output, "BCF")
        val sites = tasselCommand().test("filter-sites --help")
        assertEquals(0, sites.statusCode, sites.output)
        assertContains(sites.output, ".bcf")
        assertContains(sites.output, ".csi")
    }

    private fun assertTassel5HapMapExport(input: String) {
        val written = filterSitesHapMap(input)
        // Maintainer decision (tassel-6#14, 2026-11-04): compare from the rs# header row down. The
        // golden's first 281 lines are ##SAMPLE= taxa annotations that TASSEL 5 copies from the VCF
        // onto its Taxon objects and re-emits; Taxon is frozen as a name (#10) and HapMap.write
        // emits no ## lines, so they cannot round trip. See docs/porting-ledger.md and VcfTest.
        val goldenLines = testDataText("golden/vcf/mdp_genotype_w_depth.hmp.txt").lines()
        val skipped = goldenLines.takeWhile { it.startsWith("##") }
        assertEquals(281, skipped.size, "skipped ##SAMPLE= lines")
        assertTrue(skipped.all { it.startsWith("##SAMPLE=<ID=") }, "every skipped line is a ##SAMPLE= annotation")
        val golden = goldenLines.drop(skipped.size)
        assertTrue(golden.first().startsWith("rs#\talleles\tchrom\tpos\t"), "first compared line")
        assertEquals(golden.joinToString("\n"), written)
    }

    private fun filterSitesHapMap(input: String, region: String? = null): String {
        val directory = tempDirectory("tassel-bcf")
        val base = directory.resolve("out")
        val regionArg = region?.let { "--region $it" }.orEmpty()
        val result = tasselCommand().test(
            "--quiet filter-sites --input $input $regionArg --output $base",
        )
        assertEquals(0, result.statusCode, result.output)
        val written = directory.resolve("out.hmp.txt")
        assertTrue(written.exists(), "filter-sites wrote no out.hmp.txt for $input $regionArg")
        return outputText(written)
    }

    private fun plantedWithTbi(): Path {
        val directory = tempDirectory("tassel-bcf-tbi-only")
        val bcf = directory.resolve("sample.bcf")
        Files.copy(Path.of(tutorialBcf), bcf)
        Files.copy(Path.of(tutorialTbi), directory.resolve("sample.bcf.tbi"))
        return bcf
    }

    private fun assertUserError(args: String, vararg needles: String) {
        val result = tasselCommand().test(args)
        assertUserError(result.statusCode, result.output, *needles)
    }

    private fun assertUserError(status: Int, output: String, vararg needles: String) {
        assertEquals(EXIT_USER_ERROR, status, output)
        assertTrue(status != EXIT_INTERNAL_ERROR, output)
        for (needle in needles) assertContains(output, needle)
    }

    /** [CliLog.warn] writes to [System.err], which Clikt's test recorder does not capture. */
    private fun captureStderr(block: () -> Unit): String {
        val previous = System.err
        val buffer = ByteArrayOutputStream()
        System.setErr(PrintStream(buffer, true, Charsets.UTF_8))
        try {
            block()
            return buffer.toString(Charsets.UTF_8)
        } finally {
            System.setErr(previous)
        }
    }
}
