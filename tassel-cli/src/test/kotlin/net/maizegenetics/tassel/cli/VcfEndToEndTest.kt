package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import net.maizegenetics.tassel.io.MISSING_BGZF_EOF
import java.io.ByteArrayOutputStream
import java.io.PrintStream
import kotlin.io.path.exists
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * `--input` on the tutorial VCF, against TASSEL 5's HapMap export of the same file.
 *
 * `filter-sites` with its defaults is the identity, so this exercises the whole path the CLI takes:
 * extension dispatch, `Vcf.read`, and `HapMap.write`.
 */
class VcfEndToEndTest {
    private val tutorial = testDataPath("tutorial/mdp_genotype_w_depth.vcf")
    private val tutorialGz = testDataPath("htslib/tutorial/mdp_genotype_w_depth.vcf.gz")

    @Test
    fun writesTheTassel5HapMapExportOfTheTutorialVcf() {
        assertTassel5HapMapExport(tutorial)
    }

    @Test
    fun writesTheTassel5HapMapExportOfTheTutorialVcfGz() {
        assertTassel5HapMapExport(tutorialGz)
    }

    @Test
    fun summarizesAVcfInput() {
        val directory = tempDirectory("tassel-vcf-summary")
        val base = directory.resolve("run")
        val result = tasselCommand().test("--quiet summary --input $tutorial --output $base")

        assertEquals(0, result.statusCode, result.output)
        val overview = outputText(directory.resolve("run.overview.txt"))
        assertContains(overview, "Number of Taxa\t281")
        assertContains(overview, "Number of Sites\t3093")
    }

    @Test
    fun summarizesABgzippedVcfLikeThePlainVcf() {
        val plainDir = tempDirectory("tassel-vcf-summary-plain")
        val gzDir = tempDirectory("tassel-vcf-summary-gz")
        val plain = tasselCommand().test("--quiet summary --input $tutorial --output ${plainDir.resolve("run")}")
        val gz = tasselCommand().test("--quiet summary --input $tutorialGz --output ${gzDir.resolve("run")}")

        assertEquals(0, plain.statusCode, plain.output)
        assertEquals(0, gz.statusCode, gz.output)
        listOf("run.overview.txt", "run.alleles.txt", "run.sites.txt", "run.taxa.txt").forEach { name ->
            assertEquals(outputText(plainDir.resolve(name)), outputText(gzDir.resolve(name)), name)
        }
    }

    @Test
    fun readsPlainGzipVcfOnTheSequentialPath() {
        val directory = tempDirectory("tassel-vcf-plain-gzip")
        val input = testDataPath("htslib/variants/plain.vcf.gz")
        val result = tasselCommand().test("--quiet summary --input $input --output ${directory.resolve("run")}")
        assertEquals(0, result.statusCode, result.output)
        val overview = outputText(directory.resolve("run.overview.txt"))
        assertContains(overview, "Number of Taxa\t2")
        assertContains(overview, "Number of Sites\t2")
    }

    @Test
    fun warnsWhenABgzippedVcfHasNoEofBlock() {
        val directory = tempDirectory("tassel-vcf-no-eof")
        val input = testDataPath("htslib/variants/no-eof.vcf.gz")
        val stderr = captureStderr {
            val result = tasselCommand().test(
                "--quiet summary --input $input --output ${directory.resolve("run")}",
            )
            assertEquals(0, result.statusCode, result.output)
        }
        assertContains(stderr, MISSING_BGZF_EOF)
    }

    private fun assertTassel5HapMapExport(input: String) {
        val directory = tempDirectory("tassel-vcf")
        val base = directory.resolve("mdp_genotype_w_depth")
        val result = tasselCommand().test("--quiet filter-sites --input $input --output $base")

        assertEquals(0, result.statusCode, result.output)
        val written = directory.resolve("mdp_genotype_w_depth.hmp.txt")
        assertTrue(written.exists(), "filter-sites wrote no mdp_genotype_w_depth.hmp.txt")

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
        assertEquals(golden.joinToString("\n"), outputText(written))
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
