package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlin.io.path.createFile
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

    @Test
    fun writesTheTassel5HapMapExportOfTheTutorialVcf() {
        val directory = tempDirectory("tassel-vcf")
        val base = directory.resolve("mdp_genotype_w_depth")
        val result = tasselCommand().test("--quiet filter-sites --input $tutorial --output $base")

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
    fun rejectsABgzipCompressedVcf() {
        val directory = tempDirectory("tassel-vcf-gz")
        val input = directory.resolve("genotypes.vcf.gz").createFile()
        val result = tasselCommand().test("--quiet summary --input $input --output ${directory.resolve("run")}")

        assertEquals(EXIT_USER_ERROR, result.statusCode, result.output)
        assertContains(result.output, "bgzip")
    }
}
