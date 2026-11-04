package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.core.CliktError
import com.github.ajalt.clikt.testing.test
import net.maizegenetics.tassel.analysis.genotypeSummary
import net.maizegenetics.tassel.distance.centeredIbs
import kotlin.io.path.exists
import kotlin.io.path.readText
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class CliFilesTest {
    private val table = simulateGenotypes(taxa = 6, sites = 12, seed = 5)

    @Test
    fun theExtensionIsAppendedUnlessTheBaseAlreadyHasIt() {
        assertEquals("kinship.txt", outputPath("kinship", null, TEXT_EXTENSION))
        assertEquals("kinship.txt", outputPath("kinship.txt", null, TEXT_EXTENSION))
        assertEquals("out/filtered.hmp.txt", outputPath("out/filtered", null, HAPMAP_EXTENSION))
        assertEquals("out/filtered.hmp.txt", outputPath("out/filtered.hmp.txt", null, HAPMAP_EXTENSION))
    }

    @Test
    fun aResultWithSeveralPartsNamesEachOne() {
        assertEquals("run.sites.txt", outputPath("run", "sites", TEXT_EXTENSION))
        assertEquals("run.sites.txt", outputPath("run.txt", "sites", TEXT_EXTENSION))
    }

    @Test
    fun eachResultTypeKnowsWhatFilesItBecomes() {
        assertEquals(listOf(null to HAPMAP_EXTENSION), outputParts(table).map { it.part to it.extension })
        assertEquals(
            listOf(null to TEXT_EXTENSION),
            outputParts(centeredIbs(table)).map { it.part to it.extension },
        )
        assertEquals(
            listOf("overview", "alleles", "sites", "taxa"),
            outputParts(genotypeSummary(table)).map { it.part },
        )
    }

    @Test
    fun aResultNothingCanWriteIsAnInternalError() {
        val failure = assertFailsWith<CliktError> { outputParts("not a result") }
        assertContains(failure.message!!, "no writer for a String result")
        assertEquals(EXIT_INTERNAL_ERROR, failure.statusCode)
    }

    @Test
    fun onlyTheRequestedSummaryTablesAreWritten() {
        val directory = tempDirectory("tassel-summary-parts")
        val input = directory.resolve("in.hmp.txt")
        writeHapMap(table, input)
        val result = tasselCommand().test(
            "--quiet summary --input $input --no-site-summary --output ${directory.resolve("run")}",
        )

        assertEquals(0, result.statusCode, result.output)
        assertTrue(directory.resolve("run.overview.txt").exists())
        assertTrue(directory.resolve("run.alleles.txt").exists())
        assertTrue(directory.resolve("run.taxa.txt").exists())
        assertTrue(!directory.resolve("run.sites.txt").exists(), "--no-site-summary still wrote the site table")
    }

    @Test
    fun aMissingOutputDirectoryIsCreated() {
        val directory = tempDirectory("tassel-mkdir")
        val input = directory.resolve("in.hmp.txt")
        writeHapMap(table, input)
        val output = directory.resolve("new").resolve("deeper").resolve("filtered")
        val result = tasselCommand().test("--quiet filter-taxa --input $input --output $output")

        assertEquals(0, result.statusCode, result.output)
        val written = directory.resolve("new").resolve("deeper").resolve("filtered.hmp.txt")
        assertTrue(written.exists())
        assertContains(written.readText(), "taxon_1")
    }
}
