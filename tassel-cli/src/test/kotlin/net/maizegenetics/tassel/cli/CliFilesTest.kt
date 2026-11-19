package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.core.CliktError
import com.github.ajalt.clikt.testing.test
import kotlinx.coroutines.runBlocking
import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import kotlinx.io.readString
import net.maizegenetics.tassel.analysis.genotypeSummary
import net.maizegenetics.tassel.distance.centeredIbs
import net.maizegenetics.tassel.io.Bcf
import net.maizegenetics.tassel.io.decompressing
import java.io.ByteArrayOutputStream
import java.io.PrintStream
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
    fun theExtensionIsAppendedUnlessTheBaseAlreadyHasIt() = runBlocking {
        assertEquals("kinship.txt", outputPath("kinship", null, TEXT_EXTENSION))
        assertEquals("kinship.txt", outputPath("kinship.txt", null, TEXT_EXTENSION))
        assertEquals("out/filtered.hmp.txt", outputPath("out/filtered", null, HAPMAP_EXTENSION))
        assertEquals("out/filtered.hmp.txt", outputPath("out/filtered.hmp.txt", null, HAPMAP_EXTENSION))
    }

    @Test
    fun aResultWithSeveralPartsNamesEachOne() = runBlocking {
        assertEquals("run.sites.txt", outputPath("run", "sites", TEXT_EXTENSION))
        assertEquals("run.sites.txt", outputPath("run.txt", "sites", TEXT_EXTENSION))
        // `.txt` stays case-sensitive. Genotype suffixes do not.
        assertEquals("run.TXT.sites.txt", outputPath("run.TXT", "sites", TEXT_EXTENSION))
    }

    @Test
    fun aGenotypeSuffixKeepsTheCallersSpelling() = runBlocking {
        assertEquals("out/filtered.VCF", outputPath("out/filtered.VCF", null, VCF_EXTENSION))
        assertEquals("out/filtered.VCF.GZ", outputPath("out/filtered.VCF.GZ", null, VCF_GZIP_EXTENSION))
        assertEquals("out/filtered.BCF", outputPath("out/filtered.BCF", null, BCF_EXTENSION))
        assertEquals("out/filtered.HMP.TXT", outputPath("out/filtered.HMP.TXT", null, HAPMAP_EXTENSION))
    }

    @Test
    fun eachResultTypeKnowsWhatFilesItBecomes() = runBlocking {
        assertEquals(listOf(null to HAPMAP_EXTENSION), outputParts(table).map { it.part to it.extension })
        assertEquals(listOf(null to HAPMAP_EXTENSION), outputParts(table, null).map { it.part to it.extension })
        assertEquals(listOf(null to HAPMAP_EXTENSION), outputParts(table, "out").map { it.part to it.extension })
        assertEquals(
            listOf(null to HAPMAP_EXTENSION),
            outputParts(table, "out.HMP.TXT").map { it.part to it.extension },
        )
        assertEquals(
            listOf(null to VCF_EXTENSION),
            outputParts(table, "out/filtered.vcf").map { it.part to it.extension },
        )
        assertEquals(
            listOf(null to VCF_GZIP_EXTENSION),
            outputParts(table, "out/filtered.VCF.GZ").map { it.part to it.extension },
        )
        assertEquals(
            listOf(null to BCF_EXTENSION),
            outputParts(table, "out/filtered.bcf").map { it.part to it.extension },
        )
        assertEquals(
            listOf(null to TEXT_EXTENSION),
            outputParts(centeredIbs(table)).map { it.part to it.extension },
        )
        assertEquals(
            listOf(null to TEXT_EXTENSION),
            outputParts(centeredIbs(table), "kinship.vcf").map { it.part to it.extension },
        )
        assertEquals(
            listOf("overview", "alleles", "sites", "taxa"),
            outputParts(genotypeSummary(table)).map { it.part },
        )
        assertEquals(
            listOf(TEXT_EXTENSION, TEXT_EXTENSION, TEXT_EXTENSION, TEXT_EXTENSION),
            outputParts(genotypeSummary(table), "run.vcf").map { it.extension },
        )
    }

    @Test
    fun theSelectedGenotypeWriterIsTheOneTheSuffixNames() = runBlocking {
        val vcf = Buffer().also { outputParts(table, "x.vcf").single().write(it) }.readString()
        assertTrue(vcf.startsWith("##fileformat=VCFv4.0\n"))

        val gzipped = Buffer().also { outputParts(table, "x.vcf.gz").single().write(it) }.readByteArray()
        val inflated = inflate(gzipped)
        assertEquals(vcf, inflated)
        assertTrue(gzipped.size < vcf.encodeToByteArray().size, "stored blocks would be larger than the text")

        val bcf = Buffer().also { outputParts(table, "x.bcf").single().write(it) }.readByteArray()
        assertSameGenotypeTable(table, Bcf.read(Buffer().also { it.write(bcf) }).table)
    }

    @Test
    fun standardOutputStaysHapMap() = runBlocking {
        val previous = System.out
        val buffer = ByteArrayOutputStream()
        System.setOut(PrintStream(buffer, true, Charsets.UTF_8))
        try {
            writeResult(table, null)
        } finally {
            System.setOut(previous)
        }
        assertTrue(buffer.toString(Charsets.UTF_8).startsWith("rs#\t"))
    }

    @Test
    fun aMixedCaseVcfSuffixWritesThatFile() = runBlocking {
        val directory = tempDirectory("tassel-vcf-case")
        val input = directory.resolve("in.hmp.txt")
        writeHapMap(table, input)
        val output = directory.resolve("filtered.VCF")
        val result = tasselCommand().test("--quiet filter-taxa --input $input --output $output")

        assertEquals(0, result.statusCode, result.output)
        assertTrue(output.exists())
        assertTrue(!directory.resolve("filtered.VCF.vcf").exists())
        assertTrue(!directory.resolve("filtered.hmp.txt").exists())
        assertTrue(output.readText().startsWith("##fileformat=VCFv4.0"))
    }

    @Test
    fun aResultNothingCanWriteIsAnInternalError() = runBlocking {
        val failure = assertFailsWith<CliktError> { outputParts("not a result") }
        assertContains(failure.message!!, "no writer for a String result")
        assertEquals(EXIT_INTERNAL_ERROR, failure.statusCode)
    }

    @Test
    fun onlyTheRequestedSummaryTablesAreWritten() = runBlocking {
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
    fun aMissingOutputDirectoryIsCreated() = runBlocking {
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

    private fun inflate(bytes: ByteArray): String {
        val raw = Buffer().also { it.write(bytes) }
        val opened = raw.decompressing()
        return try {
            opened.source.readString()
        } finally {
            if (opened.source !== raw) opened.source.close()
        }
    }
}
