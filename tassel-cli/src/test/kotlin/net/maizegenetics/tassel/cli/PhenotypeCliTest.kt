package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlinx.coroutines.runBlocking
import net.maizegenetics.tassel.analysis.LogBase
import net.maizegenetics.tassel.analysis.PhenotypeTransform
import net.maizegenetics.tassel.analysis.summarizePhenotype
import net.maizegenetics.tassel.analysis.transformPhenotype
import net.maizegenetics.tassel.filter.PhenotypeFilter
import net.maizegenetics.tassel.filter.filterPhenotype
import net.maizegenetics.tassel.join.joinPhenotypes
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.table.ResultTable
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class PhenotypeCliTest {
    private val traits = testDataPath("tutorial/mdp_traits.txt")
    private val structure = testDataPath("tutorial/mdp_population_structure.txt")
    private val genotypes = testDataPath("tutorial/mdp_genotype.hmp.txt")

    private fun text(table: PhenotypeTable) = buildString { table.writeDelimited(this) }

    private fun text(table: ResultTable) = buildString { table.writeDelimited(this) }

    private fun run(args: String) = tasselCommand().test("--quiet $args")

    @Test
    fun filterPhenotypeWritesWhatTheKernelProduces() = runBlocking {
        val directory = tempDirectory("tassel-pheno-filter")
        val output = directory.resolve("filtered")
        val result = run("filter-phenotype --input $traits --keep-attributes EarHT,EarDia --output $output")
        assertEquals(0, result.statusCode, result.output)

        val filter = PhenotypeFilter(keepAttributes = listOf("EarHT", "EarDia"))
        val expected = filterPhenotype(readPhenotype(traits), filter)
        assertEquals(text(expected), outputText(directory.resolve("filtered.txt")))
    }

    @Test
    fun transformPhenotypeWritesWhatTheKernelProduces() = runBlocking {
        val directory = tempDirectory("tassel-pheno-transform")
        val output = directory.resolve("logged")
        val result = run("transform-phenotype --input $traits --log natural --traits EarHT --output $output")
        assertEquals(0, result.statusCode, result.output)

        val expected = transformPhenotype(
            readPhenotype(traits),
            PhenotypeTransform(traits = listOf("EarHT"), log = LogBase.natural),
        ).table
        assertEquals(text(expected), outputText(directory.resolve("logged.txt")))
    }

    @Test
    fun phenotypeSummaryWritesThreeParts() = runBlocking {
        val directory = tempDirectory("tassel-pheno-summary")
        val output = directory.resolve("summary")
        val result = run("phenotype-summary --input $traits --output $output")
        assertEquals(0, result.statusCode, result.output)

        val expected = summarizePhenotype(readPhenotype(traits))
        assertEquals(text(expected.attributes), outputText(directory.resolve("summary.attributes.txt")))
        assertEquals(text(expected.factorLevels), outputText(directory.resolve("summary.levels.txt")))
        assertEquals(text(expected.taxa), outputText(directory.resolve("summary.taxa.txt")))
    }

    @Test
    fun joinPhenotypesTakesARepeatedInputInOrder() = runBlocking {
        val directory = tempDirectory("tassel-pheno-join")
        val output = directory.resolve("joined")
        val result = run("join-phenotypes --input $traits --input $structure --output $output")
        assertEquals(0, result.statusCode, result.output)

        val expected = joinPhenotypes(listOf(readPhenotype(traits), readPhenotype(structure)))
        assertEquals(text(expected), outputText(directory.resolve("joined.txt")))
        assertContains(outputText(directory.resolve("joined.txt")), "EarHT")
        assertContains(outputText(directory.resolve("joined.txt")), "Q1")
    }

    @Test
    fun joinPhenotypesWithOneInputIsAUserError() {
        val result = run("join-phenotypes --input $traits")
        assertEquals(EXIT_USER_ERROR, result.statusCode)
        assertContains(result.output, "at least two phenotype tables")
    }

    @Test
    fun joinGenotypePhenotypeAcceptsAGenotypeAndAPhenotypeInEitherOrder() {
        val directory = tempDirectory("tassel-geno-pheno")
        val first = directory.resolve("a")
        val second = directory.resolve("b")
        assertEquals(0, run("join-genotype-phenotype --input $genotypes --input $traits --output $first").statusCode)
        assertEquals(0, run("join-genotype-phenotype --input $traits --input $genotypes --output $second").statusCode)

        val aligned = outputText(directory.resolve("a.txt"))
        assertEquals(aligned, outputText(directory.resolve("b.txt")))
        assertTrue(aligned.startsWith("<Phenotype>"), aligned.take(40))
    }

    @Test
    fun aGenotypeFileGivenToAPhenotypeOperationIsAUserError() {
        val result = run("phenotype-summary --input $genotypes")
        assertEquals(EXIT_USER_ERROR, result.statusCode)
        assertContains(result.output, "is not a phenotype file")
        assertContains(result.output, "phenotype-summary")
    }

    @Test
    fun aPhenotypeFileGivenToAGenotypeOperationIsAUserError() {
        val result = run("summary --input $traits")
        assertEquals(EXIT_USER_ERROR, result.statusCode)
        assertContains(result.output, "is a phenotype file")
    }

    @Test
    fun aMissingPhenotypeFileIsAUserError() {
        val result = run("phenotype-summary --input no/such/traits.txt")
        assertEquals(EXIT_USER_ERROR, result.statusCode)
        assertContains(result.output, "no such file")
    }

    @Test
    fun theTagDecidesNotTheSuffix() {
        val directory = tempDirectory("tassel-sniff")
        val renamed = directory.resolve("traits.hmp.txt")
        java.nio.file.Files.copy(java.nio.file.Path.of(traits), renamed)
        assertEquals(InputKind.Phenotype, sniffInput(renamed.toString()))
        assertEquals(InputKind.Genotype, sniffInput(genotypes))
        assertEquals(0, run("phenotype-summary --input $renamed --output ${directory.resolve("s")}").statusCode)
    }

    @Test
    fun regionNeedsAGenotypeAmongTheInputs() {
        val result = run("join-phenotypes --input $traits --input $structure --region chr1:1-10")
        assertEquals(EXIT_USER_ERROR, result.statusCode)
        assertContains(result.output, "--region needs a genotype")
    }
}
