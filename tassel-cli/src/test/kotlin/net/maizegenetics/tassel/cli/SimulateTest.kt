package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlinx.coroutines.runBlocking
import kotlinx.io.buffered
import kotlinx.io.files.SystemFileSystem
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.io.HapMap
import java.nio.file.Path
import kotlin.io.path.readBytes
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class SimulateTest {
    private fun simulate(arguments: String, into: Path): Path {
        val result = tasselCommand().test("--quiet dev simulate $arguments --output $into")
        assertEquals(0, result.statusCode, result.output)
        return into
    }

    @Test
    fun theSameSeedWritesTheSameBytes() = runBlocking {
        val directory = tempDirectory("tassel-simulate")
        val first = simulate("--taxa 25 --sites 200 --seed 1", directory.resolve("first.hmp.txt"))
        val second = simulate("--taxa 25 --sites 200 --seed 1", directory.resolve("second.hmp.txt"))
        assertContentEquals(first.readBytes(), second.readBytes())
    }

    @Test
    fun anotherSeedWritesAnotherFile() = runBlocking {
        val directory = tempDirectory("tassel-simulate-seeds")
        val first = simulate("--taxa 25 --sites 200 --seed 1", directory.resolve("first.hmp.txt"))
        val second = simulate("--taxa 25 --sites 200 --seed 2", directory.resolve("second.hmp.txt"))
        assertFalse(first.readBytes().contentEquals(second.readBytes()), "seeds 1 and 2 gave the same file")
    }

    @Test
    fun theOutputLoadsBackThroughTheHapMapReader() = runBlocking {
        val directory = tempDirectory("tassel-simulate-read")
        val file = simulate("--taxa 25 --sites 200 --seed 1 --missing-rate 0.1", directory.resolve("sim.hmp.txt"))

        val table = SystemFileSystem.source(kotlinx.io.files.Path(file.toString())).buffered()
            .use { HapMap.read(it) }
        assertEquals(25, table.numberOfTaxa)
        assertEquals(200, table.numberOfSites)
        assertEquals("taxon_01", table.taxa[0].name)
        assertTrue(table.positions.all { it.chromosome.name == "1" })
        assertEquals(1_000, table.positions[0].position)
        val missing = (0 until table.numberOfSites).sumOf { site ->
            table.siteGenotypes(site).count { it == Genotype.UNKNOWN }
        }
        assertTrue(missing in 300..700, "expected about 10% of 5,000 genotypes missing, got $missing")
    }

    @Test
    fun noGenotypeIsMissingByDefault() = runBlocking {
        val table = simulateGenotypes(taxa = 10, sites = 20, seed = 7)
        val missing = (0 until table.numberOfSites).sumOf { site ->
            table.siteGenotypes(site).count { it == Genotype.UNKNOWN }
        }
        assertEquals(0, missing)
    }

    @Test
    fun minorAlleleFrequenciesStayInsideTheRequestedRange() = runBlocking {
        val table = simulateGenotypes(taxa = 400, sites = 60, minMaf = 0.3, maxMaf = 0.5, seed = 11)
        val frequencies = (0 until table.numberOfSites).map { table.siteStats(it).minorAlleleFrequency }
        // 800 gametes per site, so sampling noise around a frequency in [0.3, 0.5] is small.
        assertTrue(frequencies.all { it in 0.22..0.5 }, "unexpected minor allele frequencies: $frequencies")
    }

    @Test
    fun impossibleOptionsAreUserErrors() = runBlocking {
        val result = tasselCommand().test("dev simulate --min-maf 0.4 --max-maf 0.1")
        assertEquals(EXIT_USER_ERROR, result.statusCode)
        assertContains(result.output, "--min-maf")
    }
}
