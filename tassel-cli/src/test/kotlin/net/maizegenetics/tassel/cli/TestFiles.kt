package net.maizegenetics.tassel.cli

import kotlinx.coroutines.runBlocking
import kotlinx.io.buffered
import kotlinx.io.files.SystemFileSystem
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.HapMap
import java.nio.file.Files
import java.nio.file.Path
import kotlin.io.path.readText
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals

/** `testdata/`, located by the environment variable `tassel-cli/build.gradle.kts` sets. */
private val testDataDir: String = requireNotNull(System.getenv("TASSEL_TESTDATA_DIR")) {
    "TASSEL_TESTDATA_DIR is not set; run these tests through Gradle"
}

/** Absolute path of a file under `testdata/`, for example `tutorial/mdp_genotype.hmp.txt`. */
fun testDataPath(relative: String): String = "$testDataDir/$relative"

/** Text of a file under `testdata/`, with line endings normalized as the core's golden tests do. */
fun testDataText(relative: String): String = normalizeLineEndings(Path.of(testDataPath(relative)).readText())

/** Text of a file the CLI wrote, with line endings normalized. */
fun outputText(path: Path): String = normalizeLineEndings(path.readText())

fun normalizeLineEndings(text: String): String = text.replace("\r\n", "\n").replace('\r', '\n')

/** A temporary directory that is deleted when the JVM exits. */
fun tempDirectory(name: String): Path = Files.createTempDirectory(name).also { it.toFile().deleteOnExit() }

/** Writes [table] as HapMap, for tests that need an input file. */
fun writeHapMap(table: GenotypeTable, path: Path) = runBlocking {
    SystemFileSystem.sink(kotlinx.io.files.Path(path.toString())).buffered().use { HapMap.write(table, it) }
}

/**
 * Taxa, positions, and site genotype bytes. This is the comparison [net.maizegenetics.tassel.io.Vcf]
 * uses when a written file is read back.
 */
fun assertSameGenotypeTable(expected: GenotypeTable, actual: GenotypeTable) {
    assertEquals(expected.taxa.names, actual.taxa.names)
    assertEquals(expected.numberOfSites, actual.numberOfSites)
    for (site in 0 until expected.numberOfSites) {
        assertEquals(expected.positions[site], actual.positions[site], "site $site")
        assertContentEquals(expected.siteGenotypes(site), actual.siteGenotypes(site), "site $site")
    }
}
