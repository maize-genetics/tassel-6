package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readString
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.testing.GoldenManifest
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.Tolerance
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * A whole BCF file against [Vcf.read] and, for the tutorial, the TASSEL 5 HapMap export (ADR-0022).
 *
 * The #118 files are the oracle. `small.raw.bcf` is raw BCF. `small.stored.bcf` is BGZF with stored blocks,
 * `bcftools view -Ou`, so the bytes at offset 0 are gzip rather than the BCF magic.
 */
class BcfReadTest {
    @Test
    fun tutorialMatchesTheGoldenAndTheVcf() = runTest {
        assertEquals(Tolerance.Exact(), GoldenManifest.load().tolerance(GOLDEN))
        val expected = readVcfGz(TUTORIAL_VCF)
        val fromSource = readSource(TUTORIAL_BCF)
        val fromRandom = readRandom(TUTORIAL_BCF)
        assertEquals(emptyList(), fromSource.warnings)
        assertEquals(emptyList(), fromRandom.warnings)
        assertGolden(fromSource.table)
        assertGolden(fromRandom.table)
        assertSameTable(expected, fromSource.table, "source")
        assertSameTable(expected, fromRandom.table, "random")
    }

    @Test
    fun wheatMatchesTheVcf() = runTest {
        val expected = readVcfGz(WHEAT_VCF)
        assertBothMatch(WHEAT_BCF, expected)
    }

    @Test
    fun smallFilesMatchTheView() = runTest {
        val expected = readPlainVcf(SMALL_VIEW)
        assertBothMatch(SMALL_BCF, expected)
        assertBothMatch(SMALL_STORED, expected)
        val rawBytes = TestData.readBytes(SMALL_RAW)
        assertEquals("BCF", rawBytes.decodeToString(0, 3))
        assertEquals(2, rawBytes[3].toInt() and 0xff)
        assertEquals(2, rawBytes[4].toInt() and 0xff)
        val raw = Bcf.read(Buffer().also { it.write(rawBytes) })
        assertEquals(emptyList(), raw.warnings)
        assertSameTable(expected, raw.table, SMALL_RAW)
    }

    @Test
    fun storedBcfStartsWithGzipMagic() {
        val stored = TestData.readBytes(SMALL_STORED)
        assertEquals(0x1f, stored[0].toInt() and 0xff)
        assertEquals(0x8b, stored[1].toInt() and 0xff)
    }

    @Test
    fun edgeAndBoundaryMatchTheView() = runTest {
        assertBothMatch(EDGE_BCF, readPlainVcf(EDGE_VIEW))
        assertBothMatch(BOUNDARY_BCF, readPlainVcf(BOUNDARY_VIEW))
    }

    @Test
    fun noSampleHasNoTaxa() = runTest {
        val column = TestData.readText(NO_SAMPLE_VIEW).lineSequence()
            .first { it.isNotEmpty() && !it.startsWith("#") }
            .split('\t')
        val siteName = "${column[0]}:${column[1]}"
        for (read in listOf(readSource(NO_SAMPLE_BCF), readRandom(NO_SAMPLE_BCF))) {
            assertEquals(emptyList(), read.warnings)
            assertEquals(0, read.table.numberOfTaxa)
            assertEquals(1, read.table.numberOfSites)
            val site = read.table.positions[0]
            assertEquals(Chromosome.of(column[0]), site.chromosome)
            assertEquals(column[1].toInt(), site.position)
            assertEquals(splitVariants(variantsOf(column[3], column[4]), siteName), site.knownVariants)
        }
    }

    @Test
    fun missingEofIsAWarningAndTheSitesRemain() = runTest {
        val expected = readPlainVcf(SMALL_VIEW)
        for (read in listOf(readSource(NO_EOF), readRandom(NO_EOF))) {
            assertEquals(listOf(MISSING_BGZF_EOF), read.warnings)
            assertSameTable(expected, read.table, NO_EOF)
        }
    }

    @Test
    fun truncatedRecordThrows() = runTest {
        val fromSource = assertFailsWith<IllegalArgumentException> { readSource(TRUNCATED) }
        assertEquals(TRUNCATED_BCF_RECORD, fromSource.message)
        val fromRandom = assertFailsWith<IllegalArgumentException> { readRandom(TRUNCATED) }
        assertEquals(TRUNCATED_BCF_RECORD, fromRandom.message)
    }

    @Test
    fun badMagicAndOldVersionAreRefused() = runTest {
        val sourceMagic = assertFailsWith<IllegalArgumentException> { readSource(BAD_MAGIC) }
        assertEquals(INVALID_BCF_MAGIC, sourceMagic.message)
        val randomMagic = assertFailsWith<IllegalArgumentException> { readRandom(BAD_MAGIC) }
        assertEquals(INVALID_BCF_MAGIC, randomMagic.message)
        val sourceVersion = assertFailsWith<IllegalArgumentException> { readSource(BCF21) }
        assertEquals(UNSUPPORTED_BCF_VERSION, sourceVersion.message)
        val randomVersion = assertFailsWith<IllegalArgumentException> { readRandom(BCF21) }
        assertEquals(UNSUPPORTED_BCF_VERSION, randomVersion.message)
    }

    private suspend fun assertBothMatch(path: String, expected: GenotypeTable) {
        val fromSource = readSource(path)
        val fromRandom = readRandom(path)
        assertEquals(emptyList(), fromSource.warnings, path)
        assertEquals(emptyList(), fromRandom.warnings, path)
        assertSameTable(expected, fromSource.table, "$path source")
        assertSameTable(expected, fromRandom.table, "$path random")
    }

    private suspend fun assertGolden(table: GenotypeTable) {
        val written = Buffer().also { HapMap.write(table, it) }.readString().lines()
        val goldenLines = TestData.readText(GOLDEN).replace("\r\n", "\n").lines()
        val skipped = goldenLines.takeWhile { it.startsWith("##") }
        assertEquals(table.numberOfTaxa, skipped.size, "skipped ##SAMPLE= lines")
        assertTrue(skipped.all { it.startsWith("##SAMPLE=<ID=") }, "every skipped line is a ##SAMPLE= annotation")
        val golden = goldenLines.drop(skipped.size)
        assertTrue(golden.first().startsWith("rs#\talleles\tchrom\tpos\t"), "first compared line: ${golden.first()}")
        assertEquals(golden.size, written.size, "number of lines")
        written.forEachIndexed { index, line ->
            if (line != golden[index]) reportLineDifference(skipped.size + index, golden[index], line)
        }
    }

    private fun assertSameTable(expected: GenotypeTable, actual: GenotypeTable, label: String) {
        assertEquals(expected.taxa.names, actual.taxa.names, label)
        assertEquals(expected.numberOfSites, actual.numberOfSites, label)
        for (site in 0 until expected.numberOfSites) {
            val want = expected.positions[site]
            val got = actual.positions[site]
            assertEquals(want, got, "$label site $site")
            assertEquals(want.knownVariants, got.knownVariants, "$label site $site variants")
            assertContentEquals(expected.siteGenotypes(site), actual.siteGenotypes(site), "$label site $site")
        }
    }

    private fun reportLineDifference(line: Int, expected: String, actual: String): Nothing {
        val column = expected.commonPrefixWith(actual).length
        throw AssertionError(
            "Line ${line + 1} differs at column ${column + 1}:\n" +
                "  expected: ...${expected.drop(column).take(40)}\n" +
                "  actual:   ...${actual.drop(column).take(40)}",
        )
    }

    private suspend fun readSource(path: String): BcfRead =
        Bcf.read(Buffer().also { it.write(TestData.readBytes(path)) })

    private suspend fun readRandom(path: String): BcfRead =
        ByteArrayRandomAccessSource(TestData.readBytes(path)).use { Bcf.read(it) }

    private suspend fun readPlainVcf(path: String): GenotypeTable =
        Vcf.read(Buffer().also { it.write(TestData.readBytes(path)) })

    private suspend fun readVcfGz(path: String): GenotypeTable =
        ByteArrayRandomAccessSource(TestData.readBytes(path)).use { Vcf.read(it).table }

    private companion object {
        const val TUTORIAL_BCF = "htslib/bcf/tutorial/mdp_genotype_w_depth.bcf"
        const val TUTORIAL_VCF = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
        const val GOLDEN = "golden/vcf/mdp_genotype_w_depth.hmp.txt"
        const val WHEAT_BCF = "htslib/bcf/wheat/wheat.bcf"
        const val WHEAT_VCF = "htslib/wheat/wheat.vcf.gz"
        const val SMALL_BCF = "htslib/bcf/small/small.bcf"
        const val SMALL_RAW = "htslib/bcf/small/small.raw.bcf"
        const val SMALL_STORED = "htslib/bcf/small/small.stored.bcf"
        const val SMALL_VIEW = "htslib/bcf/small/small.view.vcf"
        const val EDGE_BCF = "htslib/bcf/edge/edge.bcf"
        const val EDGE_VIEW = "htslib/bcf/edge/edge.view.vcf"
        const val BOUNDARY_BCF = "htslib/bcf/boundary/pos2147483647.bcf"
        const val BOUNDARY_VIEW = "htslib/bcf/boundary/pos2147483647.view.vcf"
        const val NO_SAMPLE_BCF = "htslib/bcf/edge/no-sample.bcf"
        const val NO_SAMPLE_VIEW = "htslib/bcf/edge/no-sample.view.vcf"
        const val NO_EOF = "htslib/bcf/corrupt/no-eof.bcf"
        const val TRUNCATED = "htslib/bcf/corrupt/truncated-record.bcf"
        const val BAD_MAGIC = "htslib/bcf/corrupt/bad-magic.bcf"
        const val BCF21 = "htslib/bcf/corrupt/bcf21.bcf"
    }
}
