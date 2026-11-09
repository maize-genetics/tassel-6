package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

/** Whole-file BGZF entry point. One worker on wasmJs; the host count on the JVM and on Native. */
class BgzfVcfReadTest {
    @Test
    fun bgzippedTutorialMatchesThePlainVcf() = runTest {
        val expected = Vcf.read(Buffer().also { it.write(TestData.readBytes(PLAIN_VCF)) })
        ByteArrayRandomAccessSource(TestData.readBytes(TUTORIAL_GZ)).use { source ->
            val read = Vcf.read(source)
            assertEquals(emptyList(), read.warnings)
            assertSameTable(expected, read.table)
        }
    }

    @Test
    fun missingEofIsAWarningAndTheSitesRemain() = runTest {
        ByteArrayRandomAccessSource(TestData.readBytes(NO_EOF)).use { source ->
            val read = Vcf.read(source)
            assertEquals(listOf(MISSING_BGZF_EOF), read.warnings)
            assertEquals(2, read.table.numberOfSites)
            assertEquals(2, read.table.numberOfTaxa)
        }
    }

    @Test
    fun bgzfIsDetectedFromTheFirstMember() = runTest {
        ByteArrayRandomAccessSource(TestData.readBytes(TUTORIAL_GZ)).use { source ->
            assertTrue(source.isBgzf())
        }
        ByteArrayRandomAccessSource(TestData.readBytes(PLAIN_GZ)).use { source ->
            assertFalse(source.isBgzf())
        }
    }

    private fun assertSameTable(expected: GenotypeTable, actual: GenotypeTable) {
        assertEquals(expected.taxa.names, actual.taxa.names)
        assertEquals(expected.numberOfSites, actual.numberOfSites)
        for (site in 0 until expected.numberOfSites) {
            assertContentEquals(expected.siteGenotypes(site), actual.siteGenotypes(site), "site $site")
        }
    }

    private companion object {
        const val PLAIN_VCF = "tutorial/mdp_genotype_w_depth.vcf"
        const val TUTORIAL_GZ = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
        const val NO_EOF = "htslib/variants/no-eof.vcf.gz"
        const val PLAIN_GZ = "htslib/variants/plain.vcf.gz"
    }
}
