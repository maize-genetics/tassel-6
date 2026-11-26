package net.maizegenetics.tassel.testing

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class TestDataTest {
    @Test
    fun readsTutorialHapMapHeader() {
        val header = TestData.readText("tutorial/mdp_genotype.hmp.txt").substringBefore('\n').removeSuffix("\r")
        val columns = header.split('\t')
        assertEquals(HAPMAP_FIXED_COLUMNS, columns.take(HAPMAP_FIXED_COLUMNS.size))
        assertEquals(281, columns.size - HAPMAP_FIXED_COLUMNS.size)
        assertEquals("33-16", columns[HAPMAP_FIXED_COLUMNS.size])
    }

    @Test
    fun readsBytesUnchanged() {
        val bytes = TestData.readBytes("tutorial/mdp_genotype.hmp.txt")
        assertEquals(2_796_284, bytes.size)
        assertEquals("\r\n", bytes.copyOfRange(HEADER_LENGTH, HEADER_LENGTH + 2).decodeToString())
    }

    @Test
    fun parsesGoldenManifest() {
        val manifest = GoldenManifest.load()
        val tassel5 = manifest.section("tassel5")
        assertEquals("5.2.98", tassel5["version"])
        assertEquals("v5.2.98", tassel5["tag"])
        assertEquals(Tolerance.Absolute(8e-6), manifest.tolerance("golden/kinship/centered_ibs.txt"))
        assertEquals(Tolerance.Exact(), manifest.tolerance("golden/summary/sites.txt"))
        assertEquals(Tolerance.Absolute(1e-9, upToSign = true), manifest.tolerance("golden/pca/eigenvectors.full.txt"))
        assertEquals(2, manifest.section("pca").commands.size)
        assertEquals(OUTPUT_COUNT, manifest.outputs.size)
    }

    @Test
    fun everyGoldenOutputIsReadable() {
        for (path in GoldenManifest.load().outputs.keys) {
            assertTrue(TestData.readBytes(path).isNotEmpty(), "$path is empty")
        }
    }

    /** The GLM goldens (#188) are reports with the header TASSEL 5 writes; a smoke test that they parse. */
    @Test
    fun glmGoldensHaveTheExpectedHeaders() {
        val stats = "Trait\tMarker\tChr\tPos\tmarker_F\tp\tmarker_Rsq\tadd_F\tadd_p\tdom_F\tdom_p\tmarker_df\t" +
            "marker_MS\terror_df\terror_MS\tmodel_df\tmodel_MS\tminorObs"
        val effects = "Trait\tMarker\tChr\tPos\tObs\tAllele\tEstimate"
        val manifest = GoldenManifest.load()
        val reports = manifest.outputs.keys.filter {
            it.startsWith("golden/glm-") && it.endsWith(".txt") && "README" !in it
        }
        assertEquals(57, reports.size, "GLM report count")
        for (path in reports) {
            val lines = TestData.readText(path).trimEnd('\n').split('\n')
            val header = lines.first().trimEnd('\r')
            when {
                "refprob" in path -> assertTrue(header.startsWith("Trait\tMarker\tChr\tPos\t"), "$path: $header")
                ".effects" in path -> assertEquals(effects, header, path)
                else -> assertTrue(header == stats || header.startsWith("$stats\t"), "$path: $header")
            }
            assertTrue(lines.size > 1, "$path has no rows")
            val width = header.split('\t').size
            assertTrue(lines.drop(1).all { it.trimEnd('\r').split('\t').size == width }, "$path has ragged rows")
        }
        assertEquals(Tolerance.Format(), manifest.tolerance("golden/glm-tutorial/mdp_traits.stats.txt"))
        assertEquals(Tolerance.Relative(1e-8), manifest.tolerance("golden/glm-tutorial/mdp_traits.stats.full.txt"))
    }

    @Test
    fun parsesTolerances() {
        assertEquals(Tolerance.Format(), Tolerance.parse("format"))
        assertEquals(Tolerance.Relative(1e-9), Tolerance.parse("rel 1e-9"))
        assertEquals(Tolerance.Absolute(1e-5, upToSign = true), Tolerance.parse("abs 1e-5 up-to-sign"))
        assertEquals(Tolerance.Exact(rowsUnordered = true), Tolerance.parse("exact rows-unordered"))
        assertEquals(Tolerance.Phenotype(), Tolerance.parse("phenotype"))
        assertEquals(Tolerance.Phenotype(rowsUnordered = true), Tolerance.parse("phenotype rows-unordered"))
        assertEquals(Tolerance.Phenotype(relative = 1e-6), Tolerance.parse("phenotype rel 1e-6"))
        assertEquals(
            Tolerance.Phenotype(rowsUnordered = true, relative = 1e-6),
            Tolerance.parse("phenotype rel 1e-6 rows-unordered"),
        )
        assertEquals(
            Tolerance.Phenotype(relative = 1e-6, absolute = 1e-5),
            Tolerance.parse("phenotype rel 1e-6 abs 1e-5"),
        )
        assertEquals(Tolerance.Phenotype(absolute = 1e-5), Tolerance.parse("phenotype abs 1e-5"))
        assertFailsWith<IllegalArgumentException> { Tolerance.parse("phenotype rel 1e-6 rel 1e-5") }
        assertFailsWith<IllegalArgumentException> { Tolerance.parse("phenotype exact 1e-6") }
        assertFailsWith<IllegalArgumentException> { Tolerance.parse("phenotype rel") }
        assertFailsWith<IllegalStateException> { Tolerance.parse("phenotype rel x") }
        assertEquals(Tolerance.Relative(1e-9, rowsUnordered = true), Tolerance.parse("rel 1e-9 rows-unordered"))
        assertEquals(
            Tolerance.Absolute(1e-5, upToSign = true, rowsUnordered = true),
            Tolerance.parse("abs 1e-5 up-to-sign rows-unordered"),
        )
        assertFailsWith<IllegalStateException> { Tolerance.parse("abs") }
        assertFailsWith<IllegalArgumentException> { Tolerance.parse("exact up-to-sign") }
        assertFailsWith<IllegalArgumentException> { Tolerance.parse("phenotype up-to-sign") }
        assertFailsWith<IllegalArgumentException> { Tolerance.parse("exact rows-unordered rows-unordered") }
        assertFailsWith<IllegalArgumentException> { GoldenManifest.parse("version = 1") }
    }

    private companion object {
        const val OUTPUT_COUNT = 105

        val HAPMAP_FIXED_COLUMNS = listOf(
            "rs#", "alleles", "chrom", "pos", "strand", "assembly#",
            "center", "protLSID", "assayLSID", "panel", "QCcode",
        )
        const val HEADER_LENGTH = 1642
    }
}
