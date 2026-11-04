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
        assertEquals(27, manifest.outputs.size)
    }

    @Test
    fun everyGoldenOutputIsReadable() {
        for (path in GoldenManifest.load().outputs.keys) {
            assertTrue(TestData.readBytes(path).isNotEmpty(), "$path is empty")
        }
    }

    @Test
    fun parsesTolerances() {
        assertEquals(Tolerance.Format(), Tolerance.parse("format"))
        assertEquals(Tolerance.Relative(1e-9), Tolerance.parse("rel 1e-9"))
        assertEquals(Tolerance.Absolute(1e-5, upToSign = true), Tolerance.parse("abs 1e-5 up-to-sign"))
        assertEquals(Tolerance.Exact(rowsUnordered = true), Tolerance.parse("exact rows-unordered"))
        assertEquals(Tolerance.Relative(1e-9, rowsUnordered = true), Tolerance.parse("rel 1e-9 rows-unordered"))
        assertEquals(
            Tolerance.Absolute(1e-5, upToSign = true, rowsUnordered = true),
            Tolerance.parse("abs 1e-5 up-to-sign rows-unordered"),
        )
        assertFailsWith<IllegalStateException> { Tolerance.parse("abs") }
        assertFailsWith<IllegalArgumentException> { Tolerance.parse("exact up-to-sign") }
        assertFailsWith<IllegalArgumentException> { Tolerance.parse("exact rows-unordered rows-unordered") }
        assertFailsWith<IllegalArgumentException> { GoldenManifest.parse("version = 1") }
    }

    private companion object {
        val HAPMAP_FIXED_COLUMNS = listOf(
            "rs#", "alleles", "chrom", "pos", "strand", "assembly#",
            "center", "protLSID", "assayLSID", "panel", "QCcode",
        )
        const val HEADER_LENGTH = 1642
    }
}
