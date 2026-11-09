package net.maizegenetics.tassel.testing

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Every file `testdata/htslib/MANIFEST` lists is readable on the JVM, on Node, and on Native.
 * Empty region oracles are valid: htslib 1.25 prints nothing for an unknown contig, an inverted
 * region, and a range with no records.
 */
class HtslibFixtureTest {
    @Test
    fun everyManifestFileIsReadable() {
        val manifest = HtslibManifest.load()
        assertEquals("1.25", manifest.section("htslib")["version"])
        assertEquals("8", manifest.section("tutorial")["csi_depth"])
        assertEquals("6", manifest.section("wheat")["csi_depth"])
        assertEquals("6", manifest.section("large-pos")["csi_depth"])
        assertEquals("14", manifest.section("tutorial")["csi_min_shift"])
        assertEquals("14", manifest.section("wheat")["csi_min_shift"])
        assertEquals("14", manifest.section("large-pos")["csi_min_shift"])
        assertEquals(54, manifest.files.size)

        for (path in manifest.files) {
            TestData.readBytes(path)
        }

        assertTrue(TestData.readBytes("htslib/tutorial/mdp_genotype_w_depth.vcf.gz").isNotEmpty())
        assertTrue(TestData.readBytes("htslib/wheat/regions/chr3B_536870002_536870010.record.txt").isNotEmpty())
        assertEquals(0, TestData.readBytes("htslib/wheat/regions/chr3B_536870002_536870010.pos.txt").size)
        assertEquals(0, TestData.readBytes("htslib/wheat/regions/chrUn_1_1000.record.txt").size)
        assertEquals(0, TestData.readBytes("htslib/wheat/regions/chr3B_5000_1000.record.txt").size)
        assertTrue(TestData.readBytes("htslib/large-pos/regions/chr1_3000000000_3000000000.record.txt").isNotEmpty())
    }

    @Test
    fun rejectsALineOutsideASection() {
        assertFailsWith<IllegalArgumentException> {
            HtslibManifest.parse("version = 1.25\n")
        }
    }
}
