package net.maizegenetics.tassel.testing

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNotEquals
import kotlin.test.assertTrue

/**
 * Every file `testdata/htslib/MANIFEST` lists is readable on the JVM, on Node, and on Native.
 * Empty region oracles are valid: htslib 1.25 prints nothing for an unknown contig, an inverted
 * region, and a range with no records. The BCF section records the wheat `chr6D` contig-list
 * comparison and that each BCF region query matched the stored `.vcf.gz` region text.
 * `[golden-bcf]` is htslib's encoding of each TASSEL 5 VCF-export golden. The test hashes those
 * goldens, so a golden byte change fails here until the fixtures are regenerated.
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
        assertEquals(84, manifest.files.size)

        val bcf = manifest.section("bcf")
        assertEquals("true", bcf["chr6d_in_bcf_header"])
        assertEquals("false", bcf["chr6d_in_bcf_index_s"])
        assertEquals("false", bcf["chr6d_in_vcf_index_s"])
        assertEquals("false", bcf["chr6d_in_vcf_csi"])
        assertEquals("false", bcf["header_matches_vcf_index"])
        assertEquals("0", bcf["wheat_bcf_csi_l_aux"])
        assertEquals("6", bcf["tutorial_bcf_csi_depth"])
        assertEquals("14", bcf["tutorial_bcf_csi_min_shift"])
        assertEquals(34, bcf.identical.size)

        val layout = manifest.section("bgzf-layout")
        assertEquals(
            listOf(
                "htslib/bgzf-layout/long-line.vcf.gz",
                "htslib/bgzf-layout/long-header.vcf.gz",
            ),
            layout.files,
        )
        assertTrue(layout["long_line_bytes"].toInt() > 65_536)
        assertTrue(layout["long_header_bytes"].toInt() > 65_536)
        assertTrue(layout["long_line_data_blocks"].toInt() >= 2)
        assertTrue(layout["long_header_data_blocks"].toInt() >= 2)

        val golden = manifest.section("golden-bcf")
        assertEquals("vcf_export_edges", golden["smallest_golden"])
        assertEquals("false", golden["level0_blocks_match"])
        assertEquals(
            listOf(
                "htslib/bcf/golden/mdp_genotype.bcf",
                "htslib/bcf/golden/mdp_genotype.blocks.txt",
                "htslib/bcf/golden/include-list.bcf",
                "htslib/bcf/golden/include-list.blocks.txt",
                "htslib/bcf/golden/vcf_export_edges.bcf",
                "htslib/bcf/golden/vcf_export_edges.blocks.txt",
                "htslib/bcf/golden/vcf_export_edges.blocks-l0.txt",
            ),
            golden.files,
        )
        val level0Path = "htslib/bcf/golden/${golden["smallest_golden"]}.blocks-l0.txt"
        if (golden["level0_blocks_match"] == "false") {
            assertTrue(level0Path in golden.files)
        } else {
            assertTrue(golden.files.none { it.endsWith(".blocks-l0.txt") })
        }
        val names = listOf("mdp_genotype", "include-list", "vcf_export_edges")
        var smallest = names.first()
        var smallestSize = Int.MAX_VALUE
        for (name in names) {
            val path = golden["${name}_golden"]
            val bytes = TestData.readBytes(path)
            assertEquals(golden["${name}_sha256"], sha256(bytes), path)
            if (bytes.size < smallestSize) {
                smallest = name
                smallestSize = bytes.size
            }
            val bcf = "htslib/bcf/golden/$name.bcf"
            val blocks = blockSizes("htslib/bcf/golden/$name.blocks.txt")
            assertEquals(28, blocks.last(), bcf)
            assertEquals(TestData.readBytes(bcf).size, blocks.sum(), bcf)
        }
        assertEquals(smallest, golden["smallest_golden"])
        val level0 = blockSizes(level0Path)
        assertEquals(28, level0.last())
        assertNotEquals(blockSizes("htslib/bcf/golden/${golden["smallest_golden"]}.blocks.txt"), level0)

        for (path in manifest.files) {
            TestData.readBytes(path)
        }

        assertTrue(TestData.readBytes("htslib/tutorial/mdp_genotype_w_depth.vcf.gz").isNotEmpty())
        assertTrue(TestData.readBytes("htslib/wheat/regions/chr3B_536870002_536870010.record.txt").isNotEmpty())
        assertEquals(0, TestData.readBytes("htslib/wheat/regions/chr3B_536870002_536870010.pos.txt").size)
        assertEquals(0, TestData.readBytes("htslib/wheat/regions/chrUn_1_1000.record.txt").size)
        assertEquals(0, TestData.readBytes("htslib/wheat/regions/chr3B_5000_1000.record.txt").size)
        assertTrue(TestData.readBytes("htslib/large-pos/regions/chr1_3000000000_3000000000.record.txt").isNotEmpty())

        val raw = TestData.readBytes("htslib/bcf/small/small.raw.bcf")
        assertEquals(0x42, raw[0].toInt() and 0xff)
        assertEquals(0x43, raw[1].toInt() and 0xff)
        assertEquals(0x46, raw[2].toInt() and 0xff)
        assertEquals(0x02, raw[3].toInt() and 0xff)
        assertEquals(0x02, raw[4].toInt() and 0xff)
        val stored = TestData.readBytes("htslib/bcf/small/small.stored.bcf")
        assertEquals(0x1f, stored[0].toInt() and 0xff)
        assertEquals(0x8b, stored[1].toInt() and 0xff)
        assertEquals(0x01, stored[18].toInt() and 0xff)
    }

    @Test
    fun rejectsALineOutsideASection() {
        assertFailsWith<IllegalArgumentException> {
            HtslibManifest.parse("version = 1.25\n")
        }
    }

    @Test
    fun identicalMayRepeat() {
        val manifest = HtslibManifest.parse(
            "[bcf]\n" +
                "identical = htslib/wheat/regions/a.txt\n" +
                "identical = htslib/wheat/regions/b.txt\n",
        )
        assertEquals(
            listOf("htslib/wheat/regions/a.txt", "htslib/wheat/regions/b.txt"),
            manifest.section("bcf").identical,
        )
    }

    private fun blockSizes(path: String): List<Int> =
        TestData.readText(path).lineSequence().filter { it.isNotEmpty() }.map { line ->
            val size = line.toInt()
            assertTrue(size > 0, path)
            size
        }.toList()
}
