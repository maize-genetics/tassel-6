package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class BgzfWriteBenchTest {
    @Test
    fun printsCompressedBytesAndAMedian() {
        val input = testDataPath("tutorial/mdp_genotype_w_depth.vcf")
        val result = tasselCommand().test(
            "--quiet dev bgzf-write-bench --input $input --layout vcf --workers 1 --runs 1",
        )
        assertEquals(0, result.statusCode, result.output)
        assertTrue(result.output.contains("bytes="), result.output)
        assertTrue(result.output.contains("workers=1"), result.output)
        assertTrue(result.output.contains("layout=vcf"), result.output)
        assertTrue(result.output.contains("median="), result.output)
    }
}
