package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class BgzfBenchTest {
    @Test
    fun printsInflatedBytesAndAMedian() {
        val input = testDataPath("htslib/tutorial/mdp_genotype_w_depth.vcf.gz")
        val result = tasselCommand().test("--quiet dev bgzf-bench --input $input --workers 1 --runs 1")
        assertEquals(0, result.statusCode, result.output)
        assertTrue(result.output.contains("bytes="), result.output)
        assertTrue(result.output.contains("workers=1"), result.output)
        assertTrue(result.output.contains("median="), result.output)
    }
}
