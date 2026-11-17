package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class BcfBenchTest {
    @Test
    fun printsBcfSitesAndAMedian() {
        val input = testDataPath("htslib/bcf/small/small.bcf")
        val result = tasselCommand().test("--quiet dev bcf-bench --bcf $input --runs 1")
        assertEquals(0, result.statusCode, result.output)
        assertTrue(result.output.contains("sites=1"), result.output)
        assertTrue(result.output.contains("reader=bcf"), result.output)
        assertTrue(result.output.contains("median="), result.output)
    }

    @Test
    fun printsVcfSitesAndAMedian() {
        val input = testDataPath("htslib/tutorial/mdp_genotype_w_depth.vcf.gz")
        val result = tasselCommand().test("--quiet dev bcf-bench --vcf $input --runs 1")
        assertEquals(0, result.statusCode, result.output)
        assertTrue(result.output.contains("sites=3093"), result.output)
        assertTrue(result.output.contains("reader=vcf"), result.output)
        assertTrue(result.output.contains("median="), result.output)
    }

    @Test
    fun requiresExactlyOneInput() {
        val missing = tasselCommand().test("--quiet dev bcf-bench --runs 1")
        assertEquals(1, missing.statusCode, missing.output)
        val bcf = testDataPath("htslib/bcf/small/small.bcf")
        val vcf = testDataPath("htslib/tutorial/mdp_genotype_w_depth.vcf.gz")
        val both = tasselCommand().test("--quiet dev bcf-bench --bcf $bcf --vcf $vcf --runs 1")
        assertEquals(1, both.statusCode, both.output)
    }
}
