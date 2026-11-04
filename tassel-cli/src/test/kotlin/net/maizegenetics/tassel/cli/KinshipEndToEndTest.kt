package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlin.io.path.exists
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * `tassel kinship` on the tutorial HapMap file against the TASSEL 5 golden.
 *
 * The comparison is exact, not the MANIFEST's `abs 8e-6`: tassel-6#18 made the written text byte
 * identical to TASSEL 5's, so comparing text catches a writer regression that a numeric tolerance
 * would hide.
 */
class KinshipEndToEndTest {
    @Test
    fun matchesTheTassel5GoldenExactly() {
        val tutorial = testDataPath("tutorial/mdp_genotype.hmp.txt")
        val directory = tempDirectory("tassel-kinship")
        val base = directory.resolve("centered_ibs")
        val result = tasselCommand().test("--quiet kinship --input $tutorial --output $base")

        assertEquals(0, result.statusCode, result.output)
        // TASSEL 5 writes this with -export <base> -exportType SqrMatrix, which appends .txt.
        val written = directory.resolve("centered_ibs.txt")
        assertTrue(written.exists(), "kinship wrote no centered_ibs.txt")
        assertEquals(testDataText("golden/kinship/centered_ibs.txt"), outputText(written))
    }
}
