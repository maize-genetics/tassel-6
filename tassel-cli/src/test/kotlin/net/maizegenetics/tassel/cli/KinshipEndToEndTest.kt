package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.testing.test
import kotlin.io.path.exists
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * `tassel kinship` on the tutorial HapMap file against the TASSEL 5 golden.
 *
 * The default kernel is the blocked dosage product (ADR-0020), so the matrix is held to the
 * MANIFEST tolerance `abs 8e-6` and `SumPk` to `abs 1e-4`. Byte-identical text is the retained
 * TASSEL 5-order kernel, covered by `CenteredIbsTest.writesTheGoldenByteForByte`.
 */
class KinshipEndToEndTest {
    @Test
    fun matchesTheTassel5GoldenWithinTolerance() {
        val tutorial = testDataPath("tutorial/mdp_genotype.hmp.txt")
        val directory = tempDirectory("tassel-kinship")
        val base = directory.resolve("centered_ibs")
        val result = tasselCommand().test("--quiet kinship --input $tutorial --output $base")

        assertEquals(0, result.statusCode, result.output)
        // TASSEL 5 writes this with -export <base> -exportType SqrMatrix, which appends .txt.
        val written = directory.resolve("centered_ibs.txt")
        assertTrue(written.exists(), "kinship wrote no centered_ibs.txt")
        assertWithinTolerance(
            testDataText("golden/kinship/centered_ibs.txt"),
            outputText(written),
        )
    }
}

private fun assertWithinTolerance(expectedText: String, actualText: String) {
    val expected = matrixLines(expectedText)
    val actual = matrixLines(actualText)
    assertEquals(expected.size, actual.size, "line count")
    assertEquals(expected[0].substringBefore('='), actual[0].substringBefore('='), "SumPk key")
    val expectedSum = expected[0].substringAfter('=').toDouble()
    val actualSum = actual[0].substringAfter('=').toDouble()
    assertTrue(abs(expectedSum - actualSum) <= 1e-4, "SumPk $actualSum vs $expectedSum")
    assertEquals(expected[1], actual[1], "Matrix_Type")
    assertEquals(expected[2], actual[2], "taxa count")
    val taxa = expected[2].toInt()
    for (row in 0 until taxa) {
        val expectedCells = expected[3 + row].split('\t')
        val actualCells = actual[3 + row].split('\t')
        assertEquals(expectedCells[0], actualCells[0], "row label ${row + 1}")
        assertEquals(taxa + 1, actualCells.size, "row ${row + 1} width")
        for (column in 1..taxa) {
            val diff = abs(expectedCells[column].toDouble() - actualCells[column].toDouble())
            assertTrue(diff <= 8e-6, "row ${row + 1} column $column differs by $diff")
        }
    }
    assertEquals("", actual.last(), "the file ends with a newline")
}

private fun matrixLines(text: String): List<String> = text.lines()
