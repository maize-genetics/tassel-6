package net.maizegenetics.tassel.testing

import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/** Each tolerance kind has to pass when it should and fail, readably, when it should not. */
class GoldenCompareTest {
    @Test
    fun exactAcceptsOnlyLineEndingDifferences() {
        assertDelimitedMatches(TABLE, TABLE.replace("\n", "\r\n"), Tolerance.Exact())
        val problem = failure { assertDelimitedMatches(TABLE, TABLE.replace("0.25", "0.26"), Tolerance.Exact()) }
        assertContains(problem, "line 2 differs")
        assertContains(problem, "0.26")
    }

    @Test
    fun exactCatchesAMissingFinalNewline() {
        val problem = failure { assertDelimitedMatches(TABLE, TABLE.trimEnd('\n'), Tolerance.Exact()) }
        assertContains(problem, "only the expected text continues")
    }

    @Test
    fun absoluteAllowsSmallDifferencesAndNamesTheCellThatIsTooFarOff() {
        val tolerance = Tolerance.Absolute(1e-3)
        assertDelimitedMatches(TABLE, TABLE.replace("0.25", "0.2505"), tolerance)
        val problem = failure { assertDelimitedMatches(TABLE, TABLE.replace("0.25", "0.26"), tolerance) }
        assertContains(problem, "1 of 6 cells differ (tolerance abs ")
        assertContains(problem, "line 2 (row \"a\"), column 3 (\"value\")")
        assertContains(problem, "off by")
    }

    @Test
    fun relativeScalesWithTheGoldenValue() {
        // 1e-3 relative is 1e-3 at value 1.0 and 1e-6 at value 0.001.
        assertDelimitedMatches(SMALL, SMALL.replace("0.001", "0.0010009"), Tolerance.Relative(1e-3))
        val problem = failure {
            assertDelimitedMatches(SMALL, SMALL.replace("0.001", "0.00101"), Tolerance.Relative(1e-3))
        }
        assertContains(problem, "tolerance rel ")
        assertContains(problem, "off by")
        // The same absolute difference is well inside the same tolerance at a value near one.
        assertDelimitedMatches(TABLE, TABLE.replace("0.5", "0.50001"), Tolerance.Relative(1e-3))
    }

    @Test
    fun upToSignAcceptsAWholeNegatedColumnButNotOneFlippedCell() {
        val flipped = negateColumn(EIGENVECTORS, column = 1)
        assertDelimitedMatches(EIGENVECTORS, flipped, Tolerance.Absolute(1e-9, upToSign = true))
        val problem = failure { assertDelimitedMatches(EIGENVECTORS, flipped, Tolerance.Absolute(1e-9)) }
        assertContains(problem, "3 of 9 cells differ")

        val oneCell = EIGENVECTORS.replace("-0.8", "0.8")
        val stillBad = failure {
            assertDelimitedMatches(EIGENVECTORS, oneCell, Tolerance.Absolute(1e-9, upToSign = true))
        }
        assertContains(stillBad, "cells differ")
    }

    @Test
    fun rowsUnorderedPairsRowsByContentRatherThanByLine() {
        assertDelimitedMatches(TABLE, SHUFFLED, Tolerance.Exact(rowsUnordered = true))
        assertDelimitedMatches(TABLE, SHUFFLED, Tolerance.Absolute(1e-9, rowsUnordered = true))
        val ordered = failure { assertDelimitedMatches(TABLE, SHUFFLED, Tolerance.Absolute(1e-9)) }
        assertContains(ordered, "6 of 6 cells differ")
    }

    @Test
    fun rowsUnorderedStillComparesEveryCell() {
        val rounded = SHUFFLED.replace("0.25", "0.250")
        val exact = failure { assertDelimitedMatches(TABLE, rounded, Tolerance.Exact(rowsUnordered = true)) }
        assertContains(exact, "tolerance exact rows-unordered")
        assertContains(exact, "line 2 (row \"a\"), column 3 (\"value\"): expected \"0.25\", actual \"0.250\"")
        // The same number spelled differently, so a numeric tolerance takes it.
        assertDelimitedMatches(TABLE, rounded, Tolerance.Absolute(1e-9, rowsUnordered = true))
        val renamed = SHUFFLED.replace("a\t", "c\t")
        val relabelled = failure { assertDelimitedMatches(TABLE, renamed, Tolerance.Exact(rowsUnordered = true)) }
        assertContains(relabelled, "6 of 6 cells differ")
    }

    @Test
    fun squareMatricesRejectRowsUnordered() {
        val matrix = SquareMatrix.parse(MATRIX)
        val problem = assertFailsWith<IllegalArgumentException> {
            assertMatrixMatches(matrix, matrix, Tolerance.Absolute(1e-9, rowsUnordered = true))
        }
        assertContains(problem.message.orEmpty(), "rows-unordered has no meaning for a square matrix")
    }

    @Test
    fun formatChecksLabelsAndLayoutButNotNumbers() {
        assertDelimitedMatches(TABLE, TABLE.replace("0.25", "9.9"), Tolerance.Format())
        val renamed = failure { assertDelimitedMatches(TABLE, TABLE.replace("\tvalue", "\tVALUE"), Tolerance.Format()) }
        assertContains(renamed, "header column 3 is \"value\", actual \"VALUE\"")
        val relabelled = failure { assertDelimitedMatches(TABLE, TABLE.replace("a\t", "z\t"), Tolerance.Format()) }
        assertContains(relabelled, "expected the text \"a\", actual \"z\"")
    }

    @Test
    fun textCellsMustBeIdenticalUnderEveryNumericTolerance() {
        val tolerance = Tolerance.Absolute(1.0)
        val changed = failure { assertDelimitedMatches(LABELLED, LABELLED.replace("TBD", "NA"), tolerance) }
        assertContains(changed, "expected the text \"TBD\", actual \"NA\"")
        val becameANumber = failure { assertDelimitedMatches(LABELLED, LABELLED.replace("TBD", "0.0"), tolerance) }
        assertContains(becameANumber, "one is a number and the other is not")
    }

    @Test
    fun reportsStructureBeforeValues() {
        val tolerance = Tolerance.Absolute(1e-9)
        val shortened = failure { assertDelimitedMatches(TABLE, TABLE.substringBeforeLast("b\t2\t0.5\n"), tolerance) }
        assertContains(shortened, "has 2 data rows, actual has 1")
        val narrow = failure { assertDelimitedMatches(TABLE, TABLE.replace("a\t1\t0.25", "a\t1"), tolerance) }
        assertContains(narrow, "has 2 cells, expected 3")
    }

    @Test
    fun showsAtMostFiveDifferencesAndCountsTheRest() {
        val allWrong = EIGENVECTORS.replace("0.", "1.")
        val problem = failure { assertDelimitedMatches(EIGENVECTORS, allWrong, Tolerance.Absolute(1e-9)) }
        assertContains(problem, "6 of 9 cells differ")
        assertContains(problem, "...and 1 more")
    }

    @Test
    fun comparesSquareMatrices() {
        val expected = SquareMatrix.parse(MATRIX)
        assertEquals(listOf("##Matrix_Type=Centered_IBS"), expected.properties)
        assertEquals(listOf("a", "b"), expected.labels)
        assertEquals(1.5, expected.values[0][0])

        val close = SquareMatrix(expected.labels, listOf(doubleArrayOf(1.5, -0.25), doubleArrayOf(-0.25, 1.25000001)))
        assertMatrixMatches(expected, close, Tolerance.Absolute(1e-6))
        val problem = failure { assertMatrixMatches(expected, close, Tolerance.Absolute(1e-12)) }
        assertContains(problem, "row 2 (\"b\"), column 2 (\"b\")")

        val renamed = SquareMatrix(listOf("a", "z"), expected.values)
        val labels = failure { assertMatrixMatches(expected, renamed, Tolerance.Absolute(1e-9)) }
        assertContains(labels, "row 2 is labelled \"b\", actual \"z\"")
    }

    @Test
    fun matrixUpToSignNegatesAWholeColumn() {
        val expected = SquareMatrix.parse(MATRIX)
        val flipped = SquareMatrix(expected.labels, listOf(doubleArrayOf(1.5, 0.25), doubleArrayOf(-0.25, -1.25)))
        assertMatrixMatches(expected, flipped, Tolerance.Absolute(1e-9, upToSign = true))
        val problem = failure { assertMatrixMatches(expected, flipped, Tolerance.Absolute(1e-9)) }
        assertContains(problem, "cells differ")
    }

    @Test
    fun readsTheManifestTolerances() {
        assertEquals(Tolerance.Exact(), goldenTolerance("golden/summary/alleles.txt"))
        assertEquals(Tolerance.Absolute(1e-12), goldenTolerance("golden/summary/sites.full.txt"))
        assertEquals(Tolerance.Absolute(1e-9, upToSign = true), goldenTolerance("golden/pca/eigenvectors.full.txt"))
        assertEquals(Tolerance.Exact(rowsUnordered = true), goldenTolerance("golden/pca/eigenvalues.txt"))
        assertEquals(
            Tolerance.Relative(1e-9, rowsUnordered = true),
            goldenTolerance("golden/pca/eigenvalues.full.txt"),
        )
        assertEquals(Tolerance.Format(), goldenTolerance("golden/pca/eigenvectors.txt"))
    }

    @Test
    fun comparesAGoldenAgainstItself() {
        for ((path, tolerance) in goldenManifest.section("summary").outputs) {
            assertMatchesGolden(path, goldenText(path), tolerance)
        }
        val kinship = "golden/kinship/centered_ibs.txt"
        val matrix = SquareMatrix.parse(goldenText(kinship))
        assertEquals(281, matrix.size)
        assertMatchesGoldenMatrix(kinship, matrix)
        // The rows-unordered goldens have to match themselves line for line as well.
        for (path in goldenManifest.section("pca").outputs.keys.filter { "eigenvalues" in it }) {
            assertMatchesGolden(path, goldenText(path))
        }
    }

    private fun failure(block: () -> Unit): String {
        val message = assertFailsWith<AssertionError>(block = block).message.orEmpty()
        assertTrue(message.isNotEmpty(), "the failure carried no message")
        return message
    }

    private fun negateColumn(table: String, column: Int): String = table.lines().mapIndexed { row, line ->
        if (row == 0 || line.isEmpty()) {
            line
        } else {
            val cells = line.split('\t').toMutableList()
            cells[column] = (-cells[column].toDouble()).toString()
            cells.joinToString("\t")
        }
    }.joinToString("\n")

    private companion object {
        const val TABLE = "name\tcount\tvalue\na\t1\t0.25\nb\t2\t0.5\n"
        const val SHUFFLED = "name\tcount\tvalue\nb\t2\t0.5\na\t1\t0.25\n"
        const val SMALL = "name\tvalue\na\t0.001\n"
        const val LABELLED = "Taxa\tInbreeding Coefficient\n0\tTBD\n"
        const val EIGENVECTORS = "Trait\tPC1\tPC2\nx\t0.5\t0.1\ny\t-0.8\t0.2\nz\t0.3\t0.3\n"
        const val MATRIX = "##Matrix_Type=Centered_IBS\n2\na\t1.5\t-0.25\nb\t-0.25\t1.25\n"
    }
}
