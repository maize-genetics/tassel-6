package net.maizegenetics.tassel.table

import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import kotlin.test.fail

/**
 * The oracle for the [DoubleFormat] port (ADR-0009).
 *
 * Every table golden was written twice by the same TASSEL 5 run: `name.txt` through TASSEL's own
 * `DoubleFormat`, and `name.full.txt` with that one class replaced so doubles are written with
 * `Double.toString`. Formatting each cell of the full-precision file must therefore reproduce the
 * rounded file exactly, cell for cell. No genotype reader is involved.
 */
class DoubleFormatGoldenTest {
    @Test
    fun reproducesEveryRoundedSummaryAndPcaTable() {
        val checked = PAIRS.sumOf { checkPair(it) }
        // Guards against a silently empty run if the goldens ever move.
        assertTrue(checked > 100_000, "only $checked cells were checked")
    }

    /**
     * `pca/pcs.txt` is written by TASSEL's `Phenotype` exporter, not by `TableReportUtils`, so it
     * never goes through `DoubleFormat`: its rounded and full-precision files are the same bytes.
     * Issue #19 owns that format.
     */
    @Test
    fun principalComponentScoresAreNotFormattedByDoubleFormat() {
        assertEquals(goldenText("golden/pca/pcs.full.txt"), goldenText("golden/pca/pcs.txt"))
    }

    /** Formats each cell of `$name.full.txt` and compares it to `$name.txt`. Returns cells checked. */
    private fun checkPair(name: String): Int {
        val full = rows(goldenText("golden/$name.full.txt"))
        val rounded = rows(goldenText("golden/$name.txt"))
        assertEquals(rounded.size, full.size, "$name: the two files have different row counts")
        var cells = 0
        full.forEachIndexed { row, fullCells ->
            val roundedCells = rounded[row]
            assertEquals(
                roundedCells.size,
                fullCells.size,
                "$name line ${row + 1}: the two files have different cell counts",
            )
            fullCells.forEachIndexed { column, cell ->
                cells++
                val expected = roundedCells[column]
                // TASSEL writes a Double cell with Double.toString, which always has a decimal
                // point; anything else is an Int, a Long, or text, and is copied through unchanged.
                val value = if ('.' in cell) cell.toDoubleOrNull() else null
                val actual = if (value == null) cell else DoubleFormat.format(value)
                if (actual != expected) {
                    fail(
                        "golden/$name.txt line ${row + 1} (row \"${fullCells.firstOrNull()}\"), " +
                            "column ${column + 1}: DoubleFormat.format($cell) is \"$actual\", " +
                            "the golden has \"$expected\"",
                    )
                }
            }
        }
        return cells
    }

    private fun rows(text: String): List<List<String>> =
        normalizeLineEndings(text).lines().dropLastWhile { it.isEmpty() }.map { it.split('\t') }

    private companion object {
        /** Every `*.full.txt` table golden that TASSEL wrote through `TableReportUtils`. */
        val PAIRS = listOf(
            "summary/overview",
            "summary/alleles",
            "summary/sites",
            "summary/taxa",
            "pca/eigenvalues",
            "pca/eigenvectors",
        )
    }
}
