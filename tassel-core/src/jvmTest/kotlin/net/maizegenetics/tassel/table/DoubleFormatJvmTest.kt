package net.maizegenetics.tassel.table

import java.io.File
import java.text.DecimalFormat
import java.text.NumberFormat
import java.util.Locale
import kotlin.math.pow
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Cross-checks the [DoubleFormat] port against the JVM classes TASSEL 5 formats with, over far more
 * values than a table of examples can hold. This is JVM only because `DecimalFormat` is; the same
 * port is checked against the TASSEL 5 goldens on every target in `commonTest`.
 */
class DoubleFormatJvmTest {
    @Test
    fun matchesTassel5OverRandomValues() {
        val random = Random(20261104)
        repeat(200_000) {
            val value = when (it % 6) {
                0 -> random.nextDouble(-1.0, 1.0)
                1 -> random.nextDouble(0.0009, 0.0011) // the plain to scientific boundary
                2 -> random.nextDouble(999_990.0, 1_000_010.0) // the other boundary
                3 -> random.nextDouble() * 10.0.pow(random.nextInt(-300, 300))
                // Whole numbers, which TASSEL 5 counts with and which round ties up.
                4 -> random.nextLong(-1_000_000_000L, 1_000_000_000L).toDouble()
                else -> random.nextLong(0L, 1L shl 53).toDouble()
            }
            assertEquals(tassel5(value), DoubleFormat.format(value), "value $value")
        }
    }

    @Test
    fun matchesTassel5OnEveryValueInTheTableGoldens() {
        val root = System.getenv("TASSEL_TESTDATA_DIR") ?: error("TASSEL_TESTDATA_DIR is not set")
        var checked = 0
        for (path in GOLDEN_TABLES) {
            for (line in File("$root/golden/$path.full.txt").readLines()) {
                for (cell in line.split('\t')) {
                    val value = if ('.' in cell) cell.toDoubleOrNull() else null
                    if (value != null) {
                        checked++
                        assertEquals(tassel5(value), DoubleFormat.format(value), "$path cell $cell")
                    }
                }
            }
        }
        assertTrue(checked > 50_000, "only $checked values were checked")
    }

    @Test
    fun differsFromTassel5OnlyForTinySubnormals() {
        // DecimalFormat rounds the shortest decimal that round-trips, not the exact value. Below
        // about 5e-319 an ulp is wider than the fifth digit, so that decimal has fewer than five
        // digits and the two disagree. No analysis produces such a value; see DoubleFormat's KDoc.
        assertEquals("4.9E-324", tassel5(Double.MIN_VALUE))
        assertEquals("4.9407E-324", DoubleFormat.format(Double.MIN_VALUE))
        val smallestNormal = java.lang.Double.MIN_NORMAL
        assertEquals(tassel5(smallestNormal), DoubleFormat.format(smallestNormal))
    }

    /** TASSEL 5 `net.maizegenetics.util.DoubleFormat`, transcribed. */
    private fun tassel5(value: Double): String {
        if (value.isNaN()) return "NaN"
        if (value.isInfinite()) return "Infinity"
        if (value == 0.0) return plain.format(value)
        if (value < .001 || value >= 1000000) return scientific.format(value)
        return plain.format(value)
    }

    private companion object {
        val GOLDEN_TABLES = listOf(
            "summary/overview",
            "summary/alleles",
            "summary/sites",
            "summary/taxa",
            "pca/eigenvalues",
            "pca/eigenvectors",
        )

        val plain: NumberFormat =
            NumberFormat.getInstance(Locale.ROOT).also { (it as DecimalFormat).applyPattern("0.#####") }
        val scientific: NumberFormat =
            NumberFormat.getInstance(Locale.ROOT).also { (it as DecimalFormat).applyPattern("0.0###E0") }
    }
}
