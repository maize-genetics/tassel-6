package net.maizegenetics.tassel.math

import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class TriangularTest {
    // Upper triangle holds U; the strict lower triangle is garbage that every solve must ignore.
    private val upper = DenseMatrix.ofRows(
        doubleArrayOf(2.0, -1.0, 3.0, 99.0),
        doubleArrayOf(77.0, 4.0, 0.5, 99.0),
        doubleArrayOf(77.0, 77.0, -5.0, 99.0),
    )

    // Lower triangle holds L; the strict upper triangle is garbage.
    private val lower = DenseMatrix.ofRows(
        doubleArrayOf(2.0, 77.0, 77.0),
        doubleArrayOf(-1.0, 4.0, 77.0),
        doubleArrayOf(3.0, 0.5, -5.0),
        doubleArrayOf(99.0, 99.0, 99.0),
    )

    private val x = doubleArrayOf(1.0, -2.0, 0.75)

    @Test
    fun solveUpper() {
        val b = DoubleArray(3) { i -> (i until 3).sumOf { upper[i, it] * x[it] } }
        Triangular.solveUpper(upper, 3, b)
        assertClose(x, b)
    }

    @Test
    fun solveUpperTransposed() {
        val b = DoubleArray(3) { j -> (0..j).sumOf { upper[it, j] * x[it] } }
        Triangular.solveUpperTransposed(upper, 3, b)
        assertClose(x, b)
    }

    @Test
    fun solveLower() {
        val b = DoubleArray(3) { i -> (0..i).sumOf { lower[i, it] * x[it] } }
        Triangular.solveLower(lower, 3, b)
        assertClose(x, b)
    }

    @Test
    fun solveLowerTransposed() {
        val b = DoubleArray(3) { j -> (j until 3).sumOf { lower[it, j] * x[it] } }
        Triangular.solveLowerTransposed(lower, 3, b)
        assertClose(x, b)
    }

    @Test
    fun solvesALeadingBlockAndLeavesTheRestAlone() {
        val b = doubleArrayOf(2.0, 8.0, 123.0)
        Triangular.solveUpper(upper, 2, b)
        // Directly: x1 = 8 / 4 = 2, x0 = (2 + 1 * 2) / 2 = 2.
        assertTrue(abs(b[0] - 2.0) < 1e-15 && abs(b[1] - 2.0) < 1e-15 && b[2] == 123.0)
    }

    @Test
    fun rejectsBlocksThatDoNotFit() {
        assertFailsWith<IllegalArgumentException> { Triangular.solveUpper(upper, 4, DoubleArray(4)) }
        assertFailsWith<IllegalArgumentException> { Triangular.solveUpper(upper, 3, DoubleArray(2)) }
    }

    private fun assertClose(expected: DoubleArray, actual: DoubleArray) {
        for (i in expected.indices) {
            assertTrue(
                abs(expected[i] - actual[i]) < 1e-14,
                "index $i: ${expected[i]} vs ${actual[i]}",
            )
        }
    }
}
