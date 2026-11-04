package net.maizegenetics.tassel.app.chart

import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/** The arithmetic behind the heatmap and the scatter plot, which is all of them that has no pixels in it. */
class DownsampleTest {
    @Test
    fun fewerItemsThanPixelsMeansOneBinEach() {
        val bins = Bins.across(count = 281, pixels = 600)

        assertEquals(281, bins.binCount)
        assertEquals(1, bins.size(0))
        assertEquals(280, bins.start(280))
    }

    @Test
    fun moreItemsThanPixelsSplitsThemEvenlyAndLosesNone() {
        val bins = Bins.across(count = 10, pixels = 3)

        assertEquals(3, bins.binCount)
        assertEquals(listOf(0, 3, 6), List(3) { bins.start(it) })
        assertEquals(listOf(3, 6, 10), List(3) { bins.end(it) })
        // Every item is in exactly one bin, which is what keeps a row of a heatmap from vanishing.
        assertEquals(10, (0 until 3).sumOf { bins.size(it) })
    }

    @Test
    fun noBinIsEverEmpty() {
        for (count in 1..64) {
            for (pixels in 1..64) {
                val bins = Bins.across(count, pixels)
                for (bin in 0 until bins.binCount) {
                    assertTrue(bins.size(bin) >= 1, "$count items in ${bins.binCount} bins, bin $bin")
                }
            }
        }
    }

    @Test
    fun aCanvasWithNoPixelsStillBinsSomething() {
        assertEquals(1, Bins.across(count = 281, pixels = 0).binCount)
        assertFailsWith<IllegalArgumentException> { Bins.across(count = 0, pixels = 10) }
    }

    @Test
    fun aHeatmapCellIsTheMeanOfTheBlockItCovers() {
        // A 4 x 4 matrix whose value is row * 10 + column, binned 2 x 2.
        val cells = heatmapCells(Bins(count = 4, binCount = 2)) { row, column -> (row * 10 + column).toFloat() }

        // Top-left block covers (0,0) (0,1) (1,0) (1,1): 0, 1, 10, 11.
        assertEquals(5.5f, cells[0])
        // Top-right covers 2, 3, 12, 13.
        assertEquals(7.5f, cells[1])
        assertEquals(25.5f, cells[2])
        assertEquals(27.5f, cells[3])
    }

    @Test
    fun aHeatmapWithOneBinPerRowIsTheMatrixItself() {
        val cells = heatmapCells(Bins(count = 3, binCount = 3)) { row, column -> (row - column).toFloat() }

        assertContentEquals(floatArrayOf(0f, -1f, -2f, 1f, 0f, -1f, 2f, 1f, 0f), cells)
    }

    @Test
    fun averagingRatherThanSamplingKeepsAnOutlierVisible() {
        // One high value in a block of zeros still moves the block's color; sampling would miss it entirely.
        val cells = heatmapCells(Bins(count = 4, binCount = 2)) { row, column ->
            if (row == 1 && column == 1) 100f else 0f
        }

        assertEquals(25f, cells[0])
    }

    @Test
    fun theRangeOfAMatrixIsItsLowestAndHighestValue() {
        assertEquals(-2f to 7f, range(floatArrayOf(3f, -2f, 7f, 0f)))
        assertEquals(0f to 1f, range(FloatArray(0)))
    }

    @Test
    fun aFlatRangePaintsOneMidColorRatherThanDividingByZero() {
        assertEquals(0.5f, fraction(2f, 2f, 2f))
        assertEquals(0f, fraction(0f, 0f, 10f))
        assertEquals(1f, fraction(10f, 0f, 10f))
        assertEquals(0.5f, fraction(5f, 0f, 10f))
        // Values outside the range clamp rather than painting off the end of the gradient.
        assertEquals(1f, fraction(99f, 0f, 10f))
    }

    @Test
    fun aScatterWithRoomForEveryPointKeepsEveryPoint() {
        val x = FloatArray(281) { it.toFloat() }
        val y = FloatArray(281) { (it * it).toFloat() }

        val kept = scatterSample(x, y, width = 600f, height = 400f)

        assertEquals(281, kept.size)
        assertContentEquals(IntArray(281) { it }, kept)
    }

    @Test
    fun pointsThatWouldLandOnTheSamePixelAreDrawnOnce() {
        // Four points in two clusters; a plot two cells wide can only show one of each cluster.
        val x = floatArrayOf(0f, 0.001f, 1f, 1.001f)
        val y = floatArrayOf(0f, 0.001f, 0f, 0.001f)

        val kept = scatterSample(x, y, width = 6f, height = 3f, cellSize = 3f)

        assertEquals(2, kept.size)
        // Index order decides which survives, so the picture does not depend on iteration order.
        assertContentEquals(intArrayOf(0, 2), kept)
    }

    @Test
    fun aScatterWithNothingToDrawDrawsNothing() {
        assertEquals(0, scatterSample(FloatArray(0), FloatArray(0), 100f, 100f).size)
        assertEquals(0, scatterSample(floatArrayOf(1f), floatArrayOf(1f), 0f, 100f).size)
    }

    @Test
    fun aScaleLabelKeepsItsSign() {
        assertEquals("0.500", scaleLabel(0.5f))
        assertEquals("-0.250", scaleLabel(-0.25f))
        assertEquals("1.000", scaleLabel(1f))
        assertEquals("2.125", scaleLabel(2.125f))
    }
}
