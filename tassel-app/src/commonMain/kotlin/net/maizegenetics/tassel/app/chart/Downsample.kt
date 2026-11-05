package net.maizegenetics.tassel.app.chart

import kotlin.math.floor

/**
 * An even split of [count] items into [binCount] bins, which is how a chart with fewer pixels than data decides
 * what each pixel stands for.
 *
 * Bins differ in size by at most one and none is empty, so no row of a heatmap is dropped and none is drawn twice.
 * `start(bin) = bin * count / binCount` is the same split `IntRange.chunked` would make if chunks could be uneven.
 */
class Bins(val count: Int, val binCount: Int) {
    init {
        require(count > 0) { "nothing to bin: $count items" }
        require(binCount in 1..count) { "cannot put $count items in $binCount bins" }
    }

    /** First item index of [bin]. */
    fun start(bin: Int): Int {
        require(bin in 0 until binCount) { "bin $bin of $binCount" }
        return (bin.toLong() * count / binCount).toInt()
    }

    /** One past the last item index of [bin]. */
    fun end(bin: Int): Int = if (bin == binCount - 1) count else start(bin + 1)

    /** How many items [bin] covers; always at least one. */
    fun size(bin: Int): Int = end(bin) - start(bin)

    companion object {
        /**
         * Bins for [count] items drawn across [pixels] pixels, at least one pixel per bin. Fewer items than pixels
         * means one bin each and no downsampling at all, which is the usual case for 281 taxa on a laptop.
         */
        fun across(count: Int, pixels: Int): Bins = Bins(count, count.coerceAtMost(pixels.coerceAtLeast(1)))
    }
}

/**
 * The mean of every block of a symmetric [size] by [size] matrix under [rows] by [rows] bins, row-major.
 *
 * Averaging rather than sampling is what keeps a downsampled kinship heatmap honest: one taxon per 10 drawn would
 * hide a related pair entirely, while the mean of the block it falls in still shifts the color. The result has
 * `rows.binCount * rows.binCount` entries and is computed once per matrix and size, not once per frame.
 */
fun heatmapCells(rows: Bins, value: (Int, Int) -> Float): FloatArray {
    val bins = rows.binCount
    val cells = FloatArray(bins * bins)
    for (row in 0 until bins) {
        for (column in 0 until bins) {
            var total = 0.0
            for (i in rows.start(row) until rows.end(row)) {
                for (j in rows.start(column) until rows.end(column)) total += value(i, j)
            }
            cells[row * bins + column] = (total / (rows.size(row) * rows.size(column))).toFloat()
        }
    }
    return cells
}

/** The smallest and largest value in [cells], as `(low, high)`; `(0, 1)` for an empty array. */
fun range(cells: FloatArray): Pair<Float, Float> {
    if (cells.isEmpty()) return 0f to 1f
    var low = cells[0]
    var high = cells[0]
    for (cell in cells) {
        if (cell < low) low = cell
        if (cell > high) high = cell
    }
    return low to high
}

/**
 * Where a value sits in `[low, high]`, as a fraction in `[0, 1]`. A flat range gives 0.5 rather than dividing by
 * zero, so a matrix of identical values paints one mid color instead of NaN.
 */
fun fraction(value: Float, low: Float, high: Float): Float =
    if (high <= low) 0.5f else ((value - low) / (high - low)).coerceIn(0f, 1f)

/**
 * The indices of [x] to draw, keeping at most one point per [cellSize]-pixel square of a [width] by [height] plot.
 *
 * A scatter plot cannot show two points that land on the same pixel, so drawing both is work with no picture to
 * show for it. Points are kept in index order, so which one survives a collision does not depend on iteration
 * order, and a plot with fewer points than cells keeps all of them.
 */
fun scatterSample(
    x: FloatArray,
    y: FloatArray,
    width: Float,
    height: Float,
    cellSize: Float = DEFAULT_SCATTER_CELL,
): IntArray {
    require(x.size == y.size) { "${x.size} x values and ${y.size} y values" }
    if (x.isEmpty() || width <= 0f || height <= 0f) return IntArray(0)
    val (lowX, highX) = range(x)
    val (lowY, highY) = range(y)
    val columns = (width / cellSize).toInt().coerceAtLeast(1)
    val rows = (height / cellSize).toInt().coerceAtLeast(1)
    if (columns.toLong() * rows >= x.size) return IntArray(x.size) { it }
    val taken = HashSet<Int>(x.size)
    val kept = ArrayList<Int>(x.size)
    for (index in x.indices) {
        val column = bucket(fraction(x[index], lowX, highX), columns)
        val row = bucket(fraction(y[index], lowY, highY), rows)
        if (taken.add(row * columns + column)) kept += index
    }
    return kept.toIntArray()
}

/** Which of [buckets] a fraction in `[0, 1]` falls in; 1.0 belongs to the last one rather than one past it. */
private fun bucket(fraction: Float, buckets: Int): Int = floor(fraction * buckets).toInt().coerceIn(0, buckets - 1)

/** Side of the square a scatter plot keeps one point in, in pixels. Below this, points overlap anyway. */
const val DEFAULT_SCATTER_CELL: Float = 3f
