package net.maizegenetics.tassel.app.chart

import androidx.compose.foundation.Canvas
import androidx.compose.foundation.layout.fillMaxSize
import androidx.compose.material3.MaterialTheme
import androidx.compose.runtime.Composable
import androidx.compose.runtime.LaunchedEffect
import androidx.compose.runtime.getValue
import androidx.compose.runtime.mutableIntStateOf
import androidx.compose.runtime.mutableStateOf
import androidx.compose.runtime.remember
import androidx.compose.runtime.setValue
import androidx.compose.ui.Modifier
import androidx.compose.ui.geometry.Offset
import androidx.compose.ui.geometry.Size
import androidx.compose.ui.graphics.Canvas
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.graphics.FilterQuality
import androidx.compose.ui.graphics.ImageBitmap
import androidx.compose.ui.graphics.Paint
import androidx.compose.ui.graphics.drawscope.DrawScope
import androidx.compose.ui.graphics.lerp
import androidx.compose.ui.layout.onSizeChanged
import androidx.compose.ui.text.TextMeasurer
import androidx.compose.ui.text.TextStyle
import androidx.compose.ui.text.drawText
import androidx.compose.ui.text.rememberTextMeasurer
import androidx.compose.ui.unit.IntSize
import androidx.compose.ui.unit.dp
import net.maizegenetics.tassel.app.analysis.Chart
import net.maizegenetics.tassel.app.yieldToEventLoop
import kotlin.math.min

/** Room left under the square for the color scale. */
private val SCALE_HEIGHT = 28.dp

/**
 * Rows of the heatmap bitmap painted between yields.
 *
 * One rect per taxon pair on the tutorial matrix is about 79,000 draws, and on wasm that is one task of
 * about 300 ms, which is the frame gap issue 70 measures. A slice stays well under 100 ms, and the bitmap
 * keeps the earlier rows so the last slice does not redraw the whole square.
 */
private const val HEATMAP_ROWS_PER_SLICE = 24

/**
 * A symmetric matrix as a heatmap, drawn on one [Canvas].
 *
 * There are no composables per matrix cell and none per drawn cell: the whole picture is one canvas, as
 * `.cursor/rules/compose-ui.mdc` asks and as `GenotypeGrid` already does for genotypes. When the matrix has more
 * rows than the square has pixels the cells are averaged down to one per pixel ([heatmapCells]), once per matrix
 * and size rather than once per frame. The bitmap is filled a few rows at a time, yielding between slices, so
 * that fill does not freeze the page.
 */
@Composable
fun KinshipHeatmap(heatmap: Chart.Heatmap, modifier: Modifier = Modifier) {
    val scheme = MaterialTheme.colorScheme
    val labelStyle = MaterialTheme.typography.labelSmall.copy(color = scheme.onSurfaceVariant)
    val measurer = rememberTextMeasurer(cacheSize = 8)
    val from = scheme.surfaceVariant
    val to = scheme.primary

    // The layout phase is where the size can be read without writing snapshot state from inside the draw phase.
    var viewport by remember(heatmap) { mutableStateOf(IntSize.Zero) }
    val bins = remember(heatmap, viewport) { Bins.across(heatmap.size, min(viewport.width, viewport.height)) }
    val cells = remember(heatmap, bins) { heatmapCells(bins) { row, column -> heatmap.value(row, column) } }
    val (low, high) = remember(cells) { range(cells) }
    val bitmap = remember(cells) { ImageBitmap(bins.binCount, bins.binCount) }
    var paintedRows by remember(cells) { mutableIntStateOf(0) }

    LaunchedEffect(cells, from, to) {
        val canvas = Canvas(bitmap)
        val paint = Paint()
        val count = bins.binCount
        var row = 0
        while (row < count) {
            val end = min(count, row + HEATMAP_ROWS_PER_SLICE)
            for (rasterRow in row until end) {
                for (column in 0 until count) {
                    paint.color = lerp(from, to, fraction(cells[rasterRow * count + column], low, high))
                    val x = column.toFloat()
                    val y = rasterRow.toFloat()
                    canvas.drawRect(x, y, x + 1f, y + 1f, paint)
                }
            }
            row = end
            paintedRows = row
            if (row < count) yieldToEventLoop()
        }
    }

    Canvas(modifier = modifier.fillMaxSize().onSizeChanged { viewport = it }) {
        val scaleHeight = SCALE_HEIGHT.toPx()
        val square = min(size.width, size.height - scaleHeight)
        if (square <= 0f || paintedRows == 0) return@Canvas
        val side = square.toInt().coerceAtLeast(1)
        drawImage(
            image = bitmap,
            dstSize = IntSize(side, side),
            filterQuality = FilterQuality.None,
        )
        if (paintedRows == bins.binCount) {
            drawScale(measurer, labelStyle, from, to, low, high, square, scaleHeight)
        }
    }
}

/** The legend: the same gradient as a bar, with the lowest and highest value of the matrix at its ends. */
@Suppress("LongParameterList")
private fun DrawScope.drawScale(
    measurer: TextMeasurer,
    style: TextStyle,
    from: Color,
    to: Color,
    low: Float,
    high: Float,
    width: Float,
    height: Float,
) {
    val barHeight = height / 3f
    val top = size.height - height
    val steps = width.toInt().coerceAtLeast(1)
    val step = width / steps
    for (index in 0 until steps) {
        drawRect(
            color = lerp(from, to, index.toFloat() / steps),
            topLeft = Offset(index * step, top),
            size = Size(step + 1f, barHeight),
        )
    }
    drawText(measurer, scaleLabel(low), Offset(0f, top + barHeight), style, maxLines = 1)
    val highest = scaleLabel(high)
    val measured = measurer.measure(highest, style, maxLines = 1)
    drawText(measurer, highest, Offset(width - measured.size.width, top + barHeight), style, maxLines = 1)
}

/**
 * A matrix value with three decimals, which is enough to tell two kinship coefficients apart on a legend.
 * `tassel-core`'s formatters spell numbers the way TASSEL 5 writes files, which this is not.
 */
internal fun scaleLabel(value: Float): String {
    val scaled = (value * 1000f).toInt()
    val sign = if (scaled < 0) "-" else ""
    val magnitude = if (scaled < 0) -scaled else scaled
    return "$sign${magnitude / 1000}.${(magnitude % 1000).toString().padStart(3, '0')}"
}
