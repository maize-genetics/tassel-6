package net.maizegenetics.tassel.app.chart

import androidx.compose.foundation.Canvas
import androidx.compose.foundation.layout.fillMaxSize
import androidx.compose.material3.MaterialTheme
import androidx.compose.runtime.Composable
import androidx.compose.runtime.remember
import androidx.compose.ui.Modifier
import androidx.compose.ui.geometry.Offset
import androidx.compose.ui.geometry.Size
import androidx.compose.ui.graphics.drawscope.Stroke
import androidx.compose.ui.text.drawText
import androidx.compose.ui.text.rememberTextMeasurer
import androidx.compose.ui.unit.dp
import net.maizegenetics.tassel.app.analysis.Chart

/** Room left on the left and at the bottom for the axis labels. */
private val AXIS_MARGIN = 36.dp

/** Radius of a plotted point, in pixels. */
private const val POINT_RADIUS = 2.5f

/** How far above the axis the lowest y label sits, in pixels. */
private const val LABEL_LIFT = 12f

/**
 * Two principal components as a scatter plot, drawn on one [Canvas].
 *
 * One canvas, no composable per point. Points that would land on the same few pixels are dropped before drawing
 * ([scatterSample]), so a plot of a million taxa draws as many points as the plot has room for and no more; for
 * the 281 taxa of the tutorial dataset nothing is dropped.
 */
@Composable
fun PcaScatter(scatter: Chart.Scatter, modifier: Modifier = Modifier) {
    val scheme = MaterialTheme.colorScheme
    val labelStyle = MaterialTheme.typography.labelSmall.copy(color = scheme.onSurfaceVariant)
    val measurer = rememberTextMeasurer(cacheSize = 8)
    val xRange = remember(scatter) { range(scatter.x) }
    val yRange = remember(scatter) { range(scatter.y) }

    Canvas(modifier = modifier.fillMaxSize()) {
        val margin = AXIS_MARGIN.toPx()
        val plot = Size(size.width - margin, size.height - margin)
        if (plot.width <= 0f || plot.height <= 0f) return@Canvas

        drawRect(
            color = scheme.outlineVariant,
            topLeft = Offset(margin, 0f),
            size = Size(plot.width, plot.height),
            style = Stroke(width = 1f),
        )
        for (index in scatterSample(scatter.x, scatter.y, plot.width, plot.height)) {
            drawCircle(
                color = scheme.primary,
                radius = POINT_RADIUS,
                center = Offset(
                    margin + fraction(scatter.x[index], xRange.first, xRange.second) * plot.width,
                    // Screen y grows downward, so the largest value belongs at the top.
                    (1f - fraction(scatter.y[index], yRange.first, yRange.second)) * plot.height,
                ),
            )
        }
        drawText(measurer, scatter.xLabel, Offset(margin, plot.height + 4f), labelStyle, maxLines = 1)
        drawText(measurer, scatter.yLabel, Offset(0f, 0f), labelStyle, maxLines = 1)
        drawText(measurer, scaleLabel(yRange.first), Offset(0f, plot.height - LABEL_LIFT), labelStyle, maxLines = 1)
    }
}
