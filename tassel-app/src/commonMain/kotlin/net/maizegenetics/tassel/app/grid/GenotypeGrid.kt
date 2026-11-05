package net.maizegenetics.tassel.app.grid

import androidx.compose.foundation.Canvas
import androidx.compose.foundation.gestures.Orientation
import androidx.compose.foundation.gestures.rememberScrollableState
import androidx.compose.foundation.gestures.scrollable
import androidx.compose.foundation.layout.fillMaxSize
import androidx.compose.material3.MaterialTheme
import androidx.compose.runtime.Composable
import androidx.compose.runtime.getValue
import androidx.compose.runtime.mutableStateOf
import androidx.compose.runtime.remember
import androidx.compose.runtime.setValue
import androidx.compose.ui.Modifier
import androidx.compose.ui.geometry.Offset
import androidx.compose.ui.geometry.Size
import androidx.compose.ui.graphics.drawscope.DrawScope
import androidx.compose.ui.graphics.drawscope.clipRect
import androidx.compose.ui.layout.onSizeChanged
import androidx.compose.ui.platform.LocalDensity
import androidx.compose.ui.text.TextMeasurer
import androidx.compose.ui.text.TextStyle
import androidx.compose.ui.text.drawText
import androidx.compose.ui.text.rememberTextMeasurer
import androidx.compose.ui.text.style.TextOverflow
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.toSize
import net.maizegenetics.tassel.app.theme.AlleleColors
import net.maizegenetics.tassel.app.theme.LocalAlleleColors
import net.maizegenetics.tassel.genotype.GenotypeTable

/** Width of a site column. */
private val CELL_WIDTH = 10.dp

/** Height of a taxon row, which also sets the size of the taxon name text. */
private val CELL_HEIGHT = 14.dp

/** Width of the left gutter holding taxon names. */
private val GUTTER_WIDTH = 150.dp

/** Height of the top gutter holding site positions. */
private val HEADER_HEIGHT = 18.dp

/** Site labels are drawn every this many columns, so they do not collide at [CELL_WIDTH]. */
private const val LABEL_EVERY = 10

/**
 * The genotype grid: taxa down, sites across, each call a rectangle colored by its alleles.
 *
 * Nothing here composes a node per cell. The whole cell layer is one [Canvas], and each frame reads only the window
 * the viewport shows (see [gridWindow] and [readWindow]), so the cost of a frame depends on the size of the window
 * and not on the size of [table]. Taxon names and site positions are drawn onto the same canvas, which keeps them in
 * step with the scroll offset for free.
 */
@Composable
fun GenotypeGrid(table: GenotypeTable, modifier: Modifier = Modifier) {
    val colors = LocalAlleleColors.current
    val scheme = MaterialTheme.colorScheme
    val labelStyle = MaterialTheme.typography.labelSmall.copy(color = scheme.onSurfaceVariant)
    // One entry per visible row and site label, so a frame measures no text it measured in the frame before.
    val textMeasurer = rememberTextMeasurer(cacheSize = 128)
    val density = LocalDensity.current

    var offset by remember(table) { mutableStateOf(Offset.Zero) }
    var viewport by remember(table) { mutableStateOf(Size.Zero) }

    val cellWidth = with(density) { CELL_WIDTH.toPx() }
    val cellHeight = with(density) { CELL_HEIGHT.toPx() }
    val gutterWidth = with(density) { GUTTER_WIDTH.toPx() }
    val headerHeight = with(density) { HEADER_HEIGHT.toPx() }

    /** Scrolls one axis and returns how much of [delta] was used, so an unscrollable axis passes the gesture on. */
    fun scroll(horizontal: Boolean, delta: Float): Float {
        val cells = if (horizontal) table.numberOfSites else table.numberOfTaxa
        val cellSize = if (horizontal) cellWidth else cellHeight
        val visible = if (horizontal) viewport.width - gutterWidth else viewport.height - headerHeight
        val limit = (cells * cellSize - visible).coerceAtLeast(0f)
        val before = if (horizontal) offset.x else offset.y
        val after = (before - delta).coerceIn(0f, limit)
        offset = if (horizontal) offset.copy(x = after) else offset.copy(y = after)
        return before - after
    }

    val horizontal = rememberScrollableState { scroll(horizontal = true, delta = it) }
    val vertical = rememberScrollableState { scroll(horizontal = false, delta = it) }

    Canvas(
        modifier = modifier
            .fillMaxSize()
            // The scroll limits need the viewport size, and the layout phase is where it can be read without
            // writing snapshot state from inside the draw phase.
            .onSizeChanged { viewport = it.toSize() }
            .scrollable(vertical, Orientation.Vertical)
            .scrollable(horizontal, Orientation.Horizontal),
    ) {
        val cells = Size(size.width - gutterWidth, size.height - headerHeight)
        val window = gridWindow(
            numberOfTaxa = table.numberOfTaxa,
            numberOfSites = table.numberOfSites,
            taxonOffset = offset.y,
            siteOffset = offset.x,
            width = cells.width,
            height = cells.height,
            cellWidth = cellWidth,
            cellHeight = cellHeight,
        )
        if (window.cellCount == 0) return@Canvas

        val genotypes = table.readWindow(window)
        // Where table cell (0, 0) would sit, so a cell's corner is its index times its size plus this.
        val origin = Offset(gutterWidth - offset.x, headerHeight - offset.y)

        clipRect(left = gutterWidth, top = headerHeight, right = size.width, bottom = size.height) {
            drawCells(genotypes, colors, origin, cellWidth, cellHeight)
        }
        clipRect(left = 0f, top = headerHeight, right = gutterWidth, bottom = size.height) {
            drawTaxonNames(table, window, textMeasurer, labelStyle, origin.y, cellHeight, gutterWidth)
        }
        clipRect(left = gutterWidth, top = 0f, right = size.width, bottom = headerHeight) {
            drawSiteLabels(table, window, textMeasurer, labelStyle, origin.x, cellWidth, headerHeight)
        }
        drawRect(color = scheme.outlineVariant, topLeft = Offset(gutterWidth, 0f), size = Size(1f, size.height))
        drawRect(color = scheme.outlineVariant, topLeft = Offset(0f, headerHeight), size = Size(size.width, 1f))
    }
}

/**
 * Paints the calls of [genotypes]. Consecutive rows of a column that share a call become one rectangle, so a column
 * of 50 visible taxa in a dataset with little variation costs a handful of draw calls instead of 50.
 */
private fun DrawScope.drawCells(
    genotypes: WindowGenotypes,
    colors: AlleleColors,
    origin: Offset,
    cellWidth: Float,
    cellHeight: Float,
) {
    val window = genotypes.window
    for (column in 0 until window.siteCount) {
        val x = origin.x + (window.firstSite + column) * cellWidth
        for (run in genotypes.columnRuns(column)) {
            drawRect(
                color = colors.of(run.genotype),
                topLeft = Offset(x, origin.y + (window.firstTaxon + run.firstRow) * cellHeight),
                size = Size(cellWidth, run.rowCount * cellHeight),
            )
        }
    }
}

/** Writes the name of each visible taxon in the left gutter, aligned with its row. */
private fun DrawScope.drawTaxonNames(
    table: GenotypeTable,
    window: GridWindow,
    measurer: TextMeasurer,
    style: TextStyle,
    originY: Float,
    cellHeight: Float,
    gutterWidth: Float,
) {
    for (taxon in window.taxa) {
        drawText(
            textMeasurer = measurer,
            text = table.taxa[taxon].name,
            topLeft = Offset(TEXT_PADDING, originY + taxon * cellHeight),
            style = style,
            maxLines = 1,
            overflow = TextOverflow.Ellipsis,
            size = Size(gutterWidth - 2 * TEXT_PADDING, cellHeight),
        )
    }
}

/** Writes the physical position of every [LABEL_EVERY]th visible site above its column. */
private fun DrawScope.drawSiteLabels(
    table: GenotypeTable,
    window: GridWindow,
    measurer: TextMeasurer,
    style: TextStyle,
    originX: Float,
    cellWidth: Float,
    headerHeight: Float,
) {
    val first = window.firstSite + (LABEL_EVERY - window.firstSite % LABEL_EVERY) % LABEL_EVERY
    for (site in first..window.sites.last step LABEL_EVERY) {
        val position = table.positions[site]
        drawText(
            textMeasurer = measurer,
            text = "${position.chromosome.name}:${position.position}",
            topLeft = Offset(originX + site * cellWidth + TEXT_PADDING, 0f),
            style = style,
            maxLines = 1,
            overflow = TextOverflow.Clip,
            size = Size(LABEL_EVERY * cellWidth, headerHeight),
        )
    }
}

/** Gap between text and the edge of the gutter it sits in, in pixels. */
private const val TEXT_PADDING = 4f
