package net.maizegenetics.tassel.app.ui

import androidx.compose.foundation.horizontalScroll
import androidx.compose.foundation.layout.Arrangement
import androidx.compose.foundation.layout.Column
import androidx.compose.foundation.layout.Row
import androidx.compose.foundation.layout.fillMaxSize
import androidx.compose.foundation.layout.fillMaxWidth
import androidx.compose.foundation.layout.padding
import androidx.compose.foundation.layout.width
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.rememberScrollState
import androidx.compose.material3.HorizontalDivider
import androidx.compose.material3.MaterialTheme
import androidx.compose.material3.OutlinedButton
import androidx.compose.material3.PrimaryTabRow
import androidx.compose.material3.Tab
import androidx.compose.material3.Text
import androidx.compose.runtime.Composable
import androidx.compose.ui.Modifier
import androidx.compose.ui.text.style.TextOverflow
import androidx.compose.ui.unit.dp
import net.maizegenetics.tassel.app.analysis.Chart
import net.maizegenetics.tassel.app.analysis.TableView
import net.maizegenetics.tassel.app.chart.KinshipHeatmap
import net.maizegenetics.tassel.app.chart.PcaScatter
import net.maizegenetics.tassel.app.state.AnalysisResult
import net.maizegenetics.tassel.app.state.AppEvent

/** Width of a table column. Every column is the same width, so a row needs no measuring pass. */
private val COLUMN_WIDTH = 120.dp

/**
 * Columns drawn at once. A kinship matrix has one column per taxon, so without a cap a row of the tutorial
 * dataset would be 282 composables; the whole matrix is still exported and still drawn by the heatmap.
 */
private const val MAX_COLUMNS = 24

/** One finished analysis: a tab per pane, the pane's table, its picture if it has one, and the export button. */
@Composable
fun ResultPane(result: AnalysisResult, part: Int, onEvent: (AppEvent) -> Unit, modifier: Modifier = Modifier) {
    val shown = result.parts.getOrNull(part) ?: return
    Column(modifier = modifier.fillMaxSize()) {
        if (result.parts.size > 1) {
            PrimaryTabRow(selectedTabIndex = part) {
                result.parts.forEachIndexed { index, each ->
                    Tab(
                        selected = index == part,
                        onClick = { onEvent(AppEvent.SelectResult(result.id, index)) },
                        text = { Text(each.title, maxLines = 1, overflow = TextOverflow.Ellipsis) },
                    )
                }
            }
        }
        Row(
            modifier = Modifier.fillMaxWidth().padding(horizontal = 12.dp, vertical = 8.dp),
            horizontalArrangement = Arrangement.spacedBy(12.dp),
        ) {
            Text(shown.title, style = MaterialTheme.typography.titleSmall)
            Text(
                "${shown.table.rowCount} rows x ${shown.table.columnNames.size} columns",
                style = MaterialTheme.typography.bodySmall,
                color = MaterialTheme.colorScheme.onSurfaceVariant,
            )
            OutlinedButton(onClick = { onEvent(AppEvent.ExportResult) }) { Text("Export ${shown.fileName}") }
        }
        val chart = shown.chart
        if (chart != null) {
            ChartPane(chart, modifier = Modifier.weight(CHART_WEIGHT).fillMaxWidth().padding(12.dp))
            HorizontalDivider()
        }
        ResultTablePane(shown.table, modifier = Modifier.weight(1f).fillMaxWidth())
    }
}

@Composable
private fun ChartPane(chart: Chart, modifier: Modifier = Modifier) {
    when (chart) {
        is Chart.Heatmap -> KinshipHeatmap(chart, modifier)
        is Chart.Scatter -> PcaScatter(chart, modifier)
    }
}

/**
 * The table, virtualized down the rows with a [LazyColumn] so a 3,093-row site summary composes only what is on
 * screen (`.cursor/rules/compose-ui.mdc`). Columns are capped rather than lazy, so the row scrolls as one.
 */
@Composable
private fun ResultTablePane(table: TableView, modifier: Modifier = Modifier) {
    val columns = table.columnNames.size.coerceAtMost(MAX_COLUMNS)
    val scroll = rememberScrollState()
    Column(modifier = modifier) {
        Row(modifier = Modifier.horizontalScroll(scroll).padding(horizontal = 12.dp)) {
            for (column in 0 until columns) Cell(table.columnNames[column], header = true)
        }
        HorizontalDivider()
        LazyColumn(modifier = Modifier.weight(1f)) {
            items(count = table.rowCount) { row ->
                Row(modifier = Modifier.horizontalScroll(scroll).padding(horizontal = 12.dp)) {
                    for (column in 0 until columns) Cell(table.cell(row, column), header = false)
                }
            }
        }
        if (columns < table.columnNames.size) {
            Text(
                "Showing the first $columns of ${table.columnNames.size} columns; the export has all of them.",
                style = MaterialTheme.typography.labelSmall,
                color = MaterialTheme.colorScheme.onSurfaceVariant,
                modifier = Modifier.padding(horizontal = 12.dp, vertical = 4.dp),
            )
        }
    }
}

@Composable
private fun Cell(text: String, header: Boolean) {
    Text(
        text,
        style = if (header) MaterialTheme.typography.labelMedium else MaterialTheme.typography.bodySmall,
        maxLines = 1,
        overflow = TextOverflow.Ellipsis,
        modifier = Modifier.width(COLUMN_WIDTH).padding(vertical = 2.dp),
    )
}

/** How much of the pane the picture takes when there is one. */
private const val CHART_WEIGHT = 1.4f
