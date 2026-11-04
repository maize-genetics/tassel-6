package net.maizegenetics.tassel.app.ui

import androidx.compose.foundation.layout.Arrangement
import androidx.compose.foundation.layout.Column
import androidx.compose.foundation.layout.PaddingValues
import androidx.compose.foundation.layout.fillMaxSize
import androidx.compose.foundation.layout.fillMaxWidth
import androidx.compose.foundation.layout.padding
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.items
import androidx.compose.material3.Button
import androidx.compose.material3.HorizontalDivider
import androidx.compose.material3.MaterialTheme
import androidx.compose.material3.NavigationDrawerItem
import androidx.compose.material3.OutlinedButton
import androidx.compose.material3.Text
import androidx.compose.material3.TextButton
import androidx.compose.runtime.Composable
import androidx.compose.ui.Modifier
import androidx.compose.ui.text.style.TextOverflow
import androidx.compose.ui.unit.dp
import net.maizegenetics.tassel.app.state.AnalysisResult
import net.maizegenetics.tassel.app.state.AppEvent
import net.maizegenetics.tassel.app.state.AppState
import net.maizegenetics.tassel.app.state.Dataset
import net.maizegenetics.tassel.app.state.Showing

/** How far a result row is indented under the dataset it belongs to. */
private val RESULT_INDENT = 16.dp

/**
 * The data tree: the datasets of this session with their results under them, the buttons that add to it, and the
 * analyses that can be run on whatever is selected.
 *
 * The analysis list is generated from [AppState.operations], which is whatever the engine says it can run, so an
 * operation added to `TasselOperations.registry` appears here with no change to this file (ADR-0002).
 */
@Composable
fun DataTree(state: AppState, onEvent: (AppEvent) -> Unit, modifier: Modifier = Modifier) {
    Column(modifier = modifier.fillMaxSize().padding(12.dp), verticalArrangement = Arrangement.spacedBy(8.dp)) {
        Text("Data", style = MaterialTheme.typography.titleMedium)
        Button(
            onClick = { onEvent(AppEvent.OpenTutorial) },
            enabled = !state.busy,
            modifier = Modifier.fillMaxWidth(),
        ) {
            Text("Load tutorial dataset")
        }
        OutlinedButton(
            onClick = { onEvent(AppEvent.OpenFile) },
            enabled = !state.busy,
            modifier = Modifier.fillMaxWidth(),
        ) {
            Text("Open file...")
        }
        HorizontalDivider()
        if (state.datasets.isEmpty()) {
            Text(
                "No datasets yet. Load the tutorial dataset, or open a HapMap or VCF file.",
                style = MaterialTheme.typography.bodySmall,
                color = MaterialTheme.colorScheme.onSurfaceVariant,
            )
            return@Column
        }
        LazyColumn(modifier = Modifier.weight(1f), contentPadding = PaddingValues(vertical = 4.dp)) {
            for (dataset in state.datasets) {
                item(key = dataset.id) { DatasetRow(dataset, state.showing, onEvent) }
                items(state.resultsOf(dataset.id), key = { it.id }) { result ->
                    ResultRow(result, state.showing, onEvent)
                }
            }
        }
        AnalysisButtons(state, onEvent)
    }
}

@Composable
private fun DatasetRow(dataset: Dataset, showing: Showing, onEvent: (AppEvent) -> Unit) {
    NavigationDrawerItem(
        selected = showing == Showing.Dataset(dataset.id),
        onClick = { onEvent(AppEvent.Select(dataset.id)) },
        label = {
            Column {
                Text(
                    dataset.name,
                    style = MaterialTheme.typography.bodyMedium,
                    maxLines = 1,
                    overflow = TextOverflow.Ellipsis,
                )
                Text(dataset.dimensions, style = MaterialTheme.typography.bodySmall)
            }
        },
        badge = {
            Column {
                TextButton(onClick = { onEvent(AppEvent.ExportDataset(dataset.id)) }) { Text("Export") }
                TextButton(onClick = { onEvent(AppEvent.Close(dataset.id)) }) { Text("Close") }
            }
        },
    )
}

@Composable
private fun ResultRow(result: AnalysisResult, showing: Showing, onEvent: (AppEvent) -> Unit) {
    NavigationDrawerItem(
        selected = (showing as? Showing.Result)?.id == result.id,
        onClick = { onEvent(AppEvent.SelectResult(result.id)) },
        label = { Text(result.operationId, style = MaterialTheme.typography.bodySmall, maxLines = 1) },
        badge = { TextButton(onClick = { onEvent(AppEvent.CloseResult(result.id)) }) { Text("Close") } },
        modifier = Modifier.padding(start = RESULT_INDENT),
    )
}

/** One button per registered operation, enabled when there is a dataset for it to run on. */
@Composable
private fun AnalysisButtons(state: AppState, onEvent: (AppEvent) -> Unit) {
    HorizontalDivider()
    Text("Analyze", style = MaterialTheme.typography.titleMedium)
    for (operation in state.operations) {
        OutlinedButton(
            onClick = { onEvent(AppEvent.ShowForm(operation.id)) },
            enabled = state.selected != null && !state.busy,
            modifier = Modifier.fillMaxWidth(),
        ) {
            Text(operation.id, maxLines = 1, overflow = TextOverflow.Ellipsis)
        }
    }
}
