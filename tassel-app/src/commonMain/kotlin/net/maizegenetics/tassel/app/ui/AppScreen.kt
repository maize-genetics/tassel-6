package net.maizegenetics.tassel.app.ui

import androidx.compose.foundation.layout.Arrangement
import androidx.compose.foundation.layout.Box
import androidx.compose.foundation.layout.Column
import androidx.compose.foundation.layout.Row
import androidx.compose.foundation.layout.fillMaxSize
import androidx.compose.foundation.layout.fillMaxWidth
import androidx.compose.foundation.layout.padding
import androidx.compose.foundation.layout.size
import androidx.compose.foundation.layout.width
import androidx.compose.material3.Card
import androidx.compose.material3.LinearProgressIndicator
import androidx.compose.material3.MaterialTheme
import androidx.compose.material3.Surface
import androidx.compose.material3.Text
import androidx.compose.material3.TextButton
import androidx.compose.material3.VerticalDivider
import androidx.compose.runtime.Composable
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.unit.dp
import net.maizegenetics.tassel.app.grid.GenotypeGrid
import net.maizegenetics.tassel.app.state.Activity
import net.maizegenetics.tassel.app.state.AppEvent
import net.maizegenetics.tassel.app.state.AppState
import net.maizegenetics.tassel.app.state.Showing
import net.maizegenetics.tassel.app.theme.LocalAlleleColors
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.tasselVersion

/** Width of the data tree pane. */
private val TREE_WIDTH = 260.dp

/** Width of the analysis form pane, when one is open. */
private val FORM_WIDTH = 320.dp

/** The one screen the app has: the data tree, the grid or result for what is selected, and the open form. */
@Composable
fun AppScreen(state: AppState, onEvent: (AppEvent) -> Unit) {
    Surface(modifier = Modifier.fillMaxSize()) {
        Row(modifier = Modifier.fillMaxSize()) {
            DataTree(state = state, onEvent = onEvent, modifier = Modifier.width(TREE_WIDTH))
            VerticalDivider()
            Column(modifier = Modifier.weight(1f).fillMaxSize()) {
                when (val activity = state.activity) {
                    is Activity.Loading ->
                        ProgressBanner(activity.progress, "Reading ${activity.name}", AppEvent.Cancel, onEvent)
                    is Activity.Running ->
                        ProgressBanner(activity.progress, "Running ${activity.operationId}", AppEvent.Cancel, onEvent)
                    is Activity.Failed -> FailureBanner(activity, onEvent)
                    Activity.Idle -> Unit
                }
                Content(state, onEvent)
            }
            val region = state.region
            if (region != null) {
                VerticalDivider()
                RegionPane(prompt = region, onEvent = onEvent, modifier = Modifier.width(FORM_WIDTH))
            }
            val form = state.form
            if (form != null) {
                VerticalDivider()
                AnalysisFormPane(
                    form = form,
                    running = state.busy,
                    onEvent = onEvent,
                    modifier = Modifier.width(FORM_WIDTH),
                )
            }
        }
    }
}

@Composable
private fun Content(state: AppState, onEvent: (AppEvent) -> Unit) {
    val showing = state.showing
    val result = state.selectedResult
    val dataset = state.selected
    when {
        showing is Showing.Result && result != null -> ResultPane(result, showing.part, onEvent)
        dataset != null -> {
            GridHeader(dataset.name, dataset.dimensions)
            GenotypeGrid(dataset.table, modifier = Modifier.fillMaxSize())
        }
        else -> EmptyPane()
    }
}

/** Progress of whatever is in flight, with a cancel button when the work can be cancelled. */
@Composable
private fun ProgressBanner(progress: Progress, title: String, cancel: AppEvent?, onEvent: (AppEvent) -> Unit) {
    Column(modifier = Modifier.fillMaxWidth().padding(12.dp)) {
        Row(verticalAlignment = Alignment.CenterVertically, horizontalArrangement = Arrangement.spacedBy(12.dp)) {
            Text(progress.message ?: title, style = MaterialTheme.typography.bodySmall)
            if (cancel != null) TextButton(onClick = { onEvent(cancel) }) { Text("Cancel") }
        }
        // Neither the reader nor every operation knows its total up front, so an unknown total is indeterminate.
        if (progress.total == 0L) {
            LinearProgressIndicator(modifier = Modifier.fillMaxWidth())
        } else {
            LinearProgressIndicator(
                progress = { progress.fraction.toFloat() },
                modifier = Modifier.fillMaxWidth(),
            )
        }
    }
}

@Composable
private fun FailureBanner(activity: Activity.Failed, onEvent: (AppEvent) -> Unit) {
    Card(modifier = Modifier.fillMaxWidth().padding(12.dp)) {
        Row(
            modifier = Modifier.padding(12.dp),
            verticalAlignment = Alignment.CenterVertically,
            horizontalArrangement = Arrangement.spacedBy(12.dp),
        ) {
            Text(
                "${activity.name} failed: ${activity.message}",
                style = MaterialTheme.typography.bodySmall,
                color = MaterialTheme.colorScheme.error,
                modifier = Modifier.fillMaxWidth(DISMISSIBLE_TEXT_FRACTION),
            )
            TextButton(onClick = { onEvent(AppEvent.DismissFailure) }) { Text("Dismiss") }
        }
    }
}

@Composable
private fun EmptyPane() {
    Box(modifier = Modifier.fillMaxSize(), contentAlignment = Alignment.Center) {
        Column(
            horizontalAlignment = Alignment.CenterHorizontally,
            verticalArrangement = Arrangement.spacedBy(8.dp),
        ) {
            Text("TASSEL ${tasselVersion()}", style = MaterialTheme.typography.headlineMedium)
            Text(
                "Load the tutorial dataset, then pick an analysis.",
                style = MaterialTheme.typography.bodyMedium,
                color = MaterialTheme.colorScheme.onSurfaceVariant,
            )
        }
    }
}

/** Name and size of the dataset in the grid, with the allele legend. */
@Composable
private fun GridHeader(name: String, dimensions: String) {
    Row(
        modifier = Modifier.fillMaxWidth().padding(horizontal = 12.dp, vertical = 8.dp),
        verticalAlignment = Alignment.CenterVertically,
        horizontalArrangement = Arrangement.spacedBy(12.dp),
    ) {
        Text(name, style = MaterialTheme.typography.titleSmall)
        Text(
            dimensions,
            style = MaterialTheme.typography.bodySmall,
            color = MaterialTheme.colorScheme.onSurfaceVariant,
        )
        // No fillMaxHeight here: inside a Column the incoming height constraint is everything left over, so a
        // legend that filled it would take the grid's space and leave nothing to draw in.
        Row(
            verticalAlignment = Alignment.CenterVertically,
            horizontalArrangement = Arrangement.spacedBy(6.dp),
        ) {
            for ((label, color) in LocalAlleleColors.current.legend()) {
                Surface(color = color, modifier = Modifier.size(12.dp)) { }
                Text(label, style = MaterialTheme.typography.labelSmall)
            }
        }
    }
}

/** How much of the failure banner the message takes, leaving room for the button. */
private const val DISMISSIBLE_TEXT_FRACTION = 0.8f
