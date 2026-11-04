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
import net.maizegenetics.tassel.app.theme.LocalAlleleColors
import net.maizegenetics.tassel.tasselVersion

/** Width of the data tree pane. */
private val TREE_WIDTH = 260.dp

/** The one screen the app has: the data tree on the left, the grid for the selected dataset on the right. */
@Composable
fun AppScreen(state: AppState, onEvent: (AppEvent) -> Unit) {
    Surface(modifier = Modifier.fillMaxSize()) {
        Row(modifier = Modifier.fillMaxSize()) {
            DataTree(
                datasets = state.datasets,
                selectedId = state.selectedId,
                busy = state.activity is Activity.Loading,
                onEvent = onEvent,
                modifier = Modifier.width(TREE_WIDTH),
            )
            VerticalDivider()
            Column(modifier = Modifier.fillMaxSize()) {
                when (val activity = state.activity) {
                    is Activity.Loading -> LoadingBanner(activity)
                    is Activity.Failed -> FailureBanner(activity, onEvent)
                    Activity.Idle -> Unit
                }
                val selected = state.selected
                if (selected == null) {
                    EmptyPane()
                } else {
                    GridHeader(selected.name, selected.dimensions)
                    GenotypeGrid(selected.table, modifier = Modifier.fillMaxSize())
                }
            }
        }
    }
}

@Composable
private fun LoadingBanner(activity: Activity.Loading) {
    Column(modifier = Modifier.fillMaxWidth().padding(12.dp)) {
        Text(
            activity.progress.message ?: "Reading ${activity.name}",
            style = MaterialTheme.typography.bodySmall,
        )
        // The reader does not know the site count until the end, so an unknown total gets an indeterminate bar.
        if (activity.progress.total == 0L) {
            LinearProgressIndicator(modifier = Modifier.fillMaxWidth())
        } else {
            LinearProgressIndicator(
                progress = { activity.progress.fraction.toFloat() },
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
                "Could not read ${activity.name}: ${activity.message}",
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
                "Load the tutorial dataset to see the genotype grid.",
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
