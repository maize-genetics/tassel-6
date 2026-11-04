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
import net.maizegenetics.tassel.app.state.AppEvent
import net.maizegenetics.tassel.app.state.Dataset

/**
 * The data tree: the datasets of this session, with the buttons that add to it.
 *
 * It is a [LazyColumn] rather than a `Column` because a session can hold many datasets, and because tassel-6#22 adds
 * result nodes under each of them.
 */
@Composable
fun DataTree(
    datasets: List<Dataset>,
    selectedId: Long?,
    busy: Boolean,
    onEvent: (AppEvent) -> Unit,
    modifier: Modifier = Modifier,
) {
    Column(modifier = modifier.fillMaxSize().padding(12.dp), verticalArrangement = Arrangement.spacedBy(8.dp)) {
        Text("Data", style = MaterialTheme.typography.titleMedium)
        Button(
            onClick = { onEvent(AppEvent.OpenTutorial) },
            enabled = !busy,
            modifier = Modifier.fillMaxWidth(),
        ) {
            Text("Load tutorial dataset")
        }
        OutlinedButton(
            onClick = { onEvent(AppEvent.OpenFile) },
            enabled = !busy,
            modifier = Modifier.fillMaxWidth(),
        ) {
            Text("Open file...")
        }
        HorizontalDivider()
        if (datasets.isEmpty()) {
            Text(
                "No datasets yet. Load the tutorial dataset, or open a HapMap file.",
                style = MaterialTheme.typography.bodySmall,
                color = MaterialTheme.colorScheme.onSurfaceVariant,
            )
            return@Column
        }
        LazyColumn(contentPadding = PaddingValues(vertical = 4.dp)) {
            items(datasets, key = { it.id }) { dataset ->
                NavigationDrawerItem(
                    selected = dataset.id == selectedId,
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
                        TextButton(onClick = { onEvent(AppEvent.Close(dataset.id)) }) { Text("Close") }
                    },
                )
            }
        }
    }
}
