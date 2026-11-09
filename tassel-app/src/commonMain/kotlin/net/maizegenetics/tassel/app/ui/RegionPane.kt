package net.maizegenetics.tassel.app.ui

import androidx.compose.foundation.layout.Arrangement
import androidx.compose.foundation.layout.Column
import androidx.compose.foundation.layout.PaddingValues
import androidx.compose.foundation.layout.Row
import androidx.compose.foundation.layout.fillMaxHeight
import androidx.compose.foundation.layout.fillMaxWidth
import androidx.compose.foundation.layout.padding
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.items
import androidx.compose.material3.Button
import androidx.compose.material3.Card
import androidx.compose.material3.MaterialTheme
import androidx.compose.material3.OutlinedTextField
import androidx.compose.material3.Text
import androidx.compose.material3.TextButton
import androidx.compose.runtime.Composable
import androidx.compose.ui.Modifier
import androidx.compose.ui.text.style.TextOverflow
import androidx.compose.ui.unit.dp
import net.maizegenetics.tassel.app.state.AppEvent
import net.maizegenetics.tassel.app.state.ContigRow
import net.maizegenetics.tassel.app.state.RegionPrompt

/**
 * The region step: the index's contigs, a region in the same syntax as `bcftools view -r`, and Load.
 *
 * Clicking a contig fills the field with that contig's name, which is the whole contig, and does not read. A whole
 * chromosome of a wheat panel is an in-memory table; Load is the explicit start. The message under the field is the
 * validation result, shown before any read.
 */
@Composable
fun RegionPane(prompt: RegionPrompt, onEvent: (AppEvent) -> Unit, modifier: Modifier = Modifier) {
    val invalid = prompt.message != null && prompt.text.isNotBlank()
    Card(modifier = modifier.fillMaxHeight().padding(12.dp)) {
        Column(
            modifier = Modifier.padding(12.dp).fillMaxHeight(),
            verticalArrangement = Arrangement.spacedBy(8.dp),
        ) {
            Text(
                prompt.fileName,
                style = MaterialTheme.typography.titleMedium,
                maxLines = 2,
                overflow = TextOverflow.Ellipsis,
            )
            Text(
                "Choose a contig, or enter a region. Load reads only that region.",
                style = MaterialTheme.typography.bodySmall,
                color = MaterialTheme.colorScheme.onSurfaceVariant,
            )
            LazyColumn(
                modifier = Modifier.weight(1f),
                contentPadding = PaddingValues(vertical = 4.dp),
            ) {
                items(prompt.contigs, key = { it.name }) { contig ->
                    TextButton(
                        onClick = { onEvent(AppEvent.UseContig(contig.name)) },
                        modifier = Modifier.fillMaxWidth(),
                    ) {
                        Text(contigLabel(contig), maxLines = 1, overflow = TextOverflow.Ellipsis)
                    }
                }
            }
            OutlinedTextField(
                value = prompt.text,
                onValueChange = { onEvent(AppEvent.SetRegion(it)) },
                modifier = Modifier.fillMaxWidth(),
                label = { Text("Region") },
                placeholder = { Text("chr3B:600000000-650000000") },
                isError = invalid,
                supportingText = {
                    val message = prompt.message ?: return@OutlinedTextField
                    Text(
                        message,
                        color = if (invalid) {
                            MaterialTheme.colorScheme.error
                        } else {
                            MaterialTheme.colorScheme.onSurfaceVariant
                        },
                    )
                },
                singleLine = true,
            )
            Row(horizontalArrangement = Arrangement.spacedBy(8.dp)) {
                Button(
                    onClick = { onEvent(AppEvent.ConfirmRegion) },
                    enabled = prompt.message == null,
                ) {
                    Text("Load")
                }
                TextButton(onClick = { onEvent(AppEvent.CancelRegion) }) { Text("Cancel") }
            }
        }
    }
}

/** `chr3B` with its length, or the name alone when the header has no `##contig` length. */
private fun contigLabel(contig: ContigRow): String = when (val length = contig.length) {
    null -> contig.name
    else -> "${contig.name}  ${group(length)}"
}

/** `851934019` as `851,934,019`. */
private fun group(value: Long): String = value.toString().reversed().chunked(3).joinToString(",").reversed()
