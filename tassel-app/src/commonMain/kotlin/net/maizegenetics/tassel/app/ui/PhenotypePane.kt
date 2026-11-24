package net.maizegenetics.tassel.app.ui

import androidx.compose.foundation.layout.Arrangement
import androidx.compose.foundation.layout.Column
import androidx.compose.foundation.layout.Row
import androidx.compose.foundation.layout.fillMaxSize
import androidx.compose.foundation.layout.fillMaxWidth
import androidx.compose.foundation.layout.padding
import androidx.compose.material3.MaterialTheme
import androidx.compose.material3.Text
import androidx.compose.runtime.Composable
import androidx.compose.runtime.remember
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.unit.dp
import net.maizegenetics.tassel.app.analysis.PhenotypeTableView
import net.maizegenetics.tassel.app.state.Dataset

/**
 * A phenotype dataset as a table, one row per observation. It is the same virtualized table a result uses
 * (`.cursor/rules/compose-ui.mdc`), so a 281-row traits file composes only the rows on screen.
 */
@Composable
fun PhenotypePane(dataset: Dataset, modifier: Modifier = Modifier) {
    val table = dataset.phenotype ?: return
    val view = remember(table) { PhenotypeTableView(table) }
    Column(modifier = modifier.fillMaxSize()) {
        Row(
            modifier = Modifier.fillMaxWidth().padding(horizontal = 12.dp, vertical = 8.dp),
            verticalAlignment = Alignment.CenterVertically,
            horizontalArrangement = Arrangement.spacedBy(12.dp),
        ) {
            Text(dataset.name, style = MaterialTheme.typography.titleSmall)
            Text(
                dataset.dimensions,
                style = MaterialTheme.typography.bodySmall,
                color = MaterialTheme.colorScheme.onSurfaceVariant,
            )
        }
        ResultTablePane(view, modifier = Modifier.weight(1f).fillMaxWidth())
    }
}
