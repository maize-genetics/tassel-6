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
import androidx.compose.foundation.text.KeyboardOptions
import androidx.compose.material3.Button
import androidx.compose.material3.Card
import androidx.compose.material3.DropdownMenu
import androidx.compose.material3.DropdownMenuItem
import androidx.compose.material3.MaterialTheme
import androidx.compose.material3.OutlinedTextField
import androidx.compose.material3.Switch
import androidx.compose.material3.Text
import androidx.compose.material3.TextButton
import androidx.compose.runtime.Composable
import androidx.compose.runtime.getValue
import androidx.compose.runtime.mutableStateOf
import androidx.compose.runtime.remember
import androidx.compose.runtime.setValue
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.text.input.KeyboardType
import androidx.compose.ui.unit.dp
import net.maizegenetics.tassel.app.analysis.AnalysisForm
import net.maizegenetics.tassel.app.analysis.FieldKind
import net.maizegenetics.tassel.app.analysis.FormField
import net.maizegenetics.tassel.app.state.AppEvent

/**
 * The analysis form, generated from one operation's spec.
 *
 * Nothing here names an operation or a parameter. Every row comes from a [FormField], which
 * `net.maizegenetics.tassel.app.analysis.formField` derives from the spec with an exhaustive `when` over the
 * sealed `ParameterDescriptor` hierarchy, so a new operation gets its form for free and a new descriptor type is
 * a compile error rather than a missing field (`.cursor/rules/compose-ui.mdc`, ADR-0002).
 */
@Composable
fun AnalysisFormPane(form: AnalysisForm, running: Boolean, onEvent: (AppEvent) -> Unit, modifier: Modifier = Modifier) {
    Card(modifier = modifier.fillMaxHeight().padding(12.dp)) {
        Column(
            modifier = Modifier.padding(12.dp).fillMaxHeight(),
            verticalArrangement = Arrangement.spacedBy(8.dp),
        ) {
            Text(form.operationId, style = MaterialTheme.typography.titleMedium)
            Text(
                form.title,
                style = MaterialTheme.typography.bodySmall,
                color = MaterialTheme.colorScheme.onSurfaceVariant,
            )
            LazyColumn(
                modifier = Modifier.weight(1f),
                contentPadding = PaddingValues(vertical = 8.dp),
                verticalArrangement = Arrangement.spacedBy(10.dp),
            ) {
                items(form.fields, key = { it.name }) { field ->
                    FormRow(
                        field = field,
                        value = form.values[field.name].orEmpty(),
                        problem = form.problems[field.name],
                        onChange = { onEvent(AppEvent.SetParameter(field.name, it)) },
                    )
                }
            }
            Row(horizontalArrangement = Arrangement.spacedBy(8.dp)) {
                Button(onClick = { onEvent(AppEvent.Run) }, enabled = !running) { Text("Run") }
                TextButton(onClick = { onEvent(AppEvent.CloseForm) }) { Text("Close") }
            }
        }
    }
}

/** One parameter. The `when` is over [FieldKind], which is the only shape the UI knows parameters in. */
@Composable
private fun FormRow(field: FormField, value: String, problem: String?, onChange: (String) -> Unit) {
    Column(verticalArrangement = Arrangement.spacedBy(2.dp)) {
        when (val kind = field.kind) {
            is FieldKind.Toggle -> Row(
                modifier = Modifier.fillMaxWidth(),
                verticalAlignment = Alignment.CenterVertically,
                horizontalArrangement = Arrangement.spacedBy(8.dp),
            ) {
                Switch(
                    checked = value.equals("true", ignoreCase = true),
                    onCheckedChange = { onChange(it.toString()) },
                )
                Text(label(field), style = MaterialTheme.typography.bodyMedium)
            }
            is FieldKind.Choice -> ChoiceRow(field, kind.choices, value, onChange)
            is FieldKind.Number -> TextRow(field, value, onChange, numeric = true)
            FieldKind.Text, FieldKind.TextList, is FieldKind.Path -> TextRow(field, value, onChange, numeric = false)
        }
        Text(
            problem ?: helpLine(field),
            style = MaterialTheme.typography.labelSmall,
            color = if (problem == null) {
                MaterialTheme.colorScheme.onSurfaceVariant
            } else {
                MaterialTheme.colorScheme.error
            },
        )
    }
}

@Composable
private fun TextRow(field: FormField, value: String, onChange: (String) -> Unit, numeric: Boolean) {
    OutlinedTextField(
        value = value,
        onValueChange = onChange,
        label = { Text(label(field)) },
        singleLine = true,
        keyboardOptions = KeyboardOptions(
            keyboardType = if (numeric) KeyboardType.Number else KeyboardType.Text,
        ),
        modifier = Modifier.fillMaxWidth(),
    )
}

@Composable
private fun ChoiceRow(field: FormField, choices: List<String>, value: String, onChange: (String) -> Unit) {
    var open by remember { mutableStateOf(false) }
    Column(modifier = Modifier.fillMaxWidth()) {
        TextButton(onClick = { open = true }, modifier = Modifier.fillMaxWidth()) {
            Text("${label(field)}: ${value.ifEmpty { "(not set)" }}")
        }
        DropdownMenu(expanded = open, onDismissRequest = { open = false }) {
            for (choice in choices) {
                DropdownMenuItem(
                    text = { Text(choice) },
                    onClick = {
                        open = false
                        onChange(choice)
                    },
                )
            }
        }
    }
}

private fun label(field: FormField): String = if (field.required) "${field.name} *" else field.name

/** The spec's help, with anything the control cannot show by itself appended. */
private fun helpLine(field: FormField): String {
    val kind = field.kind
    val extra = when (kind) {
        is FieldKind.Number -> when {
            kind.min != null && kind.max != null -> "from ${kind.min} to ${kind.max}"
            kind.min != null -> "at least ${kind.min}"
            kind.max != null -> "at most ${kind.max}"
            else -> null
        }
        FieldKind.TextList -> "comma separated"
        is FieldKind.Path -> kind.kind.name
        FieldKind.Text, FieldKind.Toggle, is FieldKind.Choice -> null
    }
    return listOfNotNull(field.help.ifEmpty { null }, extra).joinToString(" — ")
}
