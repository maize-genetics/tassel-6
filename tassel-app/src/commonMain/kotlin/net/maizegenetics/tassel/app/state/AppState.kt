package net.maizegenetics.tassel.app.state

import net.maizegenetics.tassel.app.analysis.AnalysisForm
import net.maizegenetics.tassel.app.analysis.ResultPart
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.ops.OperationDescriptor
import net.maizegenetics.tassel.ops.Progress

/** One dataset in the data tree. [id] identifies it for the lifetime of the session. */
data class Dataset(val id: Long, val name: String, val table: GenotypeTable) {
    /** `281 taxa x 3,093 sites`, the subtitle the tree shows. */
    val dimensions: String get() = "${group(table.numberOfTaxa)} taxa x ${group(table.numberOfSites)} sites"
}

/**
 * One finished analysis, hanging off the dataset it was run on.
 *
 * @property id identifies it for the lifetime of the session.
 * @property datasetId the dataset it was run on.
 * @property operationId the operation that produced it.
 * @property parts the panes it shows, in the order the operation returned them.
 */
data class AnalysisResult(val id: Long, val datasetId: Long, val operationId: String, val parts: List<ResultPart>)

/** What the app is doing besides showing the tree. */
sealed interface Activity {
    /** Nothing in flight. */
    data object Idle : Activity

    /** [name] is being read; [progress] is the reader's last report. */
    data class Loading(val name: String, val progress: Progress) : Activity

    /** [operationId] is running; [progress] is the operation's last report. */
    data class Running(val operationId: String, val progress: Progress) : Activity

    /** [name] failed with [message]; [name] is a file name or an operation id. */
    data class Failed(val name: String, val message: String) : Activity
}

/** What the content pane shows. */
sealed interface Showing {
    /** Nothing loaded, or nothing picked. */
    data object Nothing : Showing

    /** The genotype grid of the dataset with this id. */
    data class Dataset(val id: Long) : Showing

    /** Pane [part] of the result with this id. */
    data class Result(val id: Long, val part: Int) : Showing
}

/**
 * Everything the UI draws, immutable. The only mutable thing in the app is the `MutableStateFlow` inside
 * [AppViewModel] that holds one of these (`.cursor/rules/compose-ui.mdc`).
 *
 * @property operations what the engine can run, which is where forms come from. Empty until the engine answers.
 * @property form the open analysis form, or null when none is open.
 */
data class AppState(
    val datasets: List<Dataset> = emptyList(),
    val results: List<AnalysisResult> = emptyList(),
    val operations: List<OperationDescriptor> = emptyList(),
    val showing: Showing = Showing.Nothing,
    val form: AnalysisForm? = null,
    val activity: Activity = Activity.Idle,
) {
    /** The dataset the grid shows, or the one a result was run on, or null. */
    val selected: Dataset?
        get() = when (val showing = showing) {
            is Showing.Dataset -> datasets.firstOrNull { it.id == showing.id }
            is Showing.Result -> selectedResult?.let { result -> datasets.firstOrNull { it.id == result.datasetId } }
            Showing.Nothing -> null
        }

    /** The result the content pane shows, or null. */
    val selectedResult: AnalysisResult?
        get() = (showing as? Showing.Result)?.let { shown -> results.firstOrNull { it.id == shown.id } }

    /** The pane of [selectedResult] the content pane shows, or null. */
    val selectedPart: ResultPart?
        get() = (showing as? Showing.Result)?.let { shown -> selectedResult?.parts?.getOrNull(shown.part) }

    /** Whether something is in flight, which is what greys out the buttons that would start another. */
    val busy: Boolean get() = activity is Activity.Loading || activity is Activity.Running

    /** The results of [datasetId], in the order they were produced. */
    fun resultsOf(datasetId: Long): List<AnalysisResult> = results.filter { it.datasetId == datasetId }
}

/** Everything the user can do. The UI sends these to [AppViewModel.handle] and changes no state itself. */
sealed interface AppEvent {
    /** Read the bundled tutorial dataset. */
    data object OpenTutorial : AppEvent

    /** Ask the platform for a file and read it. */
    data object OpenFile : AppEvent

    /** Show the dataset with this id in the grid. */
    data class Select(val id: Long) : AppEvent

    /** Drop the dataset with this id, and its results, from the tree. */
    data class Close(val id: Long) : AppEvent

    /** Write the dataset with this id out as HapMap. */
    data class ExportDataset(val id: Long) : AppEvent

    /** Open the generated form for this operation. */
    data class ShowForm(val operationId: String) : AppEvent

    /** Put [value] in the form's [name] field. */
    data class SetParameter(val name: String, val value: String) : AppEvent

    /** Close the form without running anything. */
    data object CloseForm : AppEvent

    /** Run the open form against the selected dataset. */
    data object Run : AppEvent

    /** Cancel the running operation by cancelling the coroutine driving it (ADR-0004). */
    data object Cancel : AppEvent

    /** Show pane [part] of the result with this id. */
    data class SelectResult(val id: Long, val part: Int = 0) : AppEvent

    /** Drop the result with this id. */
    data class CloseResult(val id: Long) : AppEvent

    /** Write the shown pane of the shown result out in the TASSEL 5 text format. */
    data object ExportResult : AppEvent

    /** Acknowledge a failure, returning the app to [Activity.Idle]. */
    data object DismissFailure : AppEvent
}

/** `3093` as `3,093`. `tassel-core`'s formatters are for TASSEL 5 output formats, which this is not. */
private fun group(value: Int): String = value.toString().reversed().chunked(3).joinToString(",").reversed()
