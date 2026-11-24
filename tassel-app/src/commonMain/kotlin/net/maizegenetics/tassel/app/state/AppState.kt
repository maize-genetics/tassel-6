package net.maizegenetics.tassel.app.state

import net.maizegenetics.tassel.app.analysis.AnalysisForm
import net.maizegenetics.tassel.app.analysis.ResultPart
import net.maizegenetics.tassel.app.export.GenotypeExportFormat
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.ops.OperationDescriptor
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.phenotype.PhenotypeTable

/**
 * What a [Dataset] holds. Each kind names the operation input type it satisfies, which is how the tree decides what
 * to offer without the app naming an operation (ADR-0002).
 */
sealed interface DatasetData {
    /** The core value an operation takes as input. */
    val value: Any

    /** The `OperationDescriptor.inputType` an operation must have to run on this alone. */
    val inputType: String

    /** The subtitle the tree shows. */
    val dimensions: String

    /** A genotype table. */
    data class Genotypes(val table: GenotypeTable) : DatasetData {
        override val value: Any get() = table
        override val inputType: String get() = "GenotypeTable"
        override val dimensions: String
            get() = "${group(table.numberOfTaxa)} taxa x ${group(table.numberOfSites)} sites"
    }

    /** A phenotype table: observations are rows, attributes (taxa, factors, covariates, data) are columns. */
    data class Phenotypes(val table: PhenotypeTable) : DatasetData {
        override val value: Any get() = table
        override val inputType: String get() = "PhenotypeTable"
        override val dimensions: String
            get() = "${group(table.numberOfObservations)} observations x ${group(table.attributes.size)} attributes"
    }
}

/** One dataset in the data tree. [id] identifies it for the lifetime of the session. */
data class Dataset(val id: Long, val name: String, val data: DatasetData) {
    /** `281 taxa x 3,093 sites` or `281 observations x 4 attributes`, the subtitle the tree shows. */
    val dimensions: String get() = data.dimensions

    /** The operation input type this dataset satisfies on its own. */
    val inputType: String get() = data.inputType

    /** The genotype table, or null for a phenotype. */
    val genotypes: GenotypeTable? get() = (data as? DatasetData.Genotypes)?.table

    /** The phenotype table, or null for a genotype. */
    val phenotype: PhenotypeTable? get() = (data as? DatasetData.Phenotypes)?.table
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

/** One contig the region step can offer. [length] is null when the VCF header has no `##contig` length. */
data class ContigRow(val name: String, val length: Long?)

/**
 * The region step, shown after a bgzipped VCF or BCF and its index are open and before any region is read.
 *
 * [message] is blank-field guidance or a validation error. Null means [text] is a region Load may read.
 */
data class RegionPrompt(val fileName: String, val contigs: List<ContigRow>, val text: String, val message: String?)

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
 * @property checked the datasets ticked for a multi-selection, in the order they were ticked. Empty means the
 *   selection is just [selected].
 * @property form the open analysis form, or null when none is open.
 * @property region the open region step, or null when none is open.
 */
data class AppState(
    val datasets: List<Dataset> = emptyList(),
    val results: List<AnalysisResult> = emptyList(),
    val operations: List<OperationDescriptor> = emptyList(),
    val showing: Showing = Showing.Nothing,
    val checked: List<Long> = emptyList(),
    val form: AnalysisForm? = null,
    val region: RegionPrompt? = null,
    val activity: Activity = Activity.Idle,
) {
    /** The dataset the grid shows, or the one a result was run on, or null. */
    val selected: Dataset?
        get() = when (val showing = showing) {
            is Showing.Dataset -> datasets.firstOrNull { it.id == showing.id }
            is Showing.Result -> selectedResult?.let { result -> datasets.firstOrNull { it.id == result.datasetId } }
            Showing.Nothing -> null
        }

    /**
     * What the Analyze buttons act on: the ticked datasets in tick order, else [selected] alone. A join's inputs are
     * in this order (ADR-0024), which is what the user sees in the tree.
     */
    val selection: List<Dataset>
        get() = if (checked.isEmpty()) {
            listOfNotNull(selected)
        } else {
            checked.mapNotNull { id ->
                datasets.firstOrNull {
                    it.id ==
                        id
                }
            }
        }

    /** The operations that can run on [selection], in the engine's presentation order. */
    val availableOperations: List<OperationDescriptor> get() = operationsFor(operations, selection)

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

    /** Ask the platform for a file and read it, or open the region step when it is a bgzipped VCF with an index. */
    data object OpenFile : AppEvent

    /** Replace the region step's text and validate it. Does not read. */
    data class SetRegion(val text: String) : AppEvent

    /** Put a contig name in the region field, which reads that whole contig once Load is pressed. */
    data class UseContig(val name: String) : AppEvent

    /** Read the region in the open prompt. */
    data object ConfirmRegion : AppEvent

    /** Close the region step without reading. */
    data object CancelRegion : AppEvent

    /** Tick or untick the dataset with this id for a multi-selection. */
    data class ToggleChecked(val id: Long) : AppEvent

    /** Show the dataset with this id in the grid. */
    data class Select(val id: Long) : AppEvent

    /** Drop the dataset with this id, and its results, from the tree. */
    data class Close(val id: Long) : AppEvent

    /**
     * Write the dataset with this id in [format]. A phenotype dataset ignores [format] and is written as
     * `<Phenotype>` text.
     *
     * HapMap is the default, so a caller that does not ask still exports the text TASSEL 5 writes.
     */
    data class ExportDataset(val id: Long, val format: GenotypeExportFormat = GenotypeExportFormat.HapMap) : AppEvent

    /** Open the generated form for this operation. */
    data class ShowForm(val operationId: String) : AppEvent

    /** Put [value] in the form's [name] field. */
    data class SetParameter(val name: String, val value: String) : AppEvent

    /** Close the form without running anything. */
    data object CloseForm : AppEvent

    /** Run the open form against the selected dataset. */
    data object Run : AppEvent

    /**
     * Cancel the load or the running operation by cancelling the coroutine driving it (ADR-0004).
     * A region read stops at its next progress report, and in the browser at its next slice.
     */
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
