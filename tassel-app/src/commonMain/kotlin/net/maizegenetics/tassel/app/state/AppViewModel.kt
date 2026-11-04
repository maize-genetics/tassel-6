package net.maizegenetics.tassel.app.state

import kotlinx.coroutines.CoroutineScope
import kotlinx.coroutines.CoroutineStart
import kotlinx.coroutines.Job
import kotlinx.coroutines.flow.MutableStateFlow
import kotlinx.coroutines.flow.StateFlow
import kotlinx.coroutines.flow.asStateFlow
import kotlinx.coroutines.launch
import kotlinx.coroutines.withContext
import kotlinx.coroutines.yield
import net.maizegenetics.tassel.app.analysis.AnalysisEngine
import net.maizegenetics.tassel.app.analysis.analysisForm
import net.maizegenetics.tassel.app.analysis.genotypeFileName
import net.maizegenetics.tassel.app.analysis.resultParts
import net.maizegenetics.tassel.app.analysis.withProblems
import net.maizegenetics.tassel.app.data.GenotypeFile
import net.maizegenetics.tassel.app.data.GenotypeFiles
import net.maizegenetics.tassel.app.data.readGenotypes
import net.maizegenetics.tassel.app.export.ResultFiles
import net.maizegenetics.tassel.engine.DataHandle
import net.maizegenetics.tassel.engine.RunRequest
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.ops.InvalidParametersException
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.ops.ProgressReporter
import kotlin.coroutines.CoroutineContext
import kotlin.coroutines.EmptyCoroutineContext
import kotlin.coroutines.cancellation.CancellationException

/**
 * The view model for the one screen the app has: holds [state], takes [AppEvent]s, and owns the coroutines that
 * read a file and run an analysis. Unidirectional, as `.cursor/rules/compose-ui.mdc` requires.
 *
 * @param files where genotype files come from on this target.
 * @param engine the engine analyses run on. The UI never calls an analysis directly.
 * @param results where an exported result is written on this target.
 * @param scope the scope loads and runs happen in; it should live as long as the screen.
 * @param readContext the context HapMap parsing runs on. Desktop passes a background dispatcher. The browser has
 *   one thread and no worker here, so it passes [EmptyCoroutineContext] and the parse does occupy the UI thread;
 *   that is a known limitation of this PoC rather than a design choice (tassel-6#54).
 */
class AppViewModel(
    private val files: GenotypeFiles,
    private val engine: AnalysisEngine,
    private val results: ResultFiles,
    private val scope: CoroutineScope,
    private val readContext: CoroutineContext = EmptyCoroutineContext,
) {
    private val mutableState = MutableStateFlow(AppState())

    /** The current UI state. */
    val state: StateFlow<AppState> = mutableState.asStateFlow()

    private var load: Job? = null
    private var run: Job? = null
    private var export: Job? = null
    private var nextId = 1L

    /** Handles the engine issued for each dataset, so a second analysis of the same data uploads nothing again. */
    private val handles = mutableMapOf<Long, DataHandle>()

    init {
        // Forms are generated from whatever the engine says it can run, so the app never names an operation.
        scope.launch { mutate { copy(operations = engine.operations()) } }
    }

    /** Applies one user action. Returns at once; work continues in [scope]. */
    fun handle(event: AppEvent) {
        when (event) {
            AppEvent.OpenTutorial -> open { files.tutorial() }
            AppEvent.OpenFile -> open { files.choose() }
            is AppEvent.Select -> mutate { copy(showing = Showing.Dataset(event.id)) }
            is AppEvent.Close -> close(event.id)
            is AppEvent.ExportDataset -> exportDataset(event.id)
            is AppEvent.ShowForm -> showForm(event.operationId)
            is AppEvent.SetParameter -> mutate { copy(form = form?.with(event.name, event.value)) }
            AppEvent.CloseForm -> mutate { copy(form = null) }
            AppEvent.Run -> runAnalysis()
            AppEvent.Cancel -> run?.cancel()
            is AppEvent.SelectResult -> mutate { copy(showing = Showing.Result(event.id, event.part)) }
            is AppEvent.CloseResult -> mutate {
                copy(results = results.filterNot { it.id == event.id }, showing = showingAfterClose(event.id))
            }
            AppEvent.ExportResult -> exportResult()
            AppEvent.DismissFailure -> mutate { copy(activity = Activity.Idle) }
        }
    }

    /**
     * Reads the file [choose] returns, if any, and adds it to the tree as the shown dataset.
     *
     * One load at a time: a second request cancels the first, so a user who picks another file while a big one is
     * parsing does not leave two readers competing for the browser's single thread.
     */
    private fun open(choose: suspend () -> GenotypeFile?) {
        load?.cancel()
        // UNDISPATCHED so that clicking "Open file" shows the file dialog in the same frame as the click.
        load = scope.launch(start = CoroutineStart.UNDISPATCHED) {
            var name = ""
            try {
                val file = choose() ?: return@launch
                name = file.name
                mutate { copy(activity = Activity.Loading(name, Progress(0, 0, "Opening $name"))) }
                val table = withContext(readContext) { file.readGenotypes(loadProgress()) }
                add(name, table)
                mutate { copy(activity = Activity.Idle) }
            } catch (failure: CancellationException) {
                throw failure
            } catch (failure: Exception) {
                mutate { copy(activity = Activity.Failed(name, failure.message ?: failure.toString())) }
            }
        }
    }

    /** Adds [table] to the tree under [name] and shows it. */
    private fun add(name: String, table: GenotypeTable): Dataset {
        val dataset = Dataset(nextId++, name, table)
        mutate { copy(datasets = datasets + dataset, showing = Showing.Dataset(dataset.id)) }
        return dataset
    }

    /** Drops a dataset and the results that hang off it, falling back to the last dataset left. */
    private fun close(id: Long) {
        handles.remove(id)
        mutate {
            val remaining = datasets.filterNot { it.id == id }
            val keep = when (val shown = showing) {
                is Showing.Dataset -> shown.takeIf { shown.id != id }
                is Showing.Result -> shown.takeIf { results.firstOrNull { it.id == shown.id }?.datasetId != id }
                Showing.Nothing -> null
            }
            copy(
                datasets = remaining,
                results = results.filterNot { it.datasetId == id },
                showing = keep ?: remaining.lastOrNull()?.let { Showing.Dataset(it.id) } ?: Showing.Nothing,
            )
        }
    }

    /** Opens the form for [operationId], prefilled with the spec's defaults. */
    private fun showForm(operationId: String) = mutate {
        val descriptor = operations.firstOrNull { it.id == operationId } ?: return@mutate this
        copy(form = analysisForm(descriptor))
    }

    /** What to show once the result with this [id] is gone: the dataset it belonged to, or nothing. */
    private fun AppState.showingAfterClose(id: Long): Showing {
        val shown = showing
        if (shown !is Showing.Result || shown.id != id) return shown
        return results.firstOrNull { it.id == id }?.let { Showing.Dataset(it.datasetId) } ?: Showing.Nothing
    }

    /**
     * Runs the open form against the shown dataset.
     *
     * Cancellation is the coroutine's: [AppEvent.Cancel] cancels this job, `LocalEngine` turns every progress
     * report into a cancellation check, and the operation stops at its next report. There is no separate cancel
     * call to make (ADR-0004).
     *
     * The [yield] before the run exists for the browser, where the UI and the analysis share one thread: it lets
     * the "Running" banner reach the screen before the operation takes the thread. It does not make the operation
     * itself interruptible there; see the PR for tassel-6#22 and tassel-6#54.
     */
    private fun runAnalysis() {
        val form = state.value.form ?: return
        val dataset = state.value.selected ?: return
        run?.cancel()
        run = scope.launch {
            try {
                mutate {
                    copy(activity = Activity.Running(form.operationId, Progress(0, 0, "Running ${form.operationId}")))
                }
                yield()
                val handle = handles.getOrPut(dataset.id) { engine.put(dataset.table) }
                val request = RunRequest(form.operationId, handle, form.parameters)
                val output = engine.result(engine.run(request, runProgress()))
                yield()
                show(dataset, form.operationId, output)
                mutate { copy(activity = Activity.Idle, form = null) }
            } catch (failure: CancellationException) {
                // Only if no later run has started, so a cancelled run cannot clear its successor's banner.
                if (run === coroutineContext[Job]) mutate { copy(activity = Activity.Idle) }
                throw failure
            } catch (invalid: InvalidParametersException) {
                // Every bad field at once, back on the form the user is looking at.
                mutate { copy(activity = Activity.Idle, form = form.withProblems(invalid)) }
            } catch (failure: Exception) {
                mutate { copy(activity = Activity.Failed(form.operationId, failure.message ?: failure.toString())) }
            }
        }
    }

    /**
     * Puts [output] where the user can see it: a genotype table becomes another dataset in the tree, because that
     * is what a filter produces and what the next analysis takes as input; everything else becomes a result node.
     */
    private fun show(dataset: Dataset, operationId: String, output: Any) {
        if (output is GenotypeTable) {
            add("${dataset.name} + $operationId", output)
            return
        }
        val result = AnalysisResult(
            id = nextId++,
            datasetId = dataset.id,
            operationId = operationId,
            parts = resultParts(output, baseName(dataset.name, operationId)),
        )
        mutate { copy(results = results + result, showing = Showing.Result(result.id, 0)) }
    }

    /** Writes the shown pane of the shown result, in the TASSEL 5 text format its writer produces. */
    private fun exportResult() {
        val part = state.value.selectedPart ?: return
        export?.cancel()
        export = scope.launch(start = CoroutineStart.UNDISPATCHED) {
            saving(part.fileName) { results.save(part.fileName, part.write) }
        }
    }

    /** Writes a dataset out as HapMap, which is how TASSEL 5 exports genotypes. */
    private fun exportDataset(id: Long) {
        val dataset = state.value.datasets.firstOrNull { it.id == id } ?: return
        val fileName = genotypeFileName(dataset.name)
        export?.cancel()
        export = scope.launch(start = CoroutineStart.UNDISPATCHED) {
            saving(fileName) { results.save(fileName) { sink -> HapMap.write(dataset.table, sink) } }
        }
    }

    /** Reports a failed save without turning a cancelled dialog into one. */
    private suspend fun saving(fileName: String, save: suspend () -> Boolean) {
        try {
            save()
        } catch (failure: CancellationException) {
            throw failure
        } catch (failure: Exception) {
            mutate { copy(activity = Activity.Failed(fileName, failure.message ?: failure.toString())) }
        }
    }

    /** Pushes the reader's reports into [Activity.Loading], ignoring them once the load has moved on. */
    private fun loadProgress() = ProgressReporter { progress ->
        mutate { (activity as? Activity.Loading)?.let { copy(activity = it.copy(progress = progress)) } ?: this }
    }

    /** Pushes the operation's reports into [Activity.Running], ignoring them once the run has moved on. */
    private fun runProgress() = ProgressReporter { progress ->
        mutate { (activity as? Activity.Running)?.let { copy(activity = it.copy(progress = progress)) } ?: this }
    }

    /** Replaces the state with the result of [change] applied to it, with the old state as the receiver. */
    private inline fun mutate(change: AppState.() -> AppState) {
        mutableState.value = mutableState.value.change()
    }
}

/** `mdp_genotype.kinship`, the stem `tassel-cli`'s `--output` would be given for the same run. */
internal fun baseName(datasetName: String, operationId: String): String =
    datasetName.substringBefore('.').ifEmpty { "result" } + "." + operationId
