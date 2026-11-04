package net.maizegenetics.tassel.app.state

import kotlinx.coroutines.CoroutineScope
import kotlinx.coroutines.CoroutineStart
import kotlinx.coroutines.Job
import kotlinx.coroutines.flow.MutableStateFlow
import kotlinx.coroutines.flow.StateFlow
import kotlinx.coroutines.flow.asStateFlow
import kotlinx.coroutines.launch
import kotlinx.coroutines.withContext
import net.maizegenetics.tassel.app.data.GenotypeFile
import net.maizegenetics.tassel.app.data.GenotypeFiles
import net.maizegenetics.tassel.app.data.readGenotypes
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.ops.ProgressReporter
import kotlin.coroutines.CoroutineContext
import kotlin.coroutines.EmptyCoroutineContext
import kotlin.coroutines.cancellation.CancellationException

/**
 * The view model for the one screen the app has: holds [state], takes [AppEvent]s, and owns the coroutine that reads
 * a file. Unidirectional, as `.cursor/rules/compose-ui.mdc` requires.
 *
 * @param files where genotype files come from on this target.
 * @param scope the scope loads run in; it should live as long as the screen.
 * @param readContext the context HapMap parsing runs on. Desktop passes a background dispatcher. The browser has one
 *   thread and no worker here, so it passes [EmptyCoroutineContext] and the parse does occupy the UI thread; that is
 *   a known limitation of this PoC rather than a design choice.
 */
class AppViewModel(
    private val files: GenotypeFiles,
    private val scope: CoroutineScope,
    private val readContext: CoroutineContext = EmptyCoroutineContext,
) {
    private val mutableState = MutableStateFlow(AppState())

    /** The current UI state. */
    val state: StateFlow<AppState> = mutableState.asStateFlow()

    private var load: Job? = null
    private var nextId = 1L

    /** Applies one user action. Returns at once; loads continue in [scope]. */
    fun handle(event: AppEvent) {
        when (event) {
            AppEvent.OpenTutorial -> open { files.tutorial() }
            AppEvent.OpenFile -> open { files.choose() }
            is AppEvent.Select -> mutate { copy(selectedId = event.id) }
            is AppEvent.Close -> mutate {
                val remaining = datasets.filterNot { it.id == event.id }
                copy(
                    datasets = remaining,
                    selectedId = if (selectedId == event.id) remaining.lastOrNull()?.id else selectedId,
                )
            }
            AppEvent.DismissFailure -> mutate { copy(activity = Activity.Idle) }
        }
    }

    /**
     * Reads the file [choose] returns, if any, and adds it to the tree as the selected dataset.
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
                val table = withContext(readContext) { file.readGenotypes(progressReporter()) }
                val dataset = Dataset(nextId++, name, table)
                mutate { copy(datasets = datasets + dataset, selectedId = dataset.id, activity = Activity.Idle) }
            } catch (failure: CancellationException) {
                throw failure
            } catch (failure: Exception) {
                mutate { copy(activity = Activity.Failed(name, failure.message ?: failure.toString())) }
            }
        }
    }

    /** Pushes the reader's reports into [Activity.Loading], ignoring them once the load has moved on. */
    private fun progressReporter() = ProgressReporter { progress ->
        mutate {
            val activity = activity
            if (activity is Activity.Loading) copy(activity = activity.copy(progress = progress)) else this
        }
    }

    /** Replaces the state with the result of [change] applied to it, with the old state as the receiver. */
    private inline fun mutate(change: AppState.() -> AppState) {
        mutableState.value = mutableState.value.change()
    }
}
