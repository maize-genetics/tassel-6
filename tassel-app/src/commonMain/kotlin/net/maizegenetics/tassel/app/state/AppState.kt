package net.maizegenetics.tassel.app.state

import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.ops.Progress

/** One dataset in the data tree. [id] identifies it for the lifetime of the session. */
data class Dataset(val id: Long, val name: String, val table: GenotypeTable) {
    /** `281 taxa x 3,093 sites`, the subtitle the tree shows. */
    val dimensions: String get() = "${group(table.numberOfTaxa)} taxa x ${group(table.numberOfSites)} sites"
}

/** What the app is doing besides showing the tree. */
sealed interface Activity {
    /** Nothing in flight. */
    data object Idle : Activity

    /** [name] is being read; [progress] is the reader's last report. */
    data class Loading(val name: String, val progress: Progress) : Activity

    /** Reading [name] failed with [message]. */
    data class Failed(val name: String, val message: String) : Activity
}

/**
 * Everything the UI draws, immutable. The only mutable thing in the app is the `MutableStateFlow` inside
 * [AppViewModel] that holds one of these (`.cursor/rules/compose-ui.mdc`).
 */
data class AppState(
    val datasets: List<Dataset> = emptyList(),
    val selectedId: Long? = null,
    val activity: Activity = Activity.Idle,
) {
    /** The dataset the grid shows, or null when the tree is empty or nothing is selected. */
    val selected: Dataset? get() = datasets.firstOrNull { it.id == selectedId }
}

/** Everything the user can do. The UI sends these to [AppViewModel.handle] and changes no state itself. */
sealed interface AppEvent {
    /** Read the bundled tutorial dataset. */
    data object OpenTutorial : AppEvent

    /** Ask the platform for a file and read it. */
    data object OpenFile : AppEvent

    /** Show the dataset with this id in the grid. */
    data class Select(val id: Long) : AppEvent

    /** Drop the dataset with this id from the tree. */
    data class Close(val id: Long) : AppEvent

    /** Acknowledge a failure, returning the app to [Activity.Idle]. */
    data object DismissFailure : AppEvent
}

/** `3093` as `3,093`. `tassel-core`'s formatters are for TASSEL 5 output formats, which this is not. */
private fun group(value: Int): String = value.toString().reversed().chunked(3).joinToString(",").reversed()
