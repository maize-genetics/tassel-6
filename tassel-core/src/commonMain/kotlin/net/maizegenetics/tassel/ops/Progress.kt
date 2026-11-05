package net.maizegenetics.tassel.ops

import kotlinx.serialization.Serializable

/**
 * Progress of a running operation: [completed] of [total] work units (for example sites), with an
 * optional stage [message]. A [total] of 0 means the amount of work is unknown.
 *
 * @see net.maizegenetics.plugindef.AbstractPlugin.fireProgress
 */
@Serializable
data class Progress(val completed: Long, val total: Long, val message: String? = null) {
    init {
        require(completed >= 0 && total >= 0) { "negative progress: $completed of $total" }
    }

    /** Fraction done in `[0, 1]`; 0 when [total] is unknown. */
    val fraction: Double
        get() = if (total == 0L) 0.0 else (completed.toDouble() / total).coerceIn(0.0, 1.0)
}

/**
 * Receives [Progress] from a running operation. Engines may throw a cancellation exception from
 * [report] when the run has been cancelled, so long loops that report progress also stop promptly.
 *
 * [report] suspends, which makes every progress report a suspension point inside the reader or
 * kernel that emits it (ADR-0019). A kernel only reports; whether reporting also yields the thread
 * is the engine's choice, and [net.maizegenetics.tassel.engine.LocalEngine] yields so the browser's
 * one thread can paint between reports.
 */
fun interface ProgressReporter {
    /** Called by the operation as work completes. */
    suspend fun report(progress: Progress)

    companion object {
        /** Ignores all progress. */
        val None: ProgressReporter = ProgressReporter { }
    }
}

/** Reports [completed] of [total] work units. */
suspend fun ProgressReporter.report(completed: Int, total: Int, message: String? = null) =
    report(Progress(completed.toLong(), total.toLong(), message))

/** Reports [completed] of [total] work units. */
suspend fun ProgressReporter.report(completed: Long, total: Long, message: String? = null) =
    report(Progress(completed, total, message))
