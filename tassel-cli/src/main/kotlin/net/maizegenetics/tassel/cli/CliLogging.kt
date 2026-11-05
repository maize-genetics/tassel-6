package net.maizegenetics.tassel.cli

import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.ops.ProgressReporter

/**
 * Everything the CLI says that is not a result goes to standard error, so `tassel summary ... > file`
 * and `tassel kinship ... | less` carry the result alone (`.cursor/rules/cli.mdc`).
 */
object CliLog {
    /** Set by the root command's `--quiet` flag; silences [info] and progress, never [warn] or [error]. */
    var quiet: Boolean = false

    fun info(message: String) {
        if (!quiet) System.err.println("tassel: $message")
    }

    fun warn(message: String) = System.err.println("tassel: warning: $message")

    fun error(message: String) = System.err.println("tassel: error: $message")
}

/**
 * Reports an operation's progress to standard error, once per completed tenth. Operations whose
 * total is unknown (the HapMap reader, until it has read the last site) report nothing.
 */
fun stderrProgress(label: String = "Working"): ProgressReporter {
    var reported = -1
    return ProgressReporter { progress: Progress ->
        if (CliLog.quiet || progress.total <= 0L) return@ProgressReporter
        val tenth = (progress.fraction * 10).toInt()
        if (tenth > reported) {
            reported = tenth
            System.err.println("tassel: ${progress.message ?: label} ${(progress.fraction * 100).toInt()}%")
        }
    }
}
