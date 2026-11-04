package net.maizegenetics.tassel.cli

import com.github.ajalt.clikt.core.CliktError
import kotlinx.io.IOException
import net.maizegenetics.tassel.ops.InvalidParametersException

/** Exit code for a bad input file or bad options (`.cursor/rules/cli.mdc`). */
const val EXIT_USER_ERROR: Int = 1

/** Exit code for a failure that is not the user's fault. */
const val EXIT_INTERNAL_ERROR: Int = 2

/**
 * Stops with [EXIT_INTERNAL_ERROR]: the CLI has a gap the user cannot work around, such as an
 * operation whose input or output type nothing here can read or write yet.
 */
internal fun unsupported(message: String): Nothing =
    throw CliktError("$message; see tassel-cli", statusCode = EXIT_INTERNAL_ERROR)

/**
 * Runs the body of a command, mapping what the core throws onto the CLI's exit codes: a rejected
 * parameter value or an unreadable file is the user's (exit 1), anything else is ours (exit 2).
 *
 * The core reports every invalid parameter in one [InvalidParametersException] rather than stopping
 * at the first, so they are all printed.
 */
internal fun <T> cliRun(block: () -> T): T = try {
    block()
} catch (e: CliktError) {
    throw e
} catch (e: InvalidParametersException) {
    throw CliktError(
        e.problems.joinToString("\n") { "invalid value for --${kebabCase(it.parameter)}: ${it.message}" },
        statusCode = EXIT_USER_ERROR,
    )
} catch (e: IllegalArgumentException) {
    throw CliktError(e.message ?: "invalid input", statusCode = EXIT_USER_ERROR)
} catch (e: IOException) {
    throw CliktError(e.message ?: "I/O error", statusCode = EXIT_USER_ERROR)
} catch (e: Throwable) {
    e.printStackTrace()
    throw CliktError("internal error: ${e::class.simpleName}: ${e.message}", statusCode = EXIT_INTERNAL_ERROR)
}
