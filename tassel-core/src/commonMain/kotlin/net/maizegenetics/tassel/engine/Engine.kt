package net.maizegenetics.tassel.engine

import kotlinx.serialization.Serializable
import net.maizegenetics.tassel.ops.InvalidParametersException
import net.maizegenetics.tassel.ops.OperationDescriptor
import net.maizegenetics.tassel.ops.ProgressReporter

/**
 * Opaque reference to data held by an [Engine]. Front ends pass handles instead of data, so the same
 * screens work whether data lives in-process ([LocalEngine]) or on a server (a later `RemoteEngine`).
 * A handle is only meaningful to the engine that issued it.
 *
 * @property id engine-assigned identifier.
 */
@Serializable
data class DataHandle(val id: String)

/**
 * A request to run one operation.
 *
 * @property operationId id of the operation to run.
 * @property input handle of the input data, or `null` for operations whose input is `Unit`.
 * @property parameters parameter values as command-line text, keyed by parameter name. Text is the
 *   one form every front end (CLI, forms, MCP, R, Python) can send over any transport.
 */
@Serializable
data class RunRequest(
    val operationId: String,
    val input: DataHandle? = null,
    val parameters: Map<String, String> = emptyMap(),
)

/**
 * Runs operations on data the engine holds (ADR-0004). Every call is asynchronous. To cancel a run,
 * cancel the coroutine that called [run].
 *
 * @see net.maizegenetics.plugindef.AbstractPlugin.performFunction
 */
interface Engine {
    /** Descriptors of the operations this engine can run. */
    suspend fun operations(): List<OperationDescriptor>

    /**
     * Runs [request] and returns a handle to the output, which the engine keeps until [release].
     * [progress] receives progress from the operation; it may be called from another thread.
     *
     * @throws IllegalArgumentException for an unknown operation, an unknown input handle, or an input
     *   of the wrong type.
     * @throws InvalidParametersException if the parameter values are invalid.
     * @throws kotlinx.coroutines.CancellationException if the calling coroutine is cancelled.
     */
    suspend fun run(request: RunRequest, progress: ProgressReporter = ProgressReporter.None): DataHandle

    /** Frees the data behind [handle]. Releasing an unknown or already-released handle does nothing. */
    suspend fun release(handle: DataHandle)
}
