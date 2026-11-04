package net.maizegenetics.tassel.app.analysis

import net.maizegenetics.tassel.engine.DataHandle
import net.maizegenetics.tassel.engine.LocalEngine
import net.maizegenetics.tassel.engine.RunRequest
import net.maizegenetics.tassel.ops.OperationDescriptor
import net.maizegenetics.tassel.ops.ProgressReporter

/**
 * How the app reaches an engine.
 *
 * `Engine` itself has no way to hand it data a front end loaded; only [LocalEngine] has `put` and `get`, and
 * `.cursor/rules/compose-ui.mdc` says the UI must not assume results are computed in-process. So the UI talks to
 * this instead: the three calls an analysis screen needs, with [LocalAnalysisEngine] wrapping the only engine the
 * PoC has. A later remote engine uploads in [put] and downloads in [result]; nothing above this changes.
 */
interface AnalysisEngine {
    /** Descriptors of every operation, in the registry's presentation order. Forms are generated from these. */
    suspend fun operations(): List<OperationDescriptor>

    /** Gives [data] to the engine and returns the handle a [RunRequest] names it by. */
    suspend fun put(data: Any): DataHandle

    /**
     * Runs [request], reporting to [progress]. Cancel the calling coroutine to cancel the run; there is no separate
     * cancel call (ADR-0004).
     */
    suspend fun run(request: RunRequest, progress: ProgressReporter): DataHandle

    /** The data behind [handle]. */
    suspend fun result(handle: DataHandle): Any
}

/** [AnalysisEngine] over the in-process [LocalEngine]. */
class LocalAnalysisEngine(private val engine: LocalEngine) : AnalysisEngine {
    override suspend fun operations(): List<OperationDescriptor> = engine.operations()

    override suspend fun put(data: Any): DataHandle = engine.put(data)

    override suspend fun run(request: RunRequest, progress: ProgressReporter): DataHandle =
        engine.run(request, progress)

    override suspend fun result(handle: DataHandle): Any = engine.get(handle)
}
