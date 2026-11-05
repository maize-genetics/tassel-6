package net.maizegenetics.tassel.engine

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.ensureActive
import kotlinx.coroutines.job
import kotlinx.coroutines.sync.Mutex
import kotlinx.coroutines.sync.withLock
import kotlinx.coroutines.withContext
import kotlinx.coroutines.yield
import net.maizegenetics.tassel.ops.OperationDescriptor
import net.maizegenetics.tassel.ops.OperationRegistry
import net.maizegenetics.tassel.ops.ProgressReporter
import kotlin.coroutines.CoroutineContext

/**
 * An [Engine] that runs operations in this process and keeps data in memory. This is the only engine
 * in the PoC; it serves the CLI, the desktop app, and the browser (through the Wasm build).
 *
 * Every progress report from a running operation checks for cancellation and yields the dispatcher,
 * so a run stops at its next report and a single-threaded target gets the thread back between
 * reports (ADR-0019). Kernels report; this engine decides what reporting costs.
 *
 * @param registry operations this engine can run.
 * @param context coroutine context operations run in; [Dispatchers.Default] keeps work off the
 *   caller's thread where the platform has more than one.
 */
class LocalEngine(
    private val registry: OperationRegistry,
    private val context: CoroutineContext = Dispatchers.Default,
) : Engine {
    private val mutex = Mutex()
    private val store = HashMap<String, Any>()
    private var nextId = 1L

    /** Stores [value] and returns a new handle to it, for example data loaded by a front end. */
    suspend fun put(value: Any): DataHandle = mutex.withLock {
        DataHandle("local-${nextId++}").also { store[it.id] = value }
    }

    /**
     * The data behind [handle].
     *
     * @throws IllegalArgumentException if this engine holds no data for [handle].
     */
    suspend fun get(handle: DataHandle): Any =
        mutex.withLock { store[handle.id] } ?: throw IllegalArgumentException("unknown data handle '${handle.id}'")

    override suspend fun operations(): List<OperationDescriptor> = registry.descriptors

    override suspend fun run(request: RunRequest, progress: ProgressReporter): DataHandle {
        val operation = registry[request.operationId]
            ?: throw IllegalArgumentException("unknown operation '${request.operationId}'")
        val input = request.input?.let { get(it) } ?: Unit
        val values = operation.bind(request.parameters)
        val output = withContext(context) {
            val job = coroutineContext.job
            operation.runUnchecked(input, values) { report ->
                job.ensureActive()
                // Reporting is the operation's only suspension point (ADR-0019). Handing the
                // dispatcher back here is what lets a single-threaded target draw a frame; on the
                // JVM and Native it is a cheap re-dispatch onto the same pool.
                yield()
                progress.report(report)
            }
        }
        return put(output)
    }

    override suspend fun release(handle: DataHandle) {
        mutex.withLock { store.remove(handle.id) }
    }
}
