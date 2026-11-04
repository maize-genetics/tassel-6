package net.maizegenetics.tassel.app.testing

import kotlinx.io.Buffer
import kotlinx.io.Sink
import kotlinx.io.readString
import net.maizegenetics.tassel.app.analysis.AnalysisEngine
import net.maizegenetics.tassel.app.analysis.LocalAnalysisEngine
import net.maizegenetics.tassel.app.export.ResultFiles
import net.maizegenetics.tassel.engine.LocalEngine
import net.maizegenetics.tassel.ops.OperationRegistry
import net.maizegenetics.tassel.ops.TasselOperations
import kotlin.coroutines.EmptyCoroutineContext

/** Fixtures for the analysis half of the app. Nothing here reads a file or needs a UI toolkit. */

/**
 * An [AnalysisEngine] over a real [LocalEngine], so the tests exercise the engine contract the app relies on
 * (text parameters, progress reports that double as cancellation checks) rather than a stand-in for it.
 *
 * [EmptyCoroutineContext] keeps the operation on the test dispatcher, which is what makes a run deterministic.
 */
fun testEngine(registry: OperationRegistry = TasselOperations.registry): AnalysisEngine =
    LocalAnalysisEngine(LocalEngine(registry, EmptyCoroutineContext))

/** A [ResultFiles] that keeps what was written instead of saving it, and can pretend the user cancelled. */
class RecordingResultFiles(private val accept: Boolean = true) : ResultFiles {
    /** File names the app asked to save, in order. */
    val names: MutableList<String> = mutableListOf()

    /** What was written under each of [names]. */
    val contents: MutableList<String> = mutableListOf()

    override suspend fun save(fileName: String, write: (Sink) -> Unit): Boolean {
        names += fileName
        val buffer = Buffer()
        write(buffer)
        contents += buffer.readString()
        return accept
    }
}

/** The text a writer puts in a sink, which is how a test compares an export against the core's own writer. */
fun written(write: (Sink) -> Unit): String = Buffer().also(write).readString()
