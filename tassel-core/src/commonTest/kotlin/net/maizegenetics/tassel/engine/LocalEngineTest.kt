package net.maizegenetics.tassel.engine

import kotlinx.coroutines.CompletableDeferred
import kotlinx.coroutines.Job
import kotlinx.coroutines.cancelAndJoin
import kotlinx.coroutines.launch
import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.ops.InvalidParametersException
import net.maizegenetics.tassel.ops.OperationRegistry
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.countBlocking
import net.maizegenetics.tassel.ops.readLines
import net.maizegenetics.tassel.ops.scaleValues
import net.maizegenetics.tassel.ops.spin
import kotlin.coroutines.EmptyCoroutineContext
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNotEquals
import kotlin.test.assertTrue

class LocalEngineTest {
    private val registry = OperationRegistry(scaleValues, readLines, spin, countBlocking)

    @Test
    fun runsOperationOnHandleWithProgress() = runTest {
        val engine = LocalEngine(registry)
        val input = engine.put(listOf(1.0, 2.0))
        val reports = mutableListOf<Progress>()
        val output = engine.run(
            RunRequest("scale-values", input, mapOf("factor" to "3", "repeats" to "2")),
        ) { reports += it }

        assertNotEquals(input, output)
        assertEquals(listOf(9.0, 18.0), engine.get(output))
        assertEquals(listOf(Progress(1, 2), Progress(2, 2)), reports)
        assertEquals(registry.descriptors, engine.operations())
    }

    @Test
    fun runsOperationWithoutInput() = runTest {
        val engine = LocalEngine(registry)
        val output = engine.run(RunRequest("read-lines", parameters = mapOf("path" to "a.txt")))
        assertEquals("read a.txt", engine.get(output))
        engine.release(output)
        engine.release(output)
        assertFailsWith<IllegalArgumentException> { engine.get(output) }
    }

    @Test
    fun cancellingTheCallerStopsSuspendingOperation() = runTest {
        val engine = LocalEngine(registry, EmptyCoroutineContext)
        val started = CompletableDeferred<Unit>()
        var reports = 0
        val job = launch {
            engine.run(RunRequest("spin")) {
                reports++
                started.complete(Unit)
            }
        }
        started.await()
        job.cancelAndJoin()
        assertTrue(job.isCancelled)
        val afterCancel = reports
        testScheduler.advanceUntilIdle()
        assertEquals(afterCancel, reports)
    }

    @Test
    fun progressReportsAreCancellationPointsForBlockingKernels() = runTest {
        val engine = LocalEngine(registry, EmptyCoroutineContext)
        val seen = mutableListOf<Long>()
        lateinit var job: Job
        job = launch {
            engine.run(RunRequest("count-blocking", parameters = mapOf("total" to "100"))) {
                seen += it.completed
                if (it.completed == 3L) job.cancel()
            }
        }
        job.join()
        assertTrue(job.isCancelled)
        assertEquals(listOf(1L, 2L, 3L), seen)
    }

    @Test
    fun rejectsBadRequests() = runTest {
        val engine = LocalEngine(registry)
        val text = engine.put("not a list")
        assertFailsWith<IllegalArgumentException> { engine.run(RunRequest("missing")) }
        assertFailsWith<IllegalArgumentException> { engine.run(RunRequest("scale-values", DataHandle("nope"))) }
        assertFailsWith<IllegalArgumentException> { engine.run(RunRequest("scale-values", text)) }
        assertFailsWith<IllegalArgumentException> { engine.run(RunRequest("scale-values")) }
        val e = assertFailsWith<InvalidParametersException> {
            engine.run(RunRequest("scale-values", text, mapOf("factor" to "11")), ProgressReporter.None)
        }
        assertEquals("factor", e.problems.single().parameter)
    }
}
