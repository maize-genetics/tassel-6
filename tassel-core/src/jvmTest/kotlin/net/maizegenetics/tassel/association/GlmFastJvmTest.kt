package net.maizegenetics.tassel.association

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.runBlocking
import kotlinx.coroutines.withContext
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * The fast GLM kernel on real threads: `commonTest` runs its workers one after another on a test dispatcher, so this
 * runs them on `Dispatchers.Default` (as the JVM engine does) and checks the bytes do not move.
 */
class GlmFastJvmTest {
    @Test
    fun workersOnTheDefaultDispatcherGiveTheBytesOfOneWorker() = runBlocking {
        val data = AssociationData.of(randomGlmDesign(seed = 21, taxa = 150, markers = 600, missingRate = 0.04))
        val options = GlmOptions(blockSize = 128)
        suspend fun run(workers: Int, droppedRows: DroppedRows) = withContext(Dispatchers.Default) {
            val collector = AssociationCollector(1.0, MaxPRule.Glm, data.traits)
            glm(
                data,
                GenotypeMarkerSource.of(data.data),
                options,
                collector,
                workers = workers,
                tuning = GlmFastTuning(droppedRows = droppedRows),
            )
            val table = collector.finish(AnalysisProvenance("glm"))
            table.records.map { it.toString() }.toList() + table.effects.map { it.toString() }
        }
        assertTrue(Runtime.getRuntime().availableProcessors() >= 1)
        for (dropped in DroppedRows.entries) {
            val one = run(1, dropped)
            for (workers in listOf(2, 7, maxOf(8, associationWorkerCount()))) {
                // Several times each: a scheduling-dependent result would show up as a difference between runs.
                repeat(3) { assertEquals(one, run(workers, dropped), "$dropped with $workers workers") }
            }
        }
    }
}
