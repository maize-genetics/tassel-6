package net.maizegenetics.tassel.io

import kotlinx.coroutines.runBlocking
import kotlin.test.Test

/**
 * Real threads. [kotlinx.coroutines.test.runTest] would run the workers on one test thread.
 */
class BgzfParallelTest {
    @Test
    fun parallelDecodeIsIndependentOfWorkerCount() = runBlocking {
        assertParallelBgzf()
    }
}
