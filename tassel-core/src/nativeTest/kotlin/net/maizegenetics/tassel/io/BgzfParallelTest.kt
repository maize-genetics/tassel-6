package net.maizegenetics.tassel.io

import kotlinx.coroutines.runBlocking
import kotlin.test.Test

/**
 * Real threads. [kotlinx.coroutines.test.runTest] would run the workers on one test thread.
 * linuxX64Test is the native gate; macosArm64Test runs the same source set when it is available.
 */
class BgzfParallelTest {
    @Test
    fun parallelDecodeIsIndependentOfWorkerCount() = runBlocking {
        assertParallelBgzf()
    }
}
