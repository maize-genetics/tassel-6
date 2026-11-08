package net.maizegenetics.tassel.io

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.runBlocking
import kotlin.test.Test

/**
 * Overlapping reads on [Dispatchers.Default]. `runTest` would run them on one test thread.
 */
class RandomAccessConcurrencyTest {
    @Test
    fun concurrentReadsOnTheDefaultDispatcher() {
        runBlocking(Dispatchers.Default) {
            RandomAccessSourceTest().checkConcurrentReads()
        }
    }
}
