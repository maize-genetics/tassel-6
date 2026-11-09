@file:OptIn(ExperimentalAtomicApi::class)

package net.maizegenetics.tassel.io

import kotlinx.coroutines.CompletableDeferred
import kotlinx.coroutines.awaitCancellation
import kotlinx.coroutines.cancelAndJoin
import kotlinx.coroutines.coroutineScope
import kotlinx.coroutines.launch
import kotlinx.coroutines.withTimeout
import kotlinx.coroutines.yield
import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.io.deflate.parseBgzfPayload
import net.maizegenetics.tassel.io.deflate.parseGzipMemberHeader
import net.maizegenetics.tassel.testing.TestData
import kotlin.concurrent.atomics.AtomicInt
import kotlin.concurrent.atomics.ExperimentalAtomicApi
import kotlin.coroutines.cancellation.CancellationException
import kotlin.test.assertContains
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertFalse
import kotlin.test.assertTrue

/**
 * Worker counts 1, 2, 3, and 8 on a real dispatcher. jvmTest and nativeTest call this.
 * wasmJs stays on the one-worker reader: [bgzfWorkerCount] is 1, and this function is not a wasm test.
 */
internal suspend fun assertParallelBgzf() {
    val tutorial = TestData.readBytes(PLAIN_VCF)
    val wheat = sequential(WHEAT_GZ)
    val tutorialLines = lines(TUTORIAL_GZ, workers = 1)
    val wheatLines = lines(WHEAT_GZ, workers = 1)
    for (workers in WORKERS) {
        assertContentEquals(tutorial, scan(TUTORIAL_GZ, workers), "tutorial workers=$workers")
        assertContentEquals(wheat, scan(WHEAT_GZ, workers), "wheat workers=$workers")
        assertSeeks(TUTORIAL_GZ, workers, tutorialLines)
        assertSeeks(WHEAT_GZ, workers, wheatLines)
    }
    assertCrcNamesTheBlockAndStopsTheOtherWorkers()
    assertCancellingTheCallerStopsTheWorkers()
    assertBlocksInFlightStayAtTheBound()
}

private suspend fun assertSeeks(path: String, workers: Int, expected: List<Pair<Long, String>>) {
    openRandomAccessFile("${TestData.root}/$path").use { source ->
        BgzfReader(source, workers = workers).use { reader ->
            for ((at, line) in expected) {
                reader.seek(at)
                assertEquals(at, reader.tell(), "$path workers=$workers at $at")
                assertEquals(line, reader.readLine(), "$path workers=$workers at $at")
            }
        }
    }
}

private suspend fun assertCrcNamesTheBlockAndStopsTheOtherWorkers() {
    val file = TestData.readBytes(BAD_CRC)
    val bad = lastDataBlock(file)
    for (workers in listOf(1, 4)) {
        val waiting = AtomicInt(0)
        ByteArrayRandomAccessSource(file).use { source ->
            val reader = if (workers == 1) {
                BgzfReader(source, workers = 1)
            } else {
                BgzfReader(
                    source,
                    workers = workers,
                    maxBlocksInFlight = bgzfMaxBlocksInFlight(workers),
                    beforeInflate = { address ->
                        if (address != bad) {
                            waiting.addAndFetch(1)
                            awaitCancellation()
                        }
                    },
                )
            }
            reader.use {
                val failure = assertFailsWith<IllegalArgumentException> {
                    withTimeout(TIMEOUT_MS) { it.readAll() }
                }
                assertContains(failure.message ?: "", "CRC")
                assertContains(failure.message ?: "", "$BGZF_BLOCK_FAILURE $bad")
            }
        }
        if (workers > 1) assertTrue(waiting.load() >= 1, "a second worker was waiting when the CRC failed")
    }
}

private suspend fun assertCancellingTheCallerStopsTheWorkers() = coroutineScope {
    val entered = CompletableDeferred<Unit>()
    ByteArrayRandomAccessSource(TestData.readBytes(TUTORIAL_GZ)).use { source ->
        val reader = BgzfReader(
            source,
            workers = 4,
            maxBlocksInFlight = bgzfMaxBlocksInFlight(4),
            beforeInflate = {
                entered.complete(Unit)
                awaitCancellation()
            },
        )
        val job = launch {
            reader.use { it.read(ByteArray(8), 0, 8) }
        }
        try {
            withTimeout(TIMEOUT_MS) {
                entered.await()
                job.cancelAndJoin()
            }
        } catch (e: CancellationException) {
            reader.close()
            throw e
        }
        assertTrue(job.isCancelled)
        assertFalse(job.isActive)
    }
}

private suspend fun assertBlocksInFlightStayAtTheBound() = coroutineScope {
    val workers = 4
    val max = bgzfMaxBlocksInFlight(workers)
    val hold = CompletableDeferred<Unit>()
    ByteArrayRandomAccessSource(TestData.readBytes(TUTORIAL_GZ)).use { source ->
        val reader = BgzfReader(
            source,
            workers = workers,
            maxBlocksInFlight = max,
            beforeInflate = { hold.await() },
        )
        val job = launch { reader.use { it.readAll() } }
        try {
            withTimeout(TIMEOUT_MS) {
                while (reader.blocksInFlight < max) yield()
            }
            assertEquals(max, reader.blocksInFlight)
            repeat(8) { yield() }
            assertEquals(max, reader.blocksInFlight, "the window grew past $max blocks")
            assertEquals(max, reader.peakBlocksInFlight)
        } finally {
            if (!hold.isCompleted) hold.complete(Unit)
            job.join()
        }
    }
}

private suspend fun scan(path: String, workers: Int): ByteArray =
    openRandomAccessFile("${TestData.root}/$path").use { source ->
        BgzfReader(source, workers = workers).use { it.readAll() }
    }

private suspend fun lines(path: String, workers: Int): List<Pair<Long, String>> {
    val found = ArrayList<Pair<Long, String>>()
    openRandomAccessFile("${TestData.root}/$path").use { source ->
        BgzfReader(source, workers = workers).use { reader ->
            while (true) {
                val at = reader.tell()
                val line = reader.readLine() ?: break
                found.add(at to line)
            }
        }
    }
    return found
}

private fun sequential(path: String): ByteArray =
    Buffer().also { it.write(TestData.readBytes(path)) }.decompressing().source.readByteArray()

private fun lastDataBlock(file: ByteArray): Long {
    var start = 0
    var last = -1L
    while (start < file.size) {
        val header = parseGzipMemberHeader(file, start, file.size - start)
        val payload = parseBgzfPayload(file, start, header)
        if (payload.uncompressedSize > 0) last = start.toLong()
        start += header.blockSize + 1
    }
    check(last >= 0) { "no BGZF data block" }
    return last
}

private suspend fun BgzfReader.readAll(): ByteArray {
    val buffer = Buffer()
    val chunk = ByteArray(64 * 1024)
    while (true) {
        val count = read(chunk, 0, chunk.size)
        if (count == 0) break
        buffer.write(chunk, 0, count)
    }
    return buffer.readByteArray()
}

private const val TIMEOUT_MS = 5_000L
private val WORKERS = listOf(1, 2, 3, 8)
private const val PLAIN_VCF = "tutorial/mdp_genotype_w_depth.vcf"
private const val TUTORIAL_GZ = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
private const val WHEAT_GZ = "htslib/wheat/wheat.vcf.gz"
private const val BAD_CRC = "htslib/variants/bad-crc.vcf.gz"
