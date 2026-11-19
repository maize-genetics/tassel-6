package net.maizegenetics.tassel.io

import kotlinx.coroutines.CompletableDeferred
import kotlinx.coroutines.coroutineScope
import kotlinx.coroutines.launch
import kotlinx.coroutines.withTimeout
import kotlinx.coroutines.yield
import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.io.deflate.BGZF_BLOCK_SIZE
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.sha256
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Worker counts 1, 2, 3, and 8 on a real dispatcher. jvmTest and nativeTest call this.
 * wasmJs stays on the one-worker writer: [bgzfWorkerCount] is 1, and this function is not a wasm test.
 * The tutorial digest is the #131 SHA-256.
 */
internal suspend fun assertParallelBgzfWrite() {
    val tutorial = TestData.readBytes(TUTORIAL_VCF)
    val one = write(tutorial, BgzfLayout.Vcf, workers = 1)
    assertEquals(TUTORIAL_COMPRESSED_BYTES, one.size)
    assertEquals(COMPRESSED_TUTORIAL_SHA256, sha256(one), sha256(one))
    for (workers in WORKERS) {
        assertContentEquals(one, write(tutorial, BgzfLayout.Vcf, workers), "tutorial workers=$workers")
    }
    val packed = ByteArray(BGZF_BLOCK_SIZE * WORKERS.last() + 100) { (it and 0xFF).toByte() }
    val packedOne = write(packed, BgzfLayout.Bcf, workers = 1)
    for (workers in WORKERS) {
        assertContentEquals(packedOne, write(packed, BgzfLayout.Bcf, workers), "bcf workers=$workers")
    }
    assertBlocksInFlightStayAtTheBound()
    assertCompressFailureDoesNotHang()
}

private fun write(plain: ByteArray, layout: BgzfLayout, workers: Int): ByteArray {
    val sink = Buffer()
    BgzfWriter(sink, layout, workers).use { it.write(plain) }
    return sink.readByteArray()
}

private suspend fun assertBlocksInFlightStayAtTheBound() = coroutineScope {
    val workers = 4
    val max = bgzfMaxBlocksInFlight(workers)
    val hold = CompletableDeferred<Unit>()
    val window = BgzfEncodeWindow(
        workers = workers,
        maxBlocksInFlight = max,
        newCompressor = { BgzfCompressor.deflate() },
        beforeCompress = { _ -> hold.await() },
    )
    val job = launch {
        repeat(max + 4) {
            window.submit(ByteArray(64)) { _, _ -> }
        }
        window.drain { _, _ -> }
    }
    try {
        withTimeout(TIMEOUT_MS) {
            while (window.blocksInFlight < max) yield()
        }
        assertEquals(max, window.blocksInFlight)
        repeat(8) { yield() }
        assertEquals(max, window.blocksInFlight, "the window grew past $max blocks")
        assertEquals(max, window.peakBlocksInFlight)
    } finally {
        if (!hold.isCompleted) hold.complete(Unit)
        job.join()
        window.shutdown()
    }
}

private suspend fun assertCompressFailureDoesNotHang() {
    val window = BgzfEncodeWindow(
        workers = 2,
        maxBlocksInFlight = bgzfMaxBlocksInFlight(2),
        newCompressor = { BgzfCompressor { _, _, _ -> error("deflate failed") } },
        beforeCompress = { _ -> },
    )
    val thrown = try {
        withTimeout(TIMEOUT_MS) {
            repeat(4) { window.submit(ByteArray(32)) { _, _ -> } }
            window.drain { _, _ -> }
        }
        null
    } catch (e: Throwable) {
        e
    } finally {
        window.shutdown()
    }
    assertTrue(thrown is IllegalStateException, "expected deflate failure, got $thrown")
}

private const val TIMEOUT_MS = 5_000L
private val WORKERS = listOf(1, 2, 3, 8)
private const val TUTORIAL_VCF = "tutorial/mdp_genotype_w_depth.vcf"
private const val TUTORIAL_COMPRESSED_BYTES = 983_502
private const val COMPRESSED_TUTORIAL_SHA256 =
    "ce705a1f19a5b7e67d06edf06ad2d80c993a589fc0724c3c9a4957f52390550f"
