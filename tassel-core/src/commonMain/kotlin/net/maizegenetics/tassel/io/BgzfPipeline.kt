@file:OptIn(ExperimentalAtomicApi::class)

package net.maizegenetics.tassel.io

import kotlinx.coroutines.CompletableDeferred
import kotlinx.coroutines.CoroutineScope
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.Job
import kotlinx.coroutines.SupervisorJob
import kotlinx.coroutines.cancel
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.ensureActive
import kotlinx.coroutines.launch
import kotlinx.coroutines.sync.Semaphore
import net.maizegenetics.tassel.io.deflate.BGZF_MAX_BLOCK
import net.maizegenetics.tassel.io.deflate.BGZF_MAX_UNCOMPRESSED
import net.maizegenetics.tassel.io.deflate.Crc32
import net.maizegenetics.tassel.io.deflate.GZIP_FIXED_HEADER
import net.maizegenetics.tassel.io.deflate.InflateDecoder
import net.maizegenetics.tassel.io.deflate.gzipHeaderBytesRequired
import net.maizegenetics.tassel.io.deflate.inflateBgzf
import net.maizegenetics.tassel.io.deflate.parseBgzfPayload
import net.maizegenetics.tassel.io.deflate.parseGzipMemberHeader
import net.maizegenetics.tassel.io.deflate.truncatedBgzf
import kotlin.concurrent.atomics.AtomicInt
import kotlin.concurrent.atomics.AtomicReference
import kotlin.concurrent.atomics.ExperimentalAtomicApi
import kotlin.coroutines.cancellation.CancellationException
import kotlin.coroutines.coroutineContext

/** Prefix of a CRC, ISIZE, or inflate failure. The message continues with the compressed offset. */
internal const val BGZF_BLOCK_FAILURE: String = "BGZF block at offset"

internal class DecodedBgzfBlock(val address: Long, val data: ByteArray, val length: Int, val nextAddress: Long)

/**
 * Ordered inflate of BGZF blocks from [startAddress], on [Dispatchers.Default].
 *
 * One producer copies compressed members. [workers] coroutines inflate them, each with its own
 * decoder. [take] returns them in file order. At most [maxBlocksInFlight] members are copied and
 * not yet returned. A decode failure records the block offset, cancels the other workers, and
 * fails [take]. Cancelling [parent] cancels this window; this window failing does not cancel [parent].
 *
 * [beforeInflate] runs on the worker, before that block is inflated. Tests use it to stall.
 */
internal class BgzfDecodeWindow(
    private val source: RandomAccessSource,
    private val workers: Int,
    private val maxBlocksInFlight: Int,
    startAddress: Long,
    private val beforeInflate: suspend (Long) -> Unit,
    parent: Job?,
) {
    private val job = SupervisorJob(parent)
    private val scope = CoroutineScope(job + Dispatchers.Default)
    private val permits = Semaphore(maxBlocksInFlight)
    private val ordered = Channel<CompletableDeferred<DecodedBgzfBlock>>(maxBlocksInFlight)
    private val work = Channel<Task>(maxBlocksInFlight)
    private val failure = AtomicReference<Throwable?>(null)
    private val inFlight = AtomicInt(0)
    private val peak = AtomicInt(0)
    private val scratch = ByteArray(BGZF_MAX_BLOCK)

    /** Members copied and not yet returned from [take]. */
    val blocksInFlight: Int get() = inFlight.load()

    /** Largest [blocksInFlight] seen. Never above [maxBlocksInFlight]. */
    val peakBlocksInFlight: Int get() = peak.load()

    init {
        require(workers >= 1) { "workers must be positive: $workers" }
        require(maxBlocksInFlight >= 1) { "blocks in flight must be positive: $maxBlocksInFlight" }
        scope.launch {
            try {
                produce(startAddress)
            } catch (e: CancellationException) {
                throw e
            } catch (e: Throwable) {
                fail(e)
            } finally {
                work.close()
                ordered.close()
            }
        }
        repeat(workers) {
            scope.launch { worker() }
        }
    }

    fun cancel() {
        job.cancel()
    }

    suspend fun take(): DecodedBgzfBlock {
        val deferred = try {
            ordered.receive()
        } catch (e: Throwable) {
            throw preferred(e)
        }
        try {
            return deferred.await()
        } catch (e: CancellationException) {
            throw preferred(e)
        } finally {
            noteTaken()
            permits.release()
        }
    }

    private fun preferred(fallback: Throwable): Throwable {
        val cause = failure.load()
        return if (cause is IllegalArgumentException) cause else fallback
    }

    private fun fail(cause: Throwable) {
        if (cause is CancellationException) return
        failure.compareAndSet(null, cause)
        job.cancel()
    }

    private fun noteQueued() {
        val now = inFlight.addAndFetch(1)
        while (true) {
            val seen = peak.load()
            if (now <= seen || peak.compareAndSet(seen, now)) return
        }
    }

    private fun noteTaken() {
        inFlight.addAndFetch(-1)
    }

    private suspend fun produce(startAddress: Long) {
        var address = startAddress
        while (address < source.size) {
            coroutineContext.ensureActive()
            permits.acquire()
            val member = try {
                readMember(address)
            } catch (e: Throwable) {
                permits.release()
                throw e
            }
            val deferred = CompletableDeferred<DecodedBgzfBlock>()
            var published = false
            try {
                noteQueued()
                ordered.send(deferred)
                published = true
                work.send(Task(address, member, deferred))
            } catch (e: Throwable) {
                if (published) {
                    deferred.completeExceptionally(e)
                } else {
                    noteTaken()
                    permits.release()
                }
                throw e
            }
            address += member.size
        }
    }

    private suspend fun worker() {
        val decoder = InflateDecoder()
        val crc = Crc32()
        val output = ByteArray(BGZF_MAX_UNCOMPRESSED)
        for (task in work) {
            try {
                task.deferred.complete(decode(task.address, task.member, decoder, crc, output))
            } catch (e: Throwable) {
                task.deferred.completeExceptionally(e)
                if (e is CancellationException) throw e
                fail(e)
                throw CancellationException("BGZF decode failed", e)
            }
        }
    }

    private suspend fun decode(
        address: Long,
        member: ByteArray,
        decoder: InflateDecoder,
        crc: Crc32,
        output: ByteArray,
    ): DecodedBgzfBlock {
        beforeInflate(address)
        coroutineContext.ensureActive()
        val header = parseGzipMemberHeader(member, 0, member.size)
        if (header.blockSize < 0) truncatedBgzf("BGZF block size")
        val payload = parseBgzfPayload(member, 0, header)
        val written = try {
            inflateBgzf(decoder, crc, member, payload, output)
        } catch (e: IllegalArgumentException) {
            throw IllegalArgumentException("$BGZF_BLOCK_FAILURE $address: ${e.message}", e)
        }
        val data = if (written == 0) EMPTY else output.copyOf(written)
        return DecodedBgzfBlock(address, data, written, address + member.size)
    }

    private suspend fun readMember(address: Long): ByteArray {
        var needed = GZIP_FIXED_HEADER
        while (true) {
            readInto(address, needed)
            val required = gzipHeaderBytesRequired(scratch, 0, needed)
            if (required in 0..BGZF_MAX_BLOCK) {
                val header = parseGzipMemberHeader(scratch, 0, needed)
                if (header.blockSize < 0) truncatedBgzf("BGZF block size")
                val total = header.blockSize + 1
                if (total > BGZF_MAX_BLOCK || total < needed) truncatedBgzf("BGZF block size")
                if (total > needed) readRange(address + needed, needed, total - needed)
                return scratch.copyOf(total)
            }
            val want = if (required < 0) -required else required
            if (want > BGZF_MAX_BLOCK) truncatedBgzf("gzip extra field")
            needed = want
        }
    }

    private suspend fun readInto(address: Long, length: Int) {
        if (address < 0 || length < 0 || address > source.size || source.size - address < length) {
            truncatedBgzf("BGZF block")
        }
        val got = source.readAt(address, scratch, 0, length)
        if (got < length) truncatedBgzf("BGZF block")
    }

    private suspend fun readRange(address: Long, destOffset: Int, length: Int) {
        if (address < 0 || length < 0 || address > source.size || source.size - address < length) {
            truncatedBgzf("BGZF block")
        }
        val got = source.readAt(address, scratch, destOffset, length)
        if (got < length) truncatedBgzf("BGZF block")
    }

    private class Task(val address: Long, val member: ByteArray, val deferred: CompletableDeferred<DecodedBgzfBlock>)

    private companion object {
        val EMPTY = ByteArray(0)
    }
}

/** True when [this] begins with a BGZF member. Plain gzip and a short header are false. */
suspend fun RandomAccessSource.isBgzf(): Boolean {
    if (size < GZIP_FIXED_HEADER) return false
    val scratch = ByteArray(minOf(size, BGZF_MAX_BLOCK.toLong()).toInt())
    var needed = GZIP_FIXED_HEADER
    while (true) {
        if (needed > scratch.size) return false
        val got = readAt(0, scratch, 0, needed)
        if (got < needed) return false
        val required = gzipHeaderBytesRequired(scratch, 0, got)
        if (required in 0..scratch.size) {
            val header = try {
                parseGzipMemberHeader(scratch, 0, got)
            } catch (_: IllegalArgumentException) {
                return false
            }
            return header.blockSize >= 0
        }
        val want = if (required < 0) -required else required
        if (want > BGZF_MAX_BLOCK || want > scratch.size) return false
        needed = want
    }
}
