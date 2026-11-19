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
    /** Exclusive compressed offset of the chunk end, or -1 when the window may run to the end of the file. */
    private val chunkEndCoff: Long,
    /** Uncompressed offset of the chunk end. Zero means the chunk ends on a block boundary. */
    private val chunkEndUoff: Int,
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
        var spilled = false
        while (address < source.size) {
            if (chunkEndCoff >= 0) {
                val pastEnd = address > chunkEndCoff || (address == chunkEndCoff && chunkEndUoff == 0)
                if (pastEnd) {
                    if (spilled) break
                    spilled = true
                }
            }
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

/**
 * Ordered deflate of BGZF slices, on [Dispatchers.Default].
 *
 * The caller copies a closed slice and [submit]s it. [workers] coroutines compress them, each with its own
 * compressor. The caller frames the results, in submit order: [write] runs on the caller, which is the only
 * thread that may touch the sink. At most [maxBlocksInFlight] slices are copied and not yet framed.
 * [trySubmit] queues without parking; it returns false when the next slice would have to wait. A compress
 * failure records the cause, cancels the other workers, and fails the [submit] or [drain] that frames that
 * slice.
 *
 * [beforeCompress] runs on the worker, before that slice is compressed. Tests use it to stall.
 */
internal class BgzfEncodeWindow(
    workers: Int,
    private val maxBlocksInFlight: Int,
    private val newCompressor: () -> BgzfCompressor,
    private val beforeCompress: suspend (Int) -> Unit,
    parent: Job? = null,
) {
    private val job = SupervisorJob(parent)
    private val scope = CoroutineScope(job + Dispatchers.Default)
    private val permits = Semaphore(maxBlocksInFlight)
    private val ordered = Channel<Pending>(maxBlocksInFlight)
    private val work = Channel<Task>(maxBlocksInFlight)
    private val failure = AtomicReference<Throwable?>(null)
    private val inFlight = AtomicInt(0)
    private val peak = AtomicInt(0)
    private var nextIndex = 0

    /** Slices copied and not yet framed. */
    val blocksInFlight: Int get() = inFlight.load()

    /** Largest [blocksInFlight] seen. Never above [maxBlocksInFlight]. */
    val peakBlocksInFlight: Int get() = peak.load()

    /** The first compress or frame failure, if one has been recorded. */
    fun takeFailure(): Throwable? = failure.load()

    init {
        require(workers >= 1) { "workers must be positive: $workers" }
        require(maxBlocksInFlight >= 1) { "blocks in flight must be positive: $maxBlocksInFlight" }
        repeat(workers) {
            scope.launch { worker() }
        }
    }

    fun shutdown() {
        job.cancel()
    }

    /**
     * Queues [input] when a permit is free.
     *
     * Returns false when [maxBlocksInFlight] slices are already outstanding. [input] is not copied again;
     * it must stay unchanged until [write] frames it.
     */
    fun trySubmit(input: ByteArray): Boolean {
        takeFailure()?.let { throw it }
        if (!permits.tryAcquire()) return false
        enqueue(input)
        return true
    }

    /**
     * Queues [input], framing the oldest finished slice first when the window is full.
     *
     * [write] runs inline on the caller.
     */
    suspend fun submit(input: ByteArray, write: (ByteArray, ByteArray?) -> Unit) {
        takeFailure()?.let { throw it }
        if (!permits.tryAcquire()) {
            frameOne(write)
            takeFailure()?.let { throw it }
            permits.acquire()
        }
        enqueue(input)
    }

    /** Frames every slice already submitted, in order. [write] runs inline on the caller. */
    suspend fun drain(write: (ByteArray, ByteArray?) -> Unit) {
        while (inFlight.load() > 0) frameOne(write)
    }

    private fun enqueue(input: ByteArray) {
        val index = nextIndex++
        val deferred = CompletableDeferred<ByteArray?>(job)
        noteQueued()
        val published = ordered.trySend(Pending(input, deferred))
        if (!published.isSuccess) {
            noteTaken()
            permits.release()
            deferred.cancel()
            reject(index)
        }
        val accepted = work.trySend(Task(index, input, deferred))
        if (!accepted.isSuccess) {
            val cause = IllegalStateException("BGZF encode workers rejected block $index")
            deferred.completeExceptionally(cause)
            fail(cause)
            throw cause
        }
    }

    private fun reject(index: Int): Nothing {
        val cause = IllegalStateException("BGZF encode queue rejected block $index")
        fail(cause)
        throw cause
    }

    private suspend fun frameOne(write: (ByteArray, ByteArray?) -> Unit) {
        val pending = try {
            ordered.receive()
        } catch (e: Throwable) {
            throw preferred(e)
        }
        try {
            val deflate = try {
                pending.deferred.await()
            } catch (e: Throwable) {
                throw preferred(e)
            }
            try {
                write(pending.input, deflate)
            } catch (e: Throwable) {
                fail(e)
                throw e
            }
        } finally {
            noteTaken()
            permits.release()
        }
    }

    private suspend fun worker() {
        val compressor = newCompressor()
        for (task in work) {
            try {
                beforeCompress(task.index)
                coroutineContext.ensureActive()
                task.deferred.complete(compressor.compress(task.input, 0, task.input.size))
            } catch (e: Throwable) {
                task.deferred.completeExceptionally(e)
                if (e is CancellationException) throw e
                fail(e)
                throw CancellationException("BGZF compress failed", e)
            }
        }
    }

    private fun preferred(fallback: Throwable): Throwable = failure.load() ?: fallback

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

    private class Pending(val input: ByteArray, val deferred: CompletableDeferred<ByteArray?>)

    private class Task(val index: Int, val input: ByteArray, val deferred: CompletableDeferred<ByteArray?>)
}
