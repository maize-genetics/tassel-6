@file:OptIn(ExperimentalAtomicApi::class)

package net.maizegenetics.tassel.io

import kotlin.concurrent.atomics.AtomicBoolean
import kotlin.concurrent.atomics.ExperimentalAtomicApi

/**
 * Bytes that can be read at any position.
 *
 * Reads do not share a file offset, and they take no lock, so overlapping reads from several coroutines are safe.
 * That is what parallel BGZF decoding needs. The read suspends because a browser `Blob` read is asynchronous
 * (ADR-0021). The caller owns [readAt]'s destination array.
 *
 * A short [readAt] result means the source ended, including zero when [position] is at or past [size]. It does not
 * mean the platform returned a short syscall. [close] is idempotent. [size] and [readAt] fail after close.
 *
 * Readers take this type and never a path.
 */
interface RandomAccessSource : AutoCloseable {
    /** Size in bytes, fixed when the source is opened. */
    val size: Long

    /**
     * Copies up to [length] bytes starting at [position] into [destination] at [destinationOffset].
     *
     * Returns the number of bytes copied. Zero is returned when [length] is zero and [position] is not negative, and
     * when [position] is at or past [size].
     *
     * @throws IllegalArgumentException when [position] or [length] is negative, or the destination range does not fit.
     * @throws IllegalStateException when this source is closed.
     */
    suspend fun readAt(position: Long, destination: ByteArray, destinationOffset: Int, length: Int): Int

    /** Releases the source. A second call does nothing. */
    override fun close()
}

/**
 * A [RandomAccessSource] over bytes already in memory.
 *
 * The array is retained, not copied, and must not be mutated while this source is open. A [ByteArray] cannot hold a
 * file past [Int.MAX_VALUE]; a position at or past [size] returns a short count, including a position past 2^31.
 */
class ByteArrayRandomAccessSource(bytes: ByteArray) : RandomAccessSource by InMemoryRandomAccessSource(bytes)

private class InMemoryRandomAccessSource(private val bytes: ByteArray) : BaseRandomAccessSource(bytes.size.toLong()) {
    override fun readAtMost(position: Long, destination: ByteArray, destinationOffset: Int, length: Int): Int {
        val start = position.toInt()
        bytes.copyInto(destination, destinationOffset, start, start + length)
        return length
    }
}

/**
 * Shared closed-flag and argument checks. File and memory sources extend this. A browser `Blob` source can implement
 * [RandomAccessSource] directly.
 */
internal abstract class BaseRandomAccessSource(private val cachedSize: Long) : RandomAccessSource {
    private val closed = AtomicBoolean(false)

    final override val size: Long
        get() {
            ensureOpen()
            return cachedSize
        }

    final override suspend fun readAt(
        position: Long,
        destination: ByteArray,
        destinationOffset: Int,
        length: Int,
    ): Int {
        ensureOpen()
        requireReadArguments(position, destination, destinationOffset, length)
        if (length == 0 || position >= cachedSize) return 0
        val count = minOf(length.toLong(), cachedSize - position).toInt()
        return readAtMost(position, destination, destinationOffset, count)
    }

    final override fun close() {
        if (closed.compareAndSet(expectedValue = false, newValue = true)) onClose()
    }

    /**
     * Copies exactly [length] bytes, which the caller has already limited to the bytes that exist at [position].
     * Loops through a short platform read. Returns a shorter count only when the source ends early.
     */
    protected abstract fun readAtMost(position: Long, destination: ByteArray, destinationOffset: Int, length: Int): Int

    protected open fun onClose() {}

    protected fun ensureOpen() {
        if (closed.load()) throw closedRandomAccessSource()
    }
}

internal fun requireReadArguments(position: Long, destination: ByteArray, destinationOffset: Int, length: Int) {
    if (position < 0) throw IllegalArgumentException("position is negative: $position")
    val end = destinationOffset.toLong() + length.toLong()
    if (destinationOffset < 0 || length < 0 || end > destination.size) {
        throw IllegalArgumentException(
            "destination range [$destinationOffset, $end) is outside the array of size ${destination.size}",
        )
    }
}

internal fun closedRandomAccessSource(): IllegalStateException = IllegalStateException("Random access source is closed")
