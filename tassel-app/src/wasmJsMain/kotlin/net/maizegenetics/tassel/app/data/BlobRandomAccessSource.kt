@file:OptIn(ExperimentalWasmJsInterop::class)

package net.maizegenetics.tassel.app.data

import kotlinx.coroutines.suspendCancellableCoroutine
import net.maizegenetics.tassel.io.RandomAccessSource
import net.maizegenetics.tassel.io.arrayBufferToByteArray
import org.w3c.files.Blob
import org.w3c.files.FileReader
import kotlin.coroutines.resume
import kotlin.coroutines.resumeWithException
import kotlin.js.ExperimentalWasmJsInterop
import kotlin.js.JsAny

/**
 * A [RandomAccessSource] over a browser `Blob` or `File`.
 *
 * `kotlinx-browser`'s `Blob.slice` takes `Int` offsets, so a file past 2 GiB cannot be addressed through it.
 * [sliceBlob] calls `blob.slice` with `Double` offsets instead, which stay exact through 2^53. Each [readAt] uses
 * its own `FileReader`, aborted if the coroutine is cancelled, and copies the slice with [arrayBufferToByteArray].
 * The file itself is never copied whole.
 *
 * Implements [RandomAccessSource] directly. The shared base class reads synchronously, and this read cannot.
 */
internal class BlobRandomAccessSource(private val blob: Blob) : RandomAccessSource {
    private val cachedSize: Long = blobByteSize(blob).toLong()
    private var closed: Boolean = false

    override val size: Long
        get() {
            ensureOpen()
            return cachedSize
        }

    override suspend fun readAt(position: Long, destination: ByteArray, destinationOffset: Int, length: Int): Int {
        ensureOpen()
        if (position < 0) throw IllegalArgumentException("position is negative: $position")
        val end = destinationOffset.toLong() + length.toLong()
        if (destinationOffset < 0 || length < 0 || end > destination.size) {
            throw IllegalArgumentException(
                "destination range [$destinationOffset, $end) is outside the array of size ${destination.size}",
            )
        }
        if (length == 0 || position >= cachedSize) return 0
        val count = minOf(length.toLong(), cachedSize - position).toInt()
        val slice = sliceBlob(blob, position.toDouble(), (position + count).toDouble())
        val bytes = readBlobBytes(slice, "blob")
        val copied = minOf(bytes.size, count)
        if (copied > 0) bytes.copyInto(destination, destinationOffset, 0, copied)
        return copied
    }

    override fun close() {
        closed = true
    }

    private fun ensureOpen() {
        if (closed) throw IllegalStateException("Random access source is closed")
    }
}

/** `blob.size` as a Double. The binding's size is a `JsNumber` and is not wide enough to use as an address. */
private fun blobByteSize(blob: JsAny): Double = js("Number(blob.size)")

/** `blob.slice(start, end)` with Double offsets, past the `Int` limit of the kotlinx-browser binding. */
private fun sliceBlob(blob: JsAny, start: Double, end: Double): Blob = js("blob.slice(start, end)")

/**
 * Reads [blob] as bytes. `readAsText` would replace anything that is not UTF-8.
 *
 * One reader per call. Cancelling the coroutine aborts it, which is what stops a region load between slices.
 */
internal suspend fun readBlobBytes(blob: Blob, name: String): ByteArray = suspendCancellableCoroutine { continuation ->
    val reader = FileReader()
    reader.onload = {
        val buffer = reader.result
        if (buffer == null) {
            continuation.resumeWithException(IllegalStateException("could not read $name"))
        } else {
            continuation.resume(arrayBufferToByteArray(buffer))
        }
    }
    reader.onerror = { continuation.resumeWithException(IllegalStateException("could not read $name")) }
    continuation.invokeOnCancellation { reader.abort() }
    reader.readAsArrayBuffer(blob)
}
