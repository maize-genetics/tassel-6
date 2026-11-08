@file:OptIn(ExperimentalWasmJsInterop::class, ExperimentalAtomicApi::class)

package net.maizegenetics.tassel.io

import kotlinx.io.IOException
import kotlin.concurrent.atomics.ExperimentalAtomicApi

/**
 * Opens [path] for positional reads on Node.
 *
 * Bytes cross the JavaScript boundary once per slice, as a `latin1` string, and the file position is a JS number.
 * That number is exact through 2^53, which covers the 3 GiB sparse-file test. There is no `fs` implementation in
 * `wasmJsMain`, so this function is not in the browser bundle (ADR-0021). Calling it where `process` is missing
 * fails with a message that the browser has no file paths.
 */
fun openRandomAccessSource(path: String): RandomAccessSource {
    if (!nodeProcessAvailable()) {
        throw IOException(
            "Opening a file by path needs Node. The browser has no file paths; use an in-memory RandomAccessSource.",
        )
    }
    val fd = try {
        nodeOpen(path)
    } catch (e: Throwable) {
        throw IOException("Cannot open $path", e)
    }
    val size = try {
        nodeFileSize(fd).toLong()
    } catch (e: Throwable) {
        nodeClose(fd)
        throw IOException("Cannot size $path", e)
    }
    return NodeRandomAccessSource(fd, size)
}

private class NodeRandomAccessSource(private val fd: Int, size: Long) : BaseRandomAccessSource(size) {
    override fun readAtMost(position: Long, destination: ByteArray, destinationOffset: Int, length: Int): Int {
        var read = 0
        while (read < length) {
            ensureOpen()
            val chunk = try {
                nodeRead(fd, length - read, (position + read).toDouble())
            } catch (e: Throwable) {
                if (isClosedQuietly()) throw closedRandomAccessSource()
                throw IOException("Cannot read random access source", e)
            }
            if (chunk.isEmpty()) break
            for (i in chunk.indices) {
                destination[destinationOffset + read + i] = chunk[i].code.toByte()
            }
            read += chunk.length
        }
        return read
    }

    override fun onClose() {
        try {
            nodeClose(fd)
        } catch (_: Throwable) {
            // The source is already closed from the caller's point of view.
        }
    }

    private fun isClosedQuietly(): Boolean = try {
        ensureOpen()
        false
    } catch (_: IllegalStateException) {
        true
    }
}

private fun nodeProcessAvailable(): Boolean = js(
    "typeof globalThis.process !== 'undefined' && typeof globalThis.process.getBuiltinModule === 'function'",
)

private fun nodeOpen(path: String): Int = js(
    """(() => {
        const fs = globalThis.process.getBuiltinModule('node:fs');
        try {
            return fs.openSync(path, 'r');
        } catch (error) {
            const message = (error && error.message) ? String(error.message) : String(error);
            throw new Error(message);
        }
        })()""",
)

private fun nodeFileSize(fd: Int): Double = js("globalThis.process.getBuiltinModule('node:fs').fstatSync(fd).size")

private fun nodeRead(fd: Int, length: Int, position: Double): String = js(
    """(() => {
        const fs = globalThis.process.getBuiltinModule('node:fs');
        const Buffer = globalThis.process.getBuiltinModule('node:buffer').Buffer;
        const buf = Buffer.allocUnsafe(length);
        const n = fs.readSync(fd, buf, 0, length, position);
        return buf.toString('latin1', 0, n);
        })()""",
)

private fun nodeClose(fd: Int) {
    js("globalThis.process.getBuiltinModule('node:fs').closeSync(fd)")
}
