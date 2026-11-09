@file:OptIn(ExperimentalForeignApi::class)

package net.maizegenetics.tassel.io

import kotlinx.cinterop.ExperimentalForeignApi
import kotlinx.cinterop.addressOf
import kotlinx.cinterop.alloc
import kotlinx.cinterop.memScoped
import kotlinx.cinterop.ptr
import kotlinx.cinterop.toKString
import kotlinx.cinterop.usePinned
import kotlinx.io.IOException
import platform.posix.EINTR
import platform.posix.O_RDONLY
import platform.posix.close
import platform.posix.errno
import platform.posix.fstat
import platform.posix.open
import platform.posix.pread
import platform.posix.stat
import platform.posix.strerror

/**
 * Opens [path] for concurrent positional reads.
 *
 * Uses `pread`, which does not move a shared file offset, so overlapping reads are safe. The offset is a 64-bit
 * `off_t` on macosArm64 and linuxX64. A missing file throws [IOException]. Readers still take a [RandomAccessSource],
 * never a path (ADR-0021).
 */
fun openRandomAccessSource(path: String): RandomAccessSource {
    val fd = open(path, O_RDONLY)
    if (fd < 0) throw IOException("Cannot open $path: ${posixError()}")
    try {
        return NativeRandomAccessSource(fd, fileSize(fd))
    } catch (e: Throwable) {
        close(fd)
        throw e
    }
}

private class NativeRandomAccessSource(private val fd: Int, size: Long) : BaseRandomAccessSource(size) {
    override fun readAtMost(position: Long, destination: ByteArray, destinationOffset: Int, length: Int): Int {
        var read = 0
        destination.usePinned { pinned ->
            while (read < length) {
                ensureOpen()
                val n = pread(
                    fd,
                    pinned.addressOf(destinationOffset + read),
                    (length - read).toULong(),
                    position + read,
                )
                if (n < 0) {
                    val err = errno
                    if (err == EINTR) continue
                    if (isClosedQuietly()) throw closedRandomAccessSource()
                    throw IOException("Cannot read random access source: ${posixError(err)}")
                }
                if (n == 0L) break
                read += n.toInt()
            }
        }
        return read
    }

    override fun onClose() {
        while (close(fd) != 0 && errno == EINTR) {
            // close was interrupted; the descriptor is still open
        }
    }

    private fun isClosedQuietly(): Boolean = try {
        ensureOpen()
        false
    } catch (_: IllegalStateException) {
        true
    }
}

private fun fileSize(fd: Int): Long = memScoped {
    val info = alloc<stat>()
    if (fstat(fd, info.ptr) != 0) throw IOException("Cannot size random access source: ${posixError()}")
    info.st_size
}

private fun posixError(err: Int = errno): String = strerror(err)?.toKString() ?: err.toString()
