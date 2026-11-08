package net.maizegenetics.tassel.io

import kotlinx.io.IOException
import java.nio.ByteBuffer
import java.nio.channels.ClosedChannelException
import java.nio.channels.FileChannel
import java.nio.file.Path
import java.nio.file.StandardOpenOption
import java.io.IOException as JavaIOException

/**
 * Opens [path] for concurrent positional reads.
 *
 * Uses [FileChannel.read] at a position, which does not move a shared offset and is safe for overlapping reads.
 * A missing file throws [IOException]. This is the file analogue of [kotlinx.io.files.SystemFileSystem.source].
 * Readers still take a [RandomAccessSource], never a path (ADR-0021).
 */
fun openRandomAccessSource(path: String): RandomAccessSource {
    val channel = try {
        FileChannel.open(Path.of(path), StandardOpenOption.READ)
    } catch (e: JavaIOException) {
        throw IOException("Cannot open $path", e)
    }
    val size = try {
        channel.size()
    } catch (e: JavaIOException) {
        channel.close()
        throw IOException("Cannot size $path", e)
    }
    return JvmRandomAccessSource(channel, size)
}

private class JvmRandomAccessSource(private val channel: FileChannel, size: Long) : BaseRandomAccessSource(size) {
    override fun readAtMost(position: Long, destination: ByteArray, destinationOffset: Int, length: Int): Int {
        val buffer = ByteBuffer.wrap(destination, destinationOffset, length)
        var read = 0
        try {
            while (buffer.hasRemaining()) {
                ensureOpen()
                val n = channel.read(buffer, position + read)
                if (n <= 0) break
                read += n
            }
        } catch (e: ClosedChannelException) {
            throw closedRandomAccessSource()
        } catch (e: JavaIOException) {
            if (isClosedQuietly()) throw closedRandomAccessSource()
            throw IOException("Cannot read random access source", e)
        }
        return read
    }

    override fun onClose() {
        try {
            channel.close()
        } catch (e: JavaIOException) {
            throw IOException("Cannot close random access source", e)
        }
    }

    private fun isClosedQuietly(): Boolean = try {
        ensureOpen()
        false
    } catch (_: IllegalStateException) {
        true
    }
}
