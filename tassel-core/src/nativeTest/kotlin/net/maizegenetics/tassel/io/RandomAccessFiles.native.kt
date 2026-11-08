@file:OptIn(ExperimentalForeignApi::class)

package net.maizegenetics.tassel.io

import kotlinx.cinterop.ExperimentalForeignApi
import kotlinx.cinterop.addressOf
import kotlinx.cinterop.toKString
import kotlinx.cinterop.usePinned
import platform.posix.EINTR
import platform.posix.O_CREAT
import platform.posix.O_RDWR
import platform.posix.O_TRUNC
import platform.posix.close
import platform.posix.errno
import platform.posix.ftruncate
import platform.posix.getenv
import platform.posix.open
import platform.posix.pwrite
import platform.posix.remove
import kotlin.random.Random

internal actual fun openRandomAccessFile(absolutePath: String): RandomAccessSource =
    openRandomAccessSource(absolutePath)

internal actual fun createSparseFixture(size: Long, writes: List<Pair<Long, ByteArray>>): String {
    val dir = getenv("TMPDIR")?.toKString()?.trimEnd('/') ?: "/tmp"
    val path = "$dir/tassel-sparse-${Random.nextLong().toULong().toString(16)}"
    val fd = open(path, O_RDWR or O_CREAT or O_TRUNC, MODE)
    if (fd < 0) error("Cannot create $path (errno $errno)")
    try {
        if (ftruncate(fd, size) != 0) error("Cannot size $path (errno $errno)")
        for ((offset, bytes) in writes) writeAll(fd, bytes, offset)
    } finally {
        close(fd)
    }
    return path
}

internal actual fun deleteFixture(path: String) {
    remove(path)
}

private fun writeAll(fd: Int, bytes: ByteArray, offset: Long) {
    if (bytes.isEmpty()) return
    var written = 0
    bytes.usePinned { pinned ->
        while (written < bytes.size) {
            val n = pwrite(
                fd,
                pinned.addressOf(written),
                (bytes.size - written).toULong(),
                offset + written,
            )
            if (n < 0) {
                if (errno == EINTR) continue
                error("Cannot write sparse fixture (errno $errno)")
            }
            if (n == 0L) error("Short write of sparse fixture")
            written += n.toInt()
        }
    }
}

private const val MODE: UInt = 420u
