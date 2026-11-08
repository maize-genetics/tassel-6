package net.maizegenetics.tassel.io

import java.nio.ByteBuffer
import java.nio.channels.FileChannel
import java.nio.file.Files
import java.nio.file.StandardOpenOption

internal actual fun openRandomAccessFile(absolutePath: String): RandomAccessSource =
    openRandomAccessSource(absolutePath)

internal actual fun createSparseFixture(size: Long, writes: List<Pair<Long, ByteArray>>): String {
    val file = Files.createTempFile("tassel-sparse", ".bin")
    FileChannel.open(file, StandardOpenOption.WRITE, StandardOpenOption.READ).use { channel ->
        if (size > 0) {
            channel.write(ByteBuffer.wrap(byteArrayOf(0)), size - 1)
        }
        for ((offset, bytes) in writes) {
            if (bytes.isNotEmpty()) channel.write(ByteBuffer.wrap(bytes), offset)
        }
    }
    return file.toString()
}

internal actual fun deleteFixture(path: String) {
    Files.deleteIfExists(java.nio.file.Path.of(path))
}
