package net.maizegenetics.tassel.testing

import kotlinx.cinterop.ExperimentalForeignApi
import kotlinx.cinterop.addressOf
import kotlinx.cinterop.convert
import kotlinx.cinterop.toKString
import kotlinx.cinterop.usePinned
import platform.posix.SEEK_END
import platform.posix.SEEK_SET
import platform.posix.fclose
import platform.posix.fopen
import platform.posix.fread
import platform.posix.fseek
import platform.posix.ftell
import platform.posix.getenv

@OptIn(ExperimentalForeignApi::class)
internal actual fun environmentVariable(name: String): String? = getenv(name)?.toKString()

@OptIn(ExperimentalForeignApi::class)
internal actual fun readFile(absolutePath: String): ByteArray {
    val file = fopen(absolutePath, "rb") ?: error("Cannot open $absolutePath")
    try {
        check(fseek(file, 0, SEEK_END) == 0) { "Cannot seek in $absolutePath" }
        val size = ftell(file)
        check(size >= 0) { "Cannot size $absolutePath" }
        check(fseek(file, 0, SEEK_SET) == 0) { "Cannot seek in $absolutePath" }
        val bytes = ByteArray(size.toInt())
        if (bytes.isNotEmpty()) {
            val read = bytes.usePinned { fread(it.addressOf(0), 1u, bytes.size.convert(), file) }
            check(read.toLong() == size) { "Short read of $absolutePath" }
        }
        return bytes
    } finally {
        fclose(file)
    }
}
