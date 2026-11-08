package net.maizegenetics.tassel.io

import kotlinx.coroutines.runBlocking
import kotlinx.io.buffered
import kotlinx.io.files.Path
import kotlinx.io.files.SystemFileSystem
import kotlinx.io.readAtMostTo
import java.nio.ByteBuffer
import java.nio.channels.FileChannel
import java.nio.file.Files
import java.nio.file.StandardOpenOption
import kotlin.test.Test
import kotlin.time.Duration
import kotlin.time.measureTime

/**
 * Sequential positional reads versus [SystemFileSystem.source], for the issue #86 handoff.
 *
 * Stays out of CI. Run with `TASSEL_BENCHMARK=1`.
 */
class RandomAccessBenchmarkTest {
    @Test
    fun sequentialPositionalThroughput() {
        if (System.getenv("TASSEL_BENCHMARK") != "1") return
        val file = Files.createTempFile("tassel-bench", ".bin")
        try {
            val pattern = ByteArray(CHUNK) { (it and 0xFF).toByte() }
            FileChannel.open(
                file,
                StandardOpenOption.WRITE,
                StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING,
            ).use { channel ->
                val buffer = ByteBuffer.wrap(pattern)
                var left = BYTES
                while (left > 0) {
                    buffer.clear()
                    val n = minOf(left, CHUNK.toLong()).toInt()
                    buffer.limit(n)
                    while (buffer.hasRemaining()) channel.write(buffer)
                    left -= n
                }
            }
            val path = file.toString()
            val positional = measureTime {
                runBlocking {
                    openRandomAccessSource(path).use { source ->
                        var pos = 0L
                        val dest = ByteArray(CHUNK)
                        while (pos < BYTES) {
                            val n = source.readAt(pos, dest, 0, CHUNK)
                            check(n > 0) { "short read at $pos" }
                            pos += n
                        }
                        check(pos == BYTES)
                    }
                }
            }
            val sequential = measureTime {
                SystemFileSystem.source(Path(path)).buffered().use { source ->
                    var total = 0L
                    val dest = ByteArray(CHUNK)
                    while (true) {
                        val n = source.readAtMostTo(dest)
                        if (n < 0) break
                        total += n
                    }
                    check(total == BYTES)
                }
            }
            println("positional MiB/s=${mibPerSecond(BYTES, positional)}")
            println("SystemFileSystem.source MiB/s=${mibPerSecond(BYTES, sequential)}")
        } finally {
            Files.deleteIfExists(file)
        }
    }

    private fun mibPerSecond(bytes: Long, time: Duration): String {
        val seconds = time.inWholeNanoseconds.toDouble() / 1_000_000_000.0
        val mib = bytes.toDouble() / (1024.0 * 1024.0) / seconds
        return mib.toString()
    }

    private companion object {
        const val CHUNK = 1024 * 1024
        const val BYTES = 128L * 1024 * 1024
    }
}
