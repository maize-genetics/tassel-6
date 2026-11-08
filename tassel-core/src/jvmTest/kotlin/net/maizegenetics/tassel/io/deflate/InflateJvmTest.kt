package net.maizegenetics.tassel.io.deflate

import net.maizegenetics.tassel.testing.TestData
import java.io.ByteArrayOutputStream
import java.util.zip.CRC32
import java.util.zip.Deflater
import java.util.zip.Inflater
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Raw DEFLATE and CRC-32 against the JDK. `nowrap` is required: the default wrapper is zlib, not RFC 1951.
 */
class InflateJvmTest {
    @Test
    fun crcMatchesTheJdk() {
        val random = Random(20261108)
        val buffers = listOf(
            ByteArray(0),
            byteArrayOf(0),
            "123456789".encodeToByteArray(),
            ByteArray(7) { it.toByte() },
            ByteArray(8) { it.toByte() },
            ByteArray(9) { it.toByte() },
            ByteArray(64_000) { random.nextInt().toByte() },
        )
        for (buffer in buffers) {
            val jdk = CRC32()
            jdk.update(buffer)
            val ours = Crc32()
            ours.update(buffer, 0, buffer.size)
            assertEquals(jdk.value.toInt(), ours.value, "length ${buffer.size}")
        }
    }

    @Test
    fun deflaterStreamsDecodeLikeTheJdkInflater() {
        val inputs = sampleInputs()
        val levels = 0..9
        val strategies = intArrayOf(Deflater.DEFAULT_STRATEGY, Deflater.FILTERED, Deflater.HUFFMAN_ONLY)
        val decoder = InflateDecoder()
        for (input in inputs) {
            for (level in levels) {
                for (strategy in strategies) {
                    val compressed = deflate(input, level, strategy)
                    val expected = inflateWithJdk(compressed, input.size)
                    val actual = ByteArray(input.size)
                    val written = decoder.inflate(compressed, 0, compressed.size, actual, 0, actual.size)
                    assertEquals(input.size, written, "level $level strategy $strategy size ${input.size}")
                    assertContentEquals(expected, actual)
                }
            }
        }
    }

    @Test
    fun reportsJdkInflaterThroughputOnTheTutorialBlocks() {
        val gz = TestData.readBytes("htslib/tutorial/mdp_genotype_w_depth.vcf.gz")
        val blocks = bgzfBlocks(gz)
        val output = ByteArray(65_536)
        val inflater = Inflater(true)
        val mark = System.nanoTime()
        var total = 0L
        try {
            repeat(PASSES) {
                for (block in blocks) {
                    inflater.reset()
                    inflater.setInput(block.file, block.from, block.to - block.from)
                    var written = 0
                    while (!inflater.finished()) {
                        val n = inflater.inflate(output)
                        if (n == 0 && inflater.needsInput()) break
                        written += n
                    }
                    total += written
                }
            }
        } finally {
            inflater.end()
        }
        val seconds = (System.nanoTime() - mark) / 1_000_000_000.0
        val megabytes = total / 1_000_000.0
        val scaled = (megabytes / seconds * 10.0).toLong()
        println("java.util.zip.Inflater: ${scaled / 10}.${scaled % 10} MB/s over $total output bytes")
        assertTrue(total > 0L)
    }

    private fun sampleInputs(): List<ByteArray> {
        val random = Random(20261108)
        return listOf(
            ByteArray(0),
            byteArrayOf(0),
            ByteArray(258) { 'a'.code.toByte() },
            ByteArray(3_000) { (it and 0xFF).toByte() },
            ByteArray(32_768) { ((it * 31) + 7).toByte() },
            ByteArray(65_535) { random.nextInt().toByte() },
            ByteArray(4_000),
            ByteArray(8_000) { random.nextInt().toByte() },
        )
    }

    private fun deflate(input: ByteArray, level: Int, strategy: Int): ByteArray {
        val sink = ByteArrayOutputStream()
        Deflater(level, true).use { deflater ->
            deflater.setStrategy(strategy)
            deflater.setInput(input)
            deflater.finish()
            val buffer = ByteArray(8_192)
            while (!deflater.finished()) {
                val written = deflater.deflate(buffer)
                sink.write(buffer, 0, written)
            }
        }
        return sink.toByteArray()
    }

    private fun inflateWithJdk(compressed: ByteArray, size: Int): ByteArray {
        val output = ByteArray(size)
        Inflater(true).use { inflater ->
            inflater.setInput(compressed)
            var at = 0
            while (!inflater.finished()) {
                val written = inflater.inflate(output, at, output.size - at)
                if (written == 0 && inflater.needsInput()) break
                at += written
            }
            assertEquals(size, at)
            assertTrue(inflater.finished())
        }
        return output
    }

    private companion object {
        const val PASSES = 3
    }
}
