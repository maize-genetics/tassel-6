package net.maizegenetics.tassel.io.deflate

import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import kotlin.time.TimeSource

class BgzfInflateTest {
    @Test
    fun tutorialBlocksMatchThePlainVcfAndTheirTrailers() {
        val inflated = inflateTutorial()
        assertContentEquals(TestData.readBytes("tutorial/mdp_genotype_w_depth.vcf"), inflated.bytes)
        assertTrue(inflated.blocks > 1)
        assertEquals(0, inflated.lastUncompressed)
    }

    @Test
    fun reportsTutorialBlockThroughput() {
        val gz = TestData.readBytes("htslib/tutorial/mdp_genotype_w_depth.vcf.gz")
        val blocks = bgzfBlocks(gz)
        val output = ByteArray(blocks.sumOf { it.uncompressedSize })
        val decoder = InflateDecoder()
        val mark = TimeSource.Monotonic.markNow()
        var total = 0
        repeat(PASSES) {
            var at = 0
            for (block in blocks) {
                val written = decoder.inflate(block.file, block.from, block.to, output, at, block.uncompressedSize)
                at += written
            }
            total += at
        }
        val nanos = mark.elapsedNow().inWholeNanoseconds
        val seconds = nanos / 1_000_000_000.0
        val megabytes = total / 1_000_000.0
        println("pure Kotlin inflate: ${oneDecimal(megabytes / seconds)} MB/s over $total output bytes")
        assertEquals(output.size * PASSES, total)
    }

    private fun inflateTutorial(): Inflated {
        val gz = TestData.readBytes("htslib/tutorial/mdp_genotype_w_depth.vcf.gz")
        val blocks = bgzfBlocks(gz)
        val output = ByteArray(blocks.sumOf { it.uncompressedSize })
        val decoder = InflateDecoder()
        val crc = Crc32()
        var at = 0
        for (block in blocks) {
            val written = decoder.inflate(block.file, block.from, block.to, output, at, block.uncompressedSize)
            assertEquals(block.uncompressedSize, written)
            crc.reset()
            crc.update(output, at, at + written)
            assertEquals(block.crc, crc.value, "crc at output offset $at")
            at += written
        }
        return Inflated(output, blocks.size, blocks.last().uncompressedSize)
    }

    private fun oneDecimal(value: Double): String {
        val scaled = (value * 10.0).toLong()
        val whole = scaled / 10
        val fraction = scaled % 10
        return "$whole.$fraction"
    }

    private class Inflated(val bytes: ByteArray, val blocks: Int, val lastUncompressed: Int)

    private companion object {
        const val PASSES = 3
    }
}
