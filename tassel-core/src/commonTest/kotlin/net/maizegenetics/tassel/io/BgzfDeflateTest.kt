package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.io.deflate.bgzfBlocks
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.sha256
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Default DEFLATE BGZF against the level-0 writer (ADR-0022, #131).
 *
 * Block splits stay the stored writer's splits. Compressed bytes are the same on every target. The tutorial
 * file is 983,502 bytes, inside 25 percent of `bgzip -c` (992,493 bytes).
 */
class BgzfDeflateTest {
    @Test
    fun tutorialVcfMatchesTheSizeBarAndHash() = runTest {
        val plain = TestData.readBytes(TUTORIAL_VCF)
        val written = write(plain, BgzfLayout.Vcf)
        assertEquals(TUTORIAL_COMPRESSED_BYTES, written.size)
        assertTrue(written.size <= SIZE_BAR, "compressed tutorial is ${written.size} bytes, bar is $SIZE_BAR")
        val hash = sha256(written)
        assertEquals(COMPRESSED_TUTORIAL_SHA256, hash, hash)
        assertEquals(isizes(write(plain, BgzfLayout.Vcf, BgzfCompressor.Stored)), isizes(written))
        assertContentEquals(plain, readAll(written))
    }

    @Test
    fun sourceVcfBytesDoNotDependOnTheWriteSize() {
        val plain = TestData.readBytes(SOURCE_VCF)
        val whole = write(plain, BgzfLayout.Vcf)
        val sink = Buffer()
        BgzfWriter(sink, BgzfLayout.Vcf).use { writer ->
            for (index in plain.indices) writer.write(plain, index, index + 1)
        }
        assertContentEquals(whole, sink.readByteArray())
    }

    @Test
    fun fixturePayloadsRoundTrip() = runTest {
        val vcf = listOf(
            SOURCE_VCF to false,
            "golden/vcf-export/mdp_genotype.vcf" to false,
            "golden/vcf-export/include-list.vcf" to false,
            "golden/vcf-export/vcf_export_edges.vcf" to false,
            "htslib/bgzf-layout/long-line.vcf.gz" to true,
            "htslib/bgzf-layout/long-header.vcf.gz" to true,
            "htslib/wheat/wheat.vcf.gz" to true,
            "htslib/large-pos/large_pos.vcf.gz" to true,
        )
        for ((path, gzip) in vcf) assertRoundTrip(path, BgzfLayout.Vcf, gzip)
        val bcf = listOf(
            "htslib/bcf/tutorial/mdp_genotype_w_depth.bcf" to true,
            "htslib/bcf/wheat/wheat.bcf" to true,
            "htslib/bcf/small/small.bcf" to true,
            "htslib/bcf/small/small.stored.bcf" to true,
            "htslib/bcf/small/small.raw.bcf" to false,
            "htslib/bcf/boundary/pos2147483647.bcf" to true,
            "htslib/bcf/edge/edge.bcf" to true,
            "htslib/bcf/edge/no-sample.bcf" to true,
            "htslib/bcf/golden/mdp_genotype.bcf" to true,
            "htslib/bcf/golden/include-list.bcf" to true,
            "htslib/bcf/golden/vcf_export_edges.bcf" to true,
            "htslib/bcf/corrupt/bad-magic.bcf" to true,
            "htslib/bcf/corrupt/bcf21.bcf" to true,
            "htslib/bcf/corrupt/truncated-record.bcf" to true,
            "htslib/bcf/corrupt/no-eof.bcf" to true,
        )
        for ((path, gzip) in bcf) assertRoundTrip(path, BgzfLayout.Bcf, gzip)
    }

    private suspend fun assertRoundTrip(path: String, layout: BgzfLayout, gzip: Boolean) {
        val loaded = TestData.readBytes(path)
        val plain = if (gzip) inflate(loaded) else loaded
        val stored = write(plain, layout, BgzfCompressor.Stored)
        val compressed = write(plain, layout)
        assertEquals(isizes(stored), isizes(compressed), path)
        assertContentEquals(plain, readAll(compressed), path)
    }

    private suspend fun readAll(bytes: ByteArray): ByteArray = ByteArrayRandomAccessSource(bytes).use { source ->
        BgzfReader(source).use { reader ->
            val buffer = Buffer()
            val chunk = ByteArray(64 * 1024)
            while (true) {
                val count = reader.read(chunk, 0, chunk.size)
                if (count == 0) break
                buffer.write(chunk, 0, count)
            }
            buffer.readByteArray()
        }
    }

    private fun inflate(bytes: ByteArray): ByteArray =
        Buffer().also { it.write(bytes) }.decompressing().source.readByteArray()

    private fun isizes(file: ByteArray): List<Int> = bgzfBlocks(file).map { it.uncompressedSize }

    private companion object {
        const val TUTORIAL_VCF = "tutorial/mdp_genotype_w_depth.vcf"
        const val SOURCE_VCF = "htslib/variants/source.vcf"

        /** 25 percent above `bgzip -c` of the tutorial VCF (992,493 bytes). */
        const val SIZE_BAR = 1_240_616

        /** Measured size of [TUTORIAL_VCF] at the fixed level. */
        const val TUTORIAL_COMPRESSED_BYTES = 983_502

        /** Same bytes on every target. A mismatch prints the actual hash. */
        const val COMPRESSED_TUTORIAL_SHA256 =
            "ce705a1f19a5b7e67d06edf06ad2d80c993a589fc0724c3c9a4957f52390550f"

        fun write(plain: ByteArray, layout: BgzfLayout, compressor: BgzfCompressor? = null): ByteArray {
            val sink = Buffer()
            // A null compressor is the public default: parallel where the host has workers.
            val writer = if (compressor == null) BgzfWriter(sink, layout) else BgzfWriter(sink, layout, compressor)
            writer.use { it.write(plain) }
            return sink.readByteArray()
        }
    }
}
