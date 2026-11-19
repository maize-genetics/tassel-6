package net.maizegenetics.tassel.io

import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.io.deflate.bgzfBlocks
import net.maizegenetics.tassel.testing.TestData
import java.io.ByteArrayInputStream
import java.util.zip.GZIPInputStream
import java.util.zip.Inflater
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/** JDK inflate of BGZF. [Inflater] and [GZIPInputStream] stay out of `commonMain`. */
class BgzfWriterJvmTest {
    @Test
    fun gzipInputStreamInflatesStoredAndCompressedBlocks() {
        val small = TestData.readBytes(SOURCE_VCF)
        assertContentEquals(small, gunzip(framed(small, BgzfCompressor.Stored)))
        assertContentEquals(small, gunzip(framed(small)))

        val longLine = Buffer().also { it.write(TestData.readBytes(LONG_LINE)) }.decompressing().source.readByteArray()
        assertContentEquals(longLine, gunzip(framed(longLine)))

        val tutorial = TestData.readBytes(PLAIN_VCF)
        assertContentEquals(tutorial, gunzip(framed(tutorial)))
    }

    @Test
    fun inflaterAcceptsTutorialDeflatePayloads() {
        val plain = TestData.readBytes(PLAIN_VCF)
        val written = framed(plain)
        val blocks = bgzfBlocks(written)
        assertTrue(blocks.size > 2)
        var offset = 0
        for (block in blocks) {
            if (block.uncompressedSize == 0) continue
            val payload = written.copyOfRange(block.from, block.to)
            val inflated = inflateRaw(payload, block.uncompressedSize)
            assertContentEquals(plain.copyOfRange(offset, offset + inflated.size), inflated)
            offset += inflated.size
        }
        assertEquals(plain.size, offset)
    }

    private fun framed(plain: ByteArray, compressor: BgzfCompressor = BgzfCompressor.deflate()): ByteArray {
        val sink = Buffer()
        BgzfWriter(sink, BgzfLayout.Vcf, compressor).use { it.write(plain) }
        return sink.readByteArray()
    }

    private fun inflateRaw(payload: ByteArray, size: Int): ByteArray {
        val out = ByteArray(size)
        Inflater(true).use { inflater ->
            inflater.setInput(payload)
            var filled = 0
            while (filled < out.size) {
                val count = inflater.inflate(out, filled, out.size - filled)
                if (count == 0) break
                filled += count
            }
            assertEquals(size, filled)
            assertEquals(payload.size.toLong(), inflater.bytesRead, "deflate payload was rejected")
            assertEquals(size.toLong(), inflater.bytesWritten)
        }
        return out
    }

    private fun gunzip(bytes: ByteArray): ByteArray =
        GZIPInputStream(ByteArrayInputStream(bytes)).use { it.readBytes() }

    private companion object {
        const val PLAIN_VCF = "tutorial/mdp_genotype_w_depth.vcf"
        const val LONG_LINE = "htslib/bgzf-layout/long-line.vcf.gz"
        const val SOURCE_VCF = "htslib/variants/source.vcf"
    }
}
