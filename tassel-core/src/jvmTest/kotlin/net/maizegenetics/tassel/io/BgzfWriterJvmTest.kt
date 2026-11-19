package net.maizegenetics.tassel.io

import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.testing.TestData
import java.io.ByteArrayInputStream
import java.util.zip.GZIPInputStream
import kotlin.test.Test
import kotlin.test.assertContentEquals

/** JDK inflate of stored BGZF. [GZIPInputStream] stays out of `commonMain`. */
class BgzfWriterJvmTest {
    @Test
    fun gzipInputStreamInflatesStoredBlocks() {
        val small = TestData.readBytes(SOURCE_VCF)
        assertContentEquals(small, gunzip(stored(small)))

        val longLine = Buffer().also { it.write(TestData.readBytes(LONG_LINE)) }.decompressing().source.readByteArray()
        assertContentEquals(longLine, gunzip(stored(longLine)))

        val tutorial = TestData.readBytes(PLAIN_VCF)
        assertContentEquals(tutorial, gunzip(stored(tutorial)))
    }

    private fun stored(plain: ByteArray): ByteArray {
        val sink = Buffer()
        BgzfWriter(sink, BgzfLayout.Vcf).use { it.write(plain) }
        return sink.readByteArray()
    }

    private fun gunzip(bytes: ByteArray): ByteArray =
        GZIPInputStream(ByteArrayInputStream(bytes)).use { it.readBytes() }

    private companion object {
        const val PLAIN_VCF = "tutorial/mdp_genotype_w_depth.vcf"
        const val LONG_LINE = "htslib/bgzf-layout/long-line.vcf.gz"
        const val SOURCE_VCF = "htslib/variants/source.vcf"
    }
}
