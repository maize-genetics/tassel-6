package net.maizegenetics.tassel.io.deflate

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.RawSource
import kotlinx.io.buffered
import kotlinx.io.readByteArray
import kotlinx.io.readString
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.io.Vcf
import net.maizegenetics.tassel.io.decompressing
import net.maizegenetics.tassel.testing.TestData
import java.io.ByteArrayInputStream
import java.io.ByteArrayOutputStream
import java.util.zip.CRC32
import java.util.zip.Deflater
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream
import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Plain gzip against `java.util.zip.GZIPInputStream`, the VCF golden, and a stream too large to buffer.
 */
class GzipStreamJvmTest {
    @Test
    fun matchesGzipInputStreamOnRandomMembers() {
        val random = Random(85)
        val sizes = intArrayOf(0, 1, 3, 100, 4_096, 32_768, 50_000)
        for (size in sizes) {
            val data = ByteArray(size) { random.nextInt().toByte() }
            val gz = gzip(data)
            assertContentEquals(gunzipJdk(gz), gunzipOurs(gz), "size $size")
        }

        val block = ByteArray(32_768) { random.nextInt().toByte() }
        val repeated = ByteArray(block.size * 2 + 100)
        block.copyInto(repeated)
        block.copyInto(repeated, block.size)
        block.copyInto(repeated, block.size * 2, 0, 100)
        assertTrue(repeated.size > 32_768)
        assertContentEquals(gunzipJdk(gzip(repeated)), gunzipOurs(gzip(repeated)), "repeated 32 KiB block")

        val distance = gzipWithDistance32768()
        val distancePlain = gunzipJdk(distance)
        assertEquals(32_768 + 258, distancePlain.size)
        assertContentEquals(distancePlain, gunzipOurs(distance), "distance 32768")

        val second = ByteArray(200) { random.nextInt().toByte() }
        val concatenated = distance + gzip(second)
        assertContentEquals(gunzipJdk(concatenated), gunzipOurs(concatenated), "concatenated members")
    }

    @Test
    fun headerFieldsMatchGzipInputStream() {
        val data = "header-crc".encodeToByteArray()
        val gz = gzipWithHeader(data, name = "sample.vcf", comment = "note", extra = byteArrayOf(9, 8, 7))
        assertContentEquals(data, gunzipJdk(gz))
        assertContentEquals(data, gunzipOurs(gz))
    }

    @Test
    fun gzipOutputStreamTutorialMatchesTheVcfGolden() = runTest {
        val gz = gzip(TestData.readBytes(TUTORIAL))
        val opened = Buffer().also { it.write(gz) }.decompressing()
        val written = Buffer().also { HapMap.write(Vcf.read(opened.source), it) }.readString().lines()
        assertEquals(emptyList(), opened.warnings, "plain gzip does not warn about a BGZF EOF block")

        val goldenLines = TestData.readText(GOLDEN).replace("\r\n", "\n").lines()
        val skipped = goldenLines.takeWhile { it.startsWith("##") }
        assertEquals(281, skipped.size, "skipped ##SAMPLE= lines")
        assertTrue(skipped.all { it.startsWith("##SAMPLE=<ID=") }, "every skipped line is a ##SAMPLE= annotation")
        val golden = goldenLines.drop(skipped.size)
        assertTrue(golden.first().startsWith("rs#\talleles\tchrom\tpos\t"), "first compared line")
        assertEquals(golden.size, written.size, "number of lines")
        written.forEachIndexed { index, line ->
            if (line != golden[index]) {
                val column = golden[index].commonPrefixWith(line).length
                throw AssertionError(
                    "Line ${skipped.size + index + 1} differs at column ${column + 1}:\n" +
                        "  expected: ...${golden[index].drop(column).take(40)}\n" +
                        "  actual:   ...${line.drop(column).take(40)}",
                )
            }
        }
    }

    @Test
    fun storedZeroMemberInflates() {
        val total = 100_000
        val bytes = StoredZeroGzip(total.toLong()).buffered().decompressing().source.readByteArray()
        assertEquals(total, bytes.size)
        assertTrue(bytes.all { it == 0.toByte() })
    }

    @Test
    fun hundredsOfMegabytesDoNotStayOnTheHeap() {
        val total = 256L * 1024 * 1024
        System.gc()
        val before = usedHeap()
        var read = 0L
        var mixed = 0
        run {
            val opened = StoredZeroGzip(total).buffered().decompressing()
            val scratch = Buffer()
            val chunk = ByteArray(8 * 1024)
            while (true) {
                val n = opened.source.readAtMostTo(scratch, chunk.size.toLong())
                if (n <= 0L) break
                val got = scratch.readAtMostTo(chunk, 0, n.toInt())
                var index = 0
                while (index < got) {
                    mixed = mixed or chunk[index].toInt()
                    index++
                }
                read += got
            }
            assertEquals(emptyList(), opened.warnings)
        }
        assertEquals(total, read)
        assertEquals(0, mixed)
        System.gc()
        val delta = usedHeap() - before
        println("plain gzip 256 MiB stored stream retained $delta bytes after GC")
        assertTrue(delta < 64L * 1024 * 1024, "heap grew by $delta bytes")
    }

    @Test
    fun reportsPlainGzipAndBgzfOneShotThroughput() {
        val plain = TestData.readBytes(TUTORIAL)
        val gzipped = gzip(plain)
        val bgzf = TestData.readBytes(TUTORIAL_GZ)
        val blocks = bgzfBlocks(bgzf)
        val plainBytes = timePlain(gzipped)
        val bgzfBytes = timeBgzf(blocks)
        val plainSeconds = plainBytes.second / 1_000_000_000.0
        val bgzfSeconds = bgzfBytes.second / 1_000_000_000.0
        println(
            "plain gzip inflate: ${oneDecimal(plainBytes.first / 1_000_000.0 / plainSeconds)} MB/s " +
                "over ${plainBytes.first} output bytes",
        )
        println(
            "bgzf one-shot inflate: ${oneDecimal(bgzfBytes.first / 1_000_000.0 / bgzfSeconds)} MB/s " +
                "over ${bgzfBytes.first} output bytes",
        )
        assertEquals(plain.size.toLong() * PASSES, plainBytes.first)
        assertTrue(bgzfBytes.first > 0L)
    }

    private fun timePlain(gzipped: ByteArray): Pair<Long, Long> {
        val scratch = Buffer()
        val chunk = ByteArray(64 * 1024)
        val mark = System.nanoTime()
        var total = 0L
        repeat(PASSES) {
            val opened = Buffer().also { it.write(gzipped) }.decompressing()
            while (true) {
                val n = opened.source.readAtMostTo(scratch, chunk.size.toLong())
                if (n <= 0L) break
                scratch.readAtMostTo(chunk, 0, n.toInt())
                total += n
            }
        }
        return total to (System.nanoTime() - mark)
    }

    private fun timeBgzf(blocks: List<BgzfBlock>): Pair<Long, Long> {
        val output = ByteArray(65_536)
        val decoder = InflateDecoder()
        val mark = System.nanoTime()
        var total = 0L
        repeat(PASSES) {
            for (block in blocks) {
                total += decoder.inflate(block.file, block.from, block.to, output, 0, block.uncompressedSize)
            }
        }
        return total to (System.nanoTime() - mark)
    }

    private fun usedHeap(): Long {
        val runtime = Runtime.getRuntime()
        return runtime.totalMemory() - runtime.freeMemory()
    }

    private fun oneDecimal(value: Double): String {
        val scaled = (value * 10.0).toLong()
        val whole = scaled / 10
        val fraction = kotlin.math.abs(scaled % 10)
        return "$whole.$fraction"
    }

    private fun gzip(data: ByteArray): ByteArray {
        val out = ByteArrayOutputStream()
        GZIPOutputStream(out).use { it.write(data) }
        return out.toByteArray()
    }

    private fun gunzipJdk(gz: ByteArray): ByteArray = GZIPInputStream(ByteArrayInputStream(gz)).use { it.readBytes() }

    private fun gunzipOurs(gz: ByteArray): ByteArray =
        Buffer().also { it.write(gz) }.decompressing().source.readByteArray()

    /**
     * One plain member whose deflate body is stored blocks of zeros, generated as it is read.
     * The uncompressed bytes are never retained.
     */
    private class StoredZeroGzip(private val total: Long) : RawSource {
        private val crc = Crc32()
        private val zeros = ByteArray(8 * 1024)
        private val header = byteArrayOf(0x1f, 0x8b.toByte(), 8, 0, 0, 0, 0, 0, 0, 0)
        private val blockHeader = ByteArray(5)
        private val trailer = ByteArray(8)
        private var stage = HEADER
        private var headerPos = 0
        private var blockHeaderPos = 0
        private var blockRemaining = 0
        private var produced = 0L
        private var trailerPos = 0
        private var trailerReady = false

        override fun readAtMostTo(sink: Buffer, byteCount: Long): Long {
            if (byteCount == 0L) return 0L
            if (stage == DONE) return -1L
            var written = 0L
            while (written < byteCount && stage != DONE) {
                val before = stage
                val n = when (stage) {
                    HEADER -> copy(sink, header, headerPos, header.size, byteCount - written).also { headerPos += it }
                    BLOCK_HEADER -> copy(sink, blockHeader, blockHeaderPos, 5, byteCount - written).also {
                        blockHeaderPos += it
                    }
                    BLOCK_DATA -> writeZeros(sink, byteCount - written)
                    TRAILER -> writeTrailer(sink, byteCount - written)
                    else -> 0
                }
                if (stage == HEADER && headerPos == header.size) prepareBlock()
                if (stage == BLOCK_HEADER && blockHeaderPos == 5) stage = BLOCK_DATA
                if (n == 0 && stage == before) break
                written += n
            }
            return if (written == 0L && stage == DONE) -1L else written
        }

        override fun close() = Unit

        private fun writeZeros(sink: Buffer, space: Long): Int {
            val n = minOf(space, blockRemaining.toLong(), zeros.size.toLong()).toInt()
            if (n == 0) return 0
            sink.write(zeros, 0, n)
            crc.update(zeros, 0, n)
            blockRemaining -= n
            produced += n
            if (blockRemaining == 0) {
                if (produced == total) stage = TRAILER else prepareBlock()
            }
            return n
        }

        private fun writeTrailer(sink: Buffer, space: Long): Int {
            if (!trailerReady) {
                val checksum = crc.value
                val size = total.toInt()
                putU32(trailer, 0, checksum)
                putU32(trailer, 4, size)
                trailerReady = true
            }
            val n = copy(sink, trailer, trailerPos, 8, space)
            trailerPos += n
            if (trailerPos == 8) stage = DONE
            return n
        }

        private fun prepareBlock() {
            val size = minOf(65_535L, total - produced).toInt()
            val last = produced + size == total
            blockHeader[0] = if (last) 1 else 0
            blockHeader[1] = (size and 0xFF).toByte()
            blockHeader[2] = ((size ushr 8) and 0xFF).toByte()
            val complement = size xor 0xFFFF
            blockHeader[3] = (complement and 0xFF).toByte()
            blockHeader[4] = ((complement ushr 8) and 0xFF).toByte()
            blockHeaderPos = 0
            blockRemaining = size
            stage = BLOCK_HEADER
        }

        private fun copy(sink: Buffer, from: ByteArray, pos: Int, end: Int, space: Long): Int {
            val n = minOf(space, (end - pos).toLong()).toInt()
            if (n <= 0) return 0
            sink.write(from, pos, pos + n)
            return n
        }

        private fun putU32(dest: ByteArray, at: Int, value: Int) {
            dest[at] = (value and 0xFF).toByte()
            dest[at + 1] = ((value ushr 8) and 0xFF).toByte()
            dest[at + 2] = ((value ushr 16) and 0xFF).toByte()
            dest[at + 3] = ((value ushr 24) and 0xFF).toByte()
        }

        private companion object {
            const val HEADER = 0
            const val BLOCK_HEADER = 1
            const val BLOCK_DATA = 2
            const val TRAILER = 3
            const val DONE = 4
        }
    }

    private fun gzipWithDistance32768(): ByteArray {
        val lead = ByteArray(32_768) { ((it * 31) + 7).toByte() }
        val writer = BitWriter()
        writeStored(writer, lead, last = false)
        writeFixed(writer, ByteArray(0), last = true, match = Match(258, 32_768))
        val plain = ByteArray(lead.size + 258)
        lead.copyInto(plain)
        lead.copyInto(plain, lead.size, 0, 258)
        return wrapGzip(writer.toByteArray(), plain)
    }

    private fun wrapGzip(deflated: ByteArray, plain: ByteArray): ByteArray {
        val header = byteArrayOf(0x1f, 0x8b.toByte(), 8, 0, 0, 0, 0, 0, 0, 0)
        val dataCrc = CRC32().also { it.update(plain) }.value.toInt()
        val out = ByteArray(header.size + deflated.size + 8)
        header.copyInto(out)
        deflated.copyInto(out, header.size)
        val trailer = out.size - 8
        out[trailer] = (dataCrc and 0xFF).toByte()
        out[trailer + 1] = ((dataCrc ushr 8) and 0xFF).toByte()
        out[trailer + 2] = ((dataCrc ushr 16) and 0xFF).toByte()
        out[trailer + 3] = ((dataCrc ushr 24) and 0xFF).toByte()
        out[trailer + 4] = (plain.size and 0xFF).toByte()
        out[trailer + 5] = ((plain.size ushr 8) and 0xFF).toByte()
        out[trailer + 6] = ((plain.size ushr 16) and 0xFF).toByte()
        out[trailer + 7] = ((plain.size ushr 24) and 0xFF).toByte()
        return out
    }

    private fun gzipWithHeader(data: ByteArray, name: String, comment: String, extra: ByteArray): ByteArray {
        val deflated = ByteArrayOutputStream()
        Deflater(Deflater.DEFAULT_COMPRESSION, true).use { deflater ->
            deflater.setInput(data)
            deflater.finish()
            val buffer = ByteArray(256)
            while (!deflater.finished()) {
                val n = deflater.deflate(buffer)
                deflated.write(buffer, 0, n)
            }
        }
        val nameBytes = name.encodeToByteArray() + 0
        val commentBytes = comment.encodeToByteArray() + 0
        val header = ByteArray(10 + 2 + extra.size + nameBytes.size + commentBytes.size)
        header[0] = 0x1f
        header[1] = 0x8b.toByte()
        header[2] = 8
        header[3] = (4 or 8 or 16 or 2).toByte()
        var cursor = 10
        header[cursor] = (extra.size and 0xFF).toByte()
        header[cursor + 1] = ((extra.size ushr 8) and 0xFF).toByte()
        extra.copyInto(header, cursor + 2)
        cursor += 2 + extra.size
        nameBytes.copyInto(header, cursor)
        cursor += nameBytes.size
        commentBytes.copyInto(header, cursor)
        val checksum = CRC32().also { it.update(header) }
        val crc16 = (checksum.value and 0xFFFF).toInt()
        val payload = deflated.toByteArray()
        val dataCrc = CRC32().also { it.update(data) }.value.toInt()
        val out = ByteArray(header.size + 2 + payload.size + 8)
        header.copyInto(out)
        out[header.size] = (crc16 and 0xFF).toByte()
        out[header.size + 1] = ((crc16 ushr 8) and 0xFF).toByte()
        payload.copyInto(out, header.size + 2)
        val trailer = out.size - 8
        out[trailer] = (dataCrc and 0xFF).toByte()
        out[trailer + 1] = ((dataCrc ushr 8) and 0xFF).toByte()
        out[trailer + 2] = ((dataCrc ushr 16) and 0xFF).toByte()
        out[trailer + 3] = ((dataCrc ushr 24) and 0xFF).toByte()
        out[trailer + 4] = (data.size and 0xFF).toByte()
        out[trailer + 5] = ((data.size ushr 8) and 0xFF).toByte()
        out[trailer + 6] = ((data.size ushr 16) and 0xFF).toByte()
        out[trailer + 7] = ((data.size ushr 24) and 0xFF).toByte()
        return out
    }

    private companion object {
        const val PASSES = 3
        const val TUTORIAL = "tutorial/mdp_genotype_w_depth.vcf"
        const val TUTORIAL_GZ = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
        const val GOLDEN = "golden/vcf/mdp_genotype_w_depth.hmp.txt"
    }
}
