package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * One BCF record against [VcfRecordParser] on the htslib 1.25 view text (ADR-0022).
 *
 * Sentinels and overflow widths are built here. The edge, boundary, and no-sample files are the #118 fixtures.
 * A missing `QUAL` is the `.` in every edge row, and a missing `FILTER` is the `.` in the second edge row and in
 * the sample-less file.
 */
class BcfRecordTest {
    @Test
    fun integerSentinelsPromoteAndLeaveTheNextByte() {
        val both = byteArrayOf(0x21, 0x80.toByte(), 0x81.toByte(), MARKER)
        val cursor = BcfCursor(both)
        assertEquals(listOf(BCF_INT_MISSING, BCF_INT_VECTOR_END), cursor.readInts().toList())
        assertEquals(3, cursor.offset)
        assertEquals(MARKER.toInt(), both[cursor.offset].toInt())

        assertInt(byteArrayOf(0x11, 0x80.toByte()), BCF_INT_MISSING)
        assertInt(byteArrayOf(0x11, 0x81.toByte()), BCF_INT_VECTOR_END)
        assertInt(byteArrayOf(0x12, 0x00, 0x80.toByte()), BCF_INT_MISSING)
        assertInt(byteArrayOf(0x12, 0x01, 0x80.toByte()), BCF_INT_VECTOR_END)
        assertInt(byteArrayOf(0x13, 0x00, 0x00, 0x00, 0x80.toByte()), BCF_INT_MISSING)
        assertInt(byteArrayOf(0x13, 0x01, 0x00, 0x00, 0x80.toByte()), BCF_INT_VECTOR_END)
    }

    @Test
    fun floatSentinelsKeepTheirBits() {
        val bytes = byteArrayOf(
            0x25,
            0x01, 0x00, 0x80.toByte(), 0x7f,
            0x02, 0x00, 0x80.toByte(), 0x7f,
            MARKER,
        )
        val cursor = BcfCursor(bytes)
        val bits = cursor.readFloatBits()
        assertEquals(BCF_FLOAT_MISSING, bits[0])
        assertEquals(BCF_FLOAT_VECTOR_END, bits[1])
        assertEquals(9, cursor.offset)
        assertEquals(MARKER.toInt(), bytes[cursor.offset].toInt())
    }

    @Test
    fun overflowCountMayBeInt8Int16OrInt32() {
        assertOverflow(20, 0x11)
        assertOverflow(200, 0x12)
        assertOverflow(40_000, 0x13)
    }

    @Test
    fun edgeRecordsMatchTheViewLines() = runTest {
        val rows = viewParts(EDGE_VIEW).second.map { it.split('\t') }
        assertTrue(rows.all { it[5] == "." }, "every edge QUAL is missing")
        assertEquals("PASS", rows[0][6])
        assertEquals(".", rows[1][6])
        assertEquals(3, assertViewParity(EDGE, EDGE_VIEW).size)
    }

    @Test
    fun noSampleRecordMatchesTheViewLine() = runTest {
        val (header, source) = open(NO_SAMPLE)
        val calls = ByteArray(header.vcf.taxaNames.size)
        val position = checkNotNull(readBcfRecord(source, header, calls))
        assertNull(readBcfRecord(source, header, calls))
        val column = viewParts(NO_SAMPLE_VIEW).second.single().split('\t')
        assertEquals(".", column[5])
        assertEquals(".", column[6])
        assertEquals(emptyList<String>(), header.vcf.taxaNames)
        assertEquals(0, calls.size)
        assertEquals(Chromosome.of(column[0]), position.chromosome)
        assertEquals(column[1].toInt(), position.position)
        assertNull(position.snpId)
        val site = "${column[0]}:${column[1]}"
        assertEquals(splitVariants(variantsOf(column[3], column[4]), site), position.knownVariants)
    }

    @Test
    fun boundaryPositionIsAccepted() = runTest {
        val positions = assertViewParity(BOUNDARY, BOUNDARY_VIEW)
        assertEquals(Int.MAX_VALUE, positions.single().position)
    }

    @Test
    fun phasedGenotypeMatchesTheVcfLine() {
        val source = Buffer().also { it.write(bcf(headerText(), record(99, byteArrayOf(0x21, 0x02, 0x05)))) }
        val header = readBcfHeader(source)
        assertEquals(1, header.gtId)
        val calls = ByteArray(header.vcf.taxaNames.size)
        val position = checkNotNull(readBcfRecord(source, header, calls))
        val vcf = readVcfHeader(Buffer().also { it.write(headerText().encodeToByteArray()) })
        val line = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1".encodeToByteArray()
        val expectedCalls = ByteArray(vcf.taxaNames.size)
        val parsed = VcfRecordParser(vcf).tryParse(line, 0, line.size, true, expectedCalls) as VcfParse.Site
        assertSameSite(parsed.position, position, "0|1")
        assertContentEquals(expectedCalls, calls)
    }

    @Test
    fun zeroBasedIntMaxUsesTheVcfPositionMessage() {
        val source = Buffer().also {
            it.write(bcf(headerText(), record(Int.MAX_VALUE, byteArrayOf(0x21, 0x02, 0x04))))
        }
        val header = readBcfHeader(source)
        val failure = assertFailsWith<IllegalArgumentException> {
            readBcfRecord(source, header, ByteArray(header.vcf.taxaNames.size))
        }
        assertEquals(
            "VCF position 2147483648 is past TASSEL's position limit of 2147483647 in chr1:2147483648",
            failure.message,
        )
    }

    @Test
    fun headerWithNoRecordIsTheEnd() {
        val source = Buffer().also { it.write(bcf(headerText(), byteArrayOf())) }
        val header = readBcfHeader(source)
        assertNull(readBcfRecord(source, header, ByteArray(header.vcf.taxaNames.size)))
    }

    @Test
    fun shortRecordIsTruncated() {
        val lengths = LeBytes()
        lengths.int(31)
        lengths.int(0)
        val source = Buffer().also { it.write(bcf(headerText(), lengths.toByteArray())) }
        val header = readBcfHeader(source)
        val failure = assertFailsWith<IllegalArgumentException> {
            readBcfRecord(source, header, ByteArray(header.vcf.taxaNames.size))
        }
        assertEquals(TRUNCATED_BCF_RECORD, failure.message)
    }

    private fun assertInt(payload: ByteArray, expected: Int) {
        val bytes = payload + byteArrayOf(MARKER)
        val cursor = BcfCursor(bytes)
        assertEquals(expected, cursor.readInts().single())
        assertEquals(payload.size, cursor.offset)
        assertEquals(MARKER.toInt(), bytes[cursor.offset].toInt())
    }

    private fun assertOverflow(count: Int, sizeDescriptor: Int) {
        val bytes = overflowChars(count)
        assertEquals(0xf7, bytes[0].toInt() and 0xff)
        assertEquals(sizeDescriptor, bytes[1].toInt() and 0xff)
        val cursor = BcfCursor(bytes)
        assertEquals(count, cursor.readString().length)
        assertEquals(MARKER.toInt(), bytes[cursor.offset].toInt() and 0xff)
        val skipped = BcfCursor(bytes)
        skipped.skip()
        assertEquals(cursor.offset, skipped.offset)
    }

    private suspend fun assertViewParity(bcf: String, view: String): List<Position> {
        val (header, source) = open(bcf)
        val (headerLines, data) = viewParts(view)
        val vcf = readVcfHeaderLines(headerLines.iterator())
        assertEquals(vcf.taxaNames, header.vcf.taxaNames, bcf)
        val parser = VcfRecordParser(vcf)
        val found = ArrayList<Position>()
        for (line in data) {
            val window = line.encodeToByteArray()
            val vcfCalls = ByteArray(vcf.taxaNames.size)
            val parsed = parser.tryParse(window, 0, window.size, true, vcfCalls)
            check(parsed is VcfParse.Site) { "view line was not a site: $line" }
            val calls = ByteArray(header.vcf.taxaNames.size)
            val actual = checkNotNull(readBcfRecord(source, header, calls)) { line }
            assertSameSite(parsed.position, actual, line)
            assertContentEquals(vcfCalls, calls, line)
            found.add(actual)
        }
        assertNull(readBcfRecord(source, header, ByteArray(header.vcf.taxaNames.size)))
        return found
    }

    private fun assertSameSite(expected: Position, actual: Position, label: String) {
        assertEquals(expected.chromosome, actual.chromosome, label)
        assertEquals(expected.position, actual.position, label)
        assertEquals(expected.snpId, actual.snpId, label)
        assertEquals(expected.knownVariants, actual.knownVariants, label)
    }

    private suspend fun open(path: String): Pair<BcfHeader, Buffer> {
        val source = Buffer().also { it.write(inflated(path)) }
        return readBcfHeader(source) to source
    }

    private suspend fun inflated(path: String): ByteArray {
        val file = TestData.readBytes(path)
        if (file.size < 2 || file[0] != GZIP_ID0 || file[1] != GZIP_ID1) return file
        return ByteArrayRandomAccessSource(file).use { random ->
            BgzfReader(random).use { reader ->
                val chunks = ArrayList<ByteArray>()
                var total = 0
                val buffer = ByteArray(8192)
                while (true) {
                    val n = reader.read(buffer, 0, buffer.size)
                    if (n <= 0) break
                    chunks.add(buffer.copyOf(n))
                    total += n
                }
                val out = ByteArray(total)
                var at = 0
                for (chunk in chunks) {
                    chunk.copyInto(out, at)
                    at += chunk.size
                }
                out
            }
        }
    }

    private companion object {
        const val EDGE = "htslib/bcf/edge/edge.bcf"
        const val EDGE_VIEW = "htslib/bcf/edge/edge.view.vcf"
        const val NO_SAMPLE = "htslib/bcf/edge/no-sample.bcf"
        const val NO_SAMPLE_VIEW = "htslib/bcf/edge/no-sample.view.vcf"
        const val BOUNDARY = "htslib/bcf/boundary/pos2147483647.bcf"
        const val BOUNDARY_VIEW = "htslib/bcf/boundary/pos2147483647.view.vcf"
        val GZIP_ID0 = 0x1f.toByte()
        val GZIP_ID1 = 0x8b.toByte()
        val MARKER = 0x5a.toByte()

        fun viewParts(path: String): Pair<List<String>, List<String>> {
            val lines = linesOf(TestData.readText(path))
            val header = lines.takeWhile { it.startsWith("#") }
            return header to lines.drop(header.size)
        }

        fun linesOf(text: String): List<String> {
            val lines = text.split('\n').map { it.trimEnd('\r') }
            return if (lines.isNotEmpty() && lines.last().isEmpty()) lines.dropLast(1) else lines
        }

        fun headerText(): String = buildString {
            append("##fileformat=VCFv4.3\n")
            append("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
            append("##contig=<ID=chr1,length=10>\n")
            append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n")
        }

        fun bcf(header: String, record: ByteArray): ByteArray {
            val text = header.encodeToByteArray()
            val withNul = ByteArray(text.size + 1)
            text.copyInto(withNul)
            val out = LeBytes()
            out.bytes("BCF".encodeToByteArray())
            out.byte(2)
            out.byte(2)
            out.int(withNul.size)
            out.bytes(withNul)
            out.bytes(record)
            return out.toByteArray()
        }

        /** One diploid sample. [gt] is the value descriptor and its bytes, for example `0/1` or `0|1`. */
        fun record(pos0: Int, gt: ByteArray): ByteArray {
            val shared = LeBytes()
            shared.int(0)
            shared.int(pos0)
            shared.int(1)
            shared.int(BCF_FLOAT_MISSING)
            shared.int(2 shl 16)
            shared.int((1 shl 24) or 1)
            shared.byte(0x07)
            shared.byte(0x17)
            shared.byte('A'.code)
            shared.byte(0x17)
            shared.byte('C'.code)
            shared.byte(0x11)
            shared.byte(0x00)
            val indiv = LeBytes()
            indiv.byte(0x11)
            indiv.byte(0x01)
            indiv.bytes(gt)
            return framed(shared.toByteArray(), indiv.toByteArray())
        }

        fun framed(shared: ByteArray, indiv: ByteArray): ByteArray {
            val out = LeBytes()
            out.int(shared.size)
            out.int(indiv.size)
            out.bytes(shared)
            out.bytes(indiv)
            return out.toByteArray()
        }

        fun overflowChars(count: Int): ByteArray {
            val payload = ByteArray(count) { 'A'.code.toByte() }
            return byteArrayOf(0xf7.toByte()) + sizeBytes(count) + payload + byteArrayOf(MARKER)
        }

        fun sizeBytes(count: Int): ByteArray = when {
            count <= 127 -> byteArrayOf(0x11, count.toByte())
            count <= 32767 -> byteArrayOf(0x12, (count and 0xff).toByte(), ((count ushr 8) and 0xff).toByte())
            else -> byteArrayOf(
                0x13,
                (count and 0xff).toByte(),
                ((count ushr 8) and 0xff).toByte(),
                ((count ushr 16) and 0xff).toByte(),
                ((count ushr 24) and 0xff).toByte(),
            )
        }
    }
}

/** Little-endian bytes for a record built in the test. Not a BCF writer. */
private class LeBytes {
    private val out = ArrayList<Byte>()

    fun byte(value: Int) {
        out.add(value.toByte())
    }

    fun int(value: Int) {
        byte(value)
        byte(value ushr 8)
        byte(value ushr 16)
        byte(value ushr 24)
    }

    fun bytes(value: ByteArray) {
        for (item in value) out.add(item)
    }

    fun toByteArray(): ByteArray = ByteArray(out.size) { out[it] }
}
