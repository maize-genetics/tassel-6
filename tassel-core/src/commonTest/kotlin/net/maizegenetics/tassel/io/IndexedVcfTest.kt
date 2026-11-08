package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.io.deflate.bgzfBlocks
import net.maizegenetics.tassel.io.index.BgzfChunk
import net.maizegenetics.tassel.io.index.readGenomicIndex
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * Indexed region reads against the htslib 1.25 oracle (ADR-0021).
 *
 * Equality is [Vcf.read] of the fixture header plus the `.record.txt` body. `chr6D` is declared in
 * the wheat header and then left out of the CSI, so a query for it is an unknown contig. The empty
 * known range is `chr3B:800000000-851934019`.
 */
class IndexedVcfTest {
    @Test
    fun regionSyntaxFollowsBcftools() {
        assertEquals(VcfRegion("chr3B", null, null), parseVcfRegion("chr3B"))
        assertEquals(VcfRegion("chr3B", 1000, 2000), parseVcfRegion("chr3B:1000-2000"))
        assertEquals(VcfRegion("chr3B", 1000, null), parseVcfRegion("chr3B:1000-"))
        assertEquals(VcfRegion("chr3B", 600_000_000, 600_000_000), parseVcfRegion("chr3B:600000000"))
        assertEquals(
            VcfRegion("chr3B", 599_999_999, 600_000_001),
            parseVcfRegion("chr3B:599,999,999-600,000,001"),
        )
        assertEquals(VcfRegion("chr:3B", 1000, 2000), parseVcfRegion("{chr:3B}:1000-2000"))
        assertEquals(VcfRegion("chr:3B", null, null), parseVcfRegion("{chr:3B}"))

        assertContains(failure { parseVcfRegion("chr3B:2000-1000") }, "inverted")
        assertContains(failure { parseVcfRegion("1:10068726-157104") }, "inverted")
        assertContains(failure { parseVcfRegion("chr3B:5000-1000") }, "inverted")
        assertContains(failure { parseVcfRegion("") }, "malformed")
        assertContains(failure { parseVcfRegion("chr3B:") }, "malformed")
        assertContains(failure { parseVcfRegion(":1000") }, "malformed")
        assertContains(failure { parseVcfRegion("chr3B:1000-abc") }, "malformed")
        assertContains(failure { parseVcfRegion("{chr3B") }, "malformed")
        assertContains(failure { parseVcfRegion("chr3B:1-2,chr3B:3-4") }, "comma")
        assertContains(failure { parseVcfRegion("chr3B:" + "9".repeat(40)) }, "does not fit in a Long")
    }

    @Test
    fun contigResolutionMatchesTheIndexThenChromosomeOf() {
        assertEquals("chr3B", resolveContigName("chr3B", WHEAT_NAMES))
        assertEquals("3B", resolveContigName("3B", listOf("chr3B", "3B")))
        val ambiguous = failure { resolveContigName("3b", listOf("chr3B", "3B")) }
        assertContains(ambiguous, "ambiguous")
        assertContains(ambiguous, "chr3B")
        assertContains(ambiguous, "3B")
        val unknown = failure { resolveContigName("chrUn", WHEAT_NAMES) }
        assertContains(unknown, "chrUn")
        assertContains(unknown, "chr3B")
        assertContains(unknown, "chr1A")
    }

    @Test
    fun spanFollowsTbxParse1() {
        val sequence = vcfReferenceSpan("chr3B\t536870000\t.\tACGT\tA\t.\tPASS\t.\tGT\t0/1")
        assertEquals(536869999L, sequence.beg)
        assertEquals(536870003L, sequence.end)
        assertTrue(sequence.overlaps(parseVcfRegion("chr3B:536870002-536870010")))

        val infoEnd = vcfReferenceSpan("chr3B\t163840\t.\tN\t<DEL>\t.\tPASS\tEND=229376;SVTYPE=DEL\tGT\t0/1")
        assertEquals(163839L, infoEnd.beg)
        assertEquals(229376L, infoEnd.end)
        assertTrue(infoEnd.overlaps(parseVcfRegion("chr3B:180100-180200")))

        val svlen = vcfReferenceSpan("chr3B\t600000000\t.\tN\t<DEL>\t.\tPASS\tSVLEN=-40000;SVTYPE=DEL\tGT\t0/1")
        assertEquals(599999999L, svlen.beg)
        assertEquals(600039999L, svlen.end)
        assertTrue(svlen.overlaps(parseVcfRegion("chr3B:600010000-600020000")))

        val gvcf = vcfReferenceSpan("chr1\t10\t.\tN\t<*>\t.\t.\t.\tGT:LEN\t0/1:25")
        assertEquals(9L, gvcf.beg)
        assertEquals(34L, gvcf.end)
    }

    @Test
    fun tutorialRegionsMatchTheRecordOracleThroughCsiAndTbi() = runTest {
        val bytes = TestData.readBytes(TUTORIAL)
        val header = headerText(bytes)
        openIndexed(bytes, TestData.readBytes("$TUTORIAL.csi")).use { csi ->
            openIndexed(bytes, TestData.readBytes("$TUTORIAL.tbi")).use { tbi ->
                assertEquals(TestData.readText(TUTORIAL_NAMES).trimEnd().lines(), csi.contigs.map { it.name })
                assertEquals(csi.contigs.map { it.name }, tbi.contigs.map { it.name })
                assertTrue(csi.contigs.all { it.length == null })
                for (stem in TUTORIAL_EQUAL) {
                    val expected = oracle(header, "htslib/tutorial/regions/$stem.record.txt")
                    val region = parseVcfRegion(regionText(stem))
                    assertSameTable(expected, csi.read(region).table, stem)
                    assertSameTable(expected, tbi.read(region).table, stem)
                }
            }
        }
    }

    @Test
    fun wheatRegionsMatchTheRecordOracle() = runTest {
        val bytes = TestData.readBytes(WHEAT)
        val header = headerText(bytes)
        openIndexed(bytes, TestData.readBytes("$WHEAT.csi")).use { indexed ->
            assertEquals(WHEAT_NAMES, indexed.contigs.map { it.name })
            assertEquals(851_934_019L, indexed.contigs.first { it.name == "chr3B" }.length)
            assertNull(indexed.contigs.find { it.name == "chr6D" })
            for (stem in WHEAT_EQUAL) {
                val expected = oracle(header, "htslib/wheat/regions/$stem.record.txt")
                assertSameTable(expected, indexed.read(parseVcfRegion(regionText(stem))).table, stem)
            }
        }
    }

    @Test
    fun wholeContigReadsMatchTheWholeFile() = runTest {
        assertWholeFile(TUTORIAL, "$TUTORIAL.csi")
        assertWholeFile(WHEAT, "$WHEAT.csi")
    }

    @Test
    fun aSinglePositionReadsOnlyTheHeaderIndexAndChunkBlocks() = runTest {
        val file = TestData.readBytes(WHEAT)
        val indexFile = TestData.readBytes("$WHEAT.csi")
        val data = CountingRandomAccessSource(ByteArrayRandomAccessSource(file))
        val index = CountingRandomAccessSource(ByteArrayRandomAccessSource(indexFile))
        val indexed = IndexedVcf.open(data, index)
        val headerBytes = data.bytes
        val headerRanges = data.ranges.size
        indexed.use { it.read(parseVcfRegion("chr3B:1000")) }
        val queryRanges = data.ranges.drop(headerRanges)
        val members = membersOf(file)
        val genomic = readGenomicIndex(Buffer().also { it.write(indexFile) })
        val chunks = genomic.chunks(genomic.contigIndex("chr3B"), 999, 1000)
        val allowed = allowedMembers(members, chunks)
        for (range in queryRanges) {
            var offset = range.first
            while (offset <= range.last) {
                val member = members.indexOfFirst { offset in it }
                assertTrue(member in allowed, "offset $offset fell in BGZF member $member; allowed $allowed")
                offset = members[member].last + 1
            }
        }
        assertEquals(indexFile.size.toLong(), index.bytes)
        assertTrue(data.bytes < file.size, "chr3B:1000 read ${data.bytes} of ${file.size} bytes")
        println(
            "chr3B:1000 bytes header=$headerBytes query=${data.bytes - headerBytes} " +
                "index=${index.bytes} data=${data.bytes} file=${file.size} chunks=${chunks.size}",
        )
    }

    @Test
    fun noMatchDecisionsFollowTheAdr() = runTest {
        val tutorial = TestData.readBytes(TUTORIAL)
        openIndexed(tutorial, TestData.readBytes("$TUTORIAL.tbi")).use { tbi ->
            val past = tbi.read(parseVcfRegion("1:536870913-536871000")).table
            assertEquals(0, past.numberOfSites)
            assertEquals(281, past.numberOfTaxa)
            val unknown = assertFailsWith<IllegalArgumentException> {
                tbi.read(parseVcfRegion("99:1-1000"))
            }
            assertContains(unknown.message ?: "", "99")
            assertContains(unknown.message ?: "", "1")
        }
        openIndexed(TestData.readBytes(WHEAT), TestData.readBytes("$WHEAT.csi")).use { wheat ->
            val chrUn = assertFailsWith<IllegalArgumentException> {
                wheat.read(parseVcfRegion("chrUn:1-1000"))
            }
            assertContains(chrUn.message ?: "", "chrUn")
            assertContains(chrUn.message ?: "", "chr3B")
            val chr6d = assertFailsWith<IllegalArgumentException> {
                wheat.read(parseVcfRegion("chr6D:1-495380293"))
            }
            assertContains(chr6d.message ?: "", "chr6D")
            assertContains(chr6d.message ?: "", "chr1A")
        }
    }

    @Test
    fun aPositionPastTheIntLimitFailsAtTheTableBoundary() = runTest {
        openIndexed(TestData.readBytes(LARGE), TestData.readBytes("$LARGE.csi")).use { indexed ->
            val error = assertFailsWith<IllegalArgumentException> {
                indexed.read(parseVcfRegion("chr1:3000000000-3000000000"))
            }
            val message = error.message ?: ""
            assertContains(message, "2147483647")
            assertContains(message, "position limit")
        }
    }

    private suspend fun assertWholeFile(dataPath: String, indexPath: String) {
        val bytes = TestData.readBytes(dataPath)
        val expected = Vcf.read(ByteArrayRandomAccessSource(bytes)).table
        openIndexed(bytes, TestData.readBytes(indexPath)).use { indexed ->
            assertSameTable(expected, concatenated(indexed), dataPath)
        }
    }

    private suspend fun concatenated(indexed: IndexedVcf): GenotypeTable {
        val positions = ArrayList<net.maizegenetics.tassel.position.Position>()
        val genotypes = ArrayList<ByteArray>()
        var taxa = indexed.read(VcfRegion(indexed.contigs.first().name, null, null)).let { read ->
            collect(read.table, positions, genotypes)
            read.table.taxa
        }
        for (contig in indexed.contigs.drop(1)) {
            collect(indexed.read(VcfRegion(contig.name, null, null)).table, positions, genotypes)
        }
        val order = positions.indices.sortedBy { positions[it] }
        val width = taxa.size
        val packed = ByteArray(order.size * width)
        order.forEachIndexed { site, from -> genotypes[from].copyInto(packed, site * width) }
        return InMemoryGenotypeTable.of(taxa, PositionList(order.map { positions[it] }), packed)
    }

    private fun collect(
        table: GenotypeTable,
        positions: MutableList<net.maizegenetics.tassel.position.Position>,
        genotypes: MutableList<ByteArray>,
    ) {
        for (site in 0 until table.numberOfSites) {
            positions.add(table.positions[site])
            genotypes.add(table.siteGenotypes(site))
        }
    }

    private suspend fun oracle(header: String, recordPath: String): GenotypeTable {
        val body = TestData.readText(recordPath).trimEnd()
        val text = if (body.isEmpty()) header else header + body + "\n"
        return Vcf.read(Buffer().also { it.write(text.encodeToByteArray()) })
    }

    private fun headerText(gz: ByteArray): String {
        val plain = Buffer().also { it.write(gz) }.decompressing().source.readByteArray().decodeToString()
        val header = StringBuilder()
        var index = 0
        while (index < plain.length) {
            val newline = plain.indexOf('\n', index).let { if (it < 0) plain.length else it }
            val line = plain.substring(index, newline).removeSuffix("\r")
            if (!line.startsWith("#")) break
            header.append(line).append('\n')
            if (newline == plain.length) break
            index = newline + 1
        }
        return header.toString()
    }

    private suspend fun openIndexed(data: ByteArray, index: ByteArray): IndexedVcf =
        IndexedVcf.open(ByteArrayRandomAccessSource(data), ByteArrayRandomAccessSource(index))

    private fun regionText(stem: String): String {
        val parts = stem.split('_')
        val beg = parts.getOrNull(parts.size - 2)?.toLongOrNull()
        val end = parts.last().toLongOrNull()
        if (parts.size >= 3 && beg != null && end != null) {
            return parts.dropLast(2).joinToString("_") + ":$beg-$end"
        }
        return stem
    }

    private fun membersOf(file: ByteArray): List<LongRange> {
        var start = 0L
        return bgzfBlocks(file).map { block ->
            val end = block.to + 8L
            val range = start until end
            start = end
            range
        }
    }

    private fun allowedMembers(members: List<LongRange>, chunks: List<BgzfChunk>): Set<Int> {
        val allowed = HashSet<Int>()
        for (chunk in chunks) {
            val startCoff = chunk.start ushr 16
            val endCoff = chunk.end ushr 16
            val first = members.indexOfFirst { startCoff in it }
            val lastInside = if ((chunk.end and 0xFFFFL) == 0L) {
                members.indexOfLast { it.last + 1 <= endCoff }
            } else {
                members.indexOfFirst { endCoff in it }
            }
            assertTrue(first >= 0 && lastInside >= first, "chunk $chunk is not in the file")
            val spill = (lastInside + 1).coerceAtMost(members.lastIndex)
            for (index in first..spill) allowed.add(index)
        }
        return allowed
    }

    private fun assertSameTable(expected: GenotypeTable, actual: GenotypeTable, label: String) {
        assertEquals(expected.taxa.names, actual.taxa.names, label)
        assertEquals(expected.numberOfSites, actual.numberOfSites, label)
        for (site in 0 until expected.numberOfSites) {
            val want = expected.positions[site]
            val got = actual.positions[site]
            assertEquals(want, got, "$label site $site")
            assertEquals(want.knownVariants, got.knownVariants, "$label site $site variants")
            assertContentEquals(expected.siteGenotypes(site), actual.siteGenotypes(site), "$label site $site")
        }
    }

    private fun failure(block: () -> Unit): String =
        assertFailsWith<IllegalArgumentException>(block = block).message.orEmpty()

    private class CountingRandomAccessSource(private val inner: RandomAccessSource) : RandomAccessSource {
        var bytes: Long = 0
        val ranges = ArrayList<LongRange>()

        override val size: Long get() = inner.size

        override suspend fun readAt(position: Long, destination: ByteArray, destinationOffset: Int, length: Int): Int {
            val count = inner.readAt(position, destination, destinationOffset, length)
            if (count > 0) {
                bytes += count
                ranges.add(position until position + count)
            }
            return count
        }

        override fun close() = inner.close()
    }

    private companion object {
        const val TUTORIAL = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
        const val TUTORIAL_NAMES = "htslib/tutorial/contigs.txt"
        const val WHEAT = "htslib/wheat/wheat.vcf.gz"
        const val LARGE = "htslib/large-pos/large_pos.vcf.gz"
        val WHEAT_NAMES = listOf("chr1A", "chr2A", "chr3B", "chr7B", "chr7D")
        val TUTORIAL_EQUAL = listOf(
            "1_157104_157104",
            "1_157104_10068726",
            "1_299170078_299999999",
        )
        val WHEAT_EQUAL = listOf(
            "chr3B_1000_1000",
            "chr3B_1000_5000000",
            "chr3B",
            "chr3B_599000000_601000000",
            "chr3B_700000000_700000000",
            "chr3B_536870002_536870010",
            "chr3B_180100_180200",
            "chr3B_600010000_600020000",
            "chr3B_800000000_851934019",
        )
    }
}
