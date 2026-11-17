package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
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
import kotlin.test.assertTrue

/**
 * One CSI region of a BCF against [IndexedVcf.read] of the sibling `.vcf.gz` (ADR-0022).
 *
 * `chr6D` is in the wheat BCF header and has no records, so the query is an empty table. The VCF CSI omits it, and
 * [IndexedVcf.read] reports an unknown contig instead. The tutorial's `##contig` lines have no length. Wheat keeps
 * the lengths from the VCF, including `chr3B`.
 */
class IndexedBcfTest {
    @Test
    fun regionSyntaxRejectsTheAdrCases() {
        assertContains(failure { parseVcfRegion("chr3B:2000-1000") }, "inverted")
        assertContains(failure { parseVcfRegion("1:10068726-157104") }, "inverted")
        assertContains(failure { parseVcfRegion("") }, "malformed")
        assertContains(failure { parseVcfRegion("chr3B:") }, "malformed")
        assertContains(failure { parseVcfRegion("chr3B:1000-abc") }, "malformed")
        assertContains(failure { parseVcfRegion("chr3B:1-2,chr3B:3-4") }, "comma")
        assertContains(failure { parseVcfRegion("chr3B:" + "9".repeat(40)) }, "does not fit in a Long")
    }

    @Test
    fun tutorialRegionsMatchIndexedVcf() = runTest {
        val region = TUTORIAL_EQUAL
        openGenotypes(TUTORIAL_BCF, "$TUTORIAL_BCF.csi").use { bcf ->
            openVcf(TUTORIAL_VCF, "$TUTORIAL_VCF.csi").use { vcf ->
                assertEquals((1..10).map { it.toString() }, bcf.contigs.map { it.name })
                assertTrue(bcf.contigs.all { it.length == null })
                for (stem in region) {
                    val parsed = parseVcfRegion(regionText(stem))
                    assertSameTable(vcf.read(parsed).table, bcf.read(parsed).table, stem)
                }
            }
        }
    }

    @Test
    fun wheatRegionsMatchIndexedVcf() = runTest {
        openGenotypes(WHEAT_BCF, "$WHEAT_BCF.csi").use { bcf ->
            openVcf(WHEAT_VCF, "$WHEAT_VCF.csi").use { vcf ->
                assertEquals(WHEAT_HEADER, bcf.contigs.map { it.name })
                assertTrue(bcf.contigs.all { it.length != null })
                assertEquals(851_934_019L, bcf.contigs.first { it.name == "chr3B" }.length)
                assertEquals(495_380_293L, bcf.contigs.first { it.name == "chr6D" }.length)
                for (stem in WHEAT_EQUAL) {
                    val parsed = parseVcfRegion(regionText(stem))
                    assertSameTable(vcf.read(parsed).table, bcf.read(parsed).table, stem)
                }
            }
        }
    }

    @Test
    fun wholeContigReadsMatchTheWholeFile() = runTest {
        assertWholeFile(TUTORIAL_BCF, "$TUTORIAL_BCF.csi")
        assertWholeFile(WHEAT_BCF, "$WHEAT_BCF.csi")
    }

    @Test
    fun aSinglePositionReadsOnlyTheHeaderIndexAndChunkBlocks() = runTest {
        val file = TestData.readBytes(WHEAT_BCF)
        val indexFile = TestData.readBytes("$WHEAT_BCF.csi")
        val data = CountingRandomAccessSource(ByteArrayRandomAccessSource(file))
        val index = CountingRandomAccessSource(ByteArrayRandomAccessSource(indexFile))
        val indexed = IndexedBcf.open(data, index)
        val headerBytes = data.bytes
        val headerRanges = data.ranges.size
        indexed.use { it.read(parseVcfRegion("chr3B:1000")) }
        val queryRanges = data.ranges.drop(headerRanges)
        val members = membersOf(file)
        val genomic = readGenomicIndex(Buffer().also { it.write(indexFile) }, WHEAT_HEADER)
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
        openGenotypes(TUTORIAL_BCF, "$TUTORIAL_BCF.csi").use { tutorial ->
            val past = tutorial.read(parseVcfRegion("1:4294967297")).table
            assertEquals(0, past.numberOfSites)
            assertEquals(281, past.numberOfTaxa)
            val empty = tutorial.read(parseVcfRegion("1:1-1000")).table
            assertEquals(0, empty.numberOfSites)
            assertEquals(281, empty.numberOfTaxa)
            val unknown = assertFailsWith<IllegalArgumentException> {
                tutorial.read(parseVcfRegion("99:1-1000"))
            }
            assertContains(unknown.message ?: "", "99")
            assertContains(unknown.message ?: "", "1")
            val inverted = assertFailsWith<IllegalArgumentException> {
                tutorial.read(VcfRegion("1", 10_068_726, 157_104))
            }
            assertContains(inverted.message ?: "", "inverted")
        }
        openGenotypes(WHEAT_BCF, "$WHEAT_BCF.csi").use { wheat ->
            val taxa = wheat.read(parseVcfRegion("chr3B:1000")).table.numberOfTaxa
            val chr6d = wheat.read(parseVcfRegion("chr6D:1-495380293")).table
            assertEquals(0, chr6d.numberOfSites)
            assertEquals(taxa, chr6d.numberOfTaxa)
            val emptyRange = wheat.read(parseVcfRegion("chr3B:800000000-851934019")).table
            assertEquals(0, emptyRange.numberOfSites)
            assertEquals(taxa, emptyRange.numberOfTaxa)
            val chrUn = assertFailsWith<IllegalArgumentException> {
                wheat.read(parseVcfRegion("chrUn:1-1000"))
            }
            assertContains(chrUn.message ?: "", "chrUn")
            assertContains(chrUn.message ?: "", "chr1A")
        }
    }

    @Test
    fun theIntMaxPositionReads() = runTest {
        val expected = Vcf.read(Buffer().also { it.write(TestData.readText(BOUNDARY_VIEW).encodeToByteArray()) })
        openGenotypes(BOUNDARY_BCF, "$BOUNDARY_BCF.csi").use { indexed ->
            assertSameTable(expected, indexed.read(parseVcfRegion("chr1:2147483647")).table, BOUNDARY_BCF)
        }
    }

    @Test
    fun vcfMagicStaysOnIndexedVcf() = runTest {
        val region = parseVcfRegion("1:157104-157104")
        openGenotypes(TUTORIAL_VCF, "$TUTORIAL_VCF.csi").use { shared ->
            openVcf(TUTORIAL_VCF, "$TUTORIAL_VCF.csi").use { vcf ->
                assertEquals(vcf.contigs.map { it.name }, shared.contigs.map { it.name })
                assertSameTable(vcf.read(region).table, shared.read(region).table, "vcf")
            }
        }
    }

    @Test
    fun anIndexThatAlreadyHasNamesIsRejected() = runTest {
        val bcf = TestData.readBytes(TUTORIAL_BCF)
        val tbi = assertFailsWith<IllegalArgumentException> {
            IndexedGenotypes.open(
                ByteArrayRandomAccessSource(bcf),
                ByteArrayRandomAccessSource(TestData.readBytes("$TUTORIAL_VCF.tbi")),
            )
        }
        assertEquals("TBI index already contains contig names", tbi.message)
        val csi = assertFailsWith<IllegalArgumentException> {
            IndexedGenotypes.open(
                ByteArrayRandomAccessSource(bcf),
                ByteArrayRandomAccessSource(TestData.readBytes("$TUTORIAL_VCF.csi")),
            )
        }
        assertEquals("CSI index already contains contig names", csi.message)
    }

    private suspend fun assertWholeFile(dataPath: String, indexPath: String) {
        val bytes = TestData.readBytes(dataPath)
        val expected = Bcf.read(ByteArrayRandomAccessSource(bytes)).table
        openGenotypes(bytes, TestData.readBytes(indexPath)).use { indexed ->
            assertSameTable(expected, concatenated(indexed), dataPath)
        }
    }

    private suspend fun concatenated(indexed: IndexedGenotypes): GenotypeTable {
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

    private suspend fun openGenotypes(data: ByteArray, index: ByteArray): IndexedGenotypes =
        IndexedGenotypes.open(ByteArrayRandomAccessSource(data), ByteArrayRandomAccessSource(index))

    private suspend fun openGenotypes(dataPath: String, indexPath: String): IndexedGenotypes =
        openGenotypes(TestData.readBytes(dataPath), TestData.readBytes(indexPath))

    private suspend fun openVcf(dataPath: String, indexPath: String): IndexedVcf = IndexedVcf.open(
        ByteArrayRandomAccessSource(TestData.readBytes(dataPath)),
        ByteArrayRandomAccessSource(TestData.readBytes(indexPath)),
    )

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
        const val TUTORIAL_BCF = "htslib/bcf/tutorial/mdp_genotype_w_depth.bcf"
        const val TUTORIAL_VCF = "htslib/tutorial/mdp_genotype_w_depth.vcf.gz"
        const val WHEAT_BCF = "htslib/bcf/wheat/wheat.bcf"
        const val WHEAT_VCF = "htslib/wheat/wheat.vcf.gz"
        const val BOUNDARY_BCF = "htslib/bcf/boundary/pos2147483647.bcf"
        const val BOUNDARY_VIEW = "htslib/bcf/boundary/pos2147483647.view.vcf"
        val WHEAT_HEADER = listOf(
            "chr1A", "chr1B", "chr1D",
            "chr2A", "chr2B", "chr2D",
            "chr3A", "chr3B", "chr3D",
            "chr4A", "chr4B", "chr4D",
            "chr5A", "chr5B", "chr5D",
            "chr6A", "chr6B", "chr6D",
            "chr7A", "chr7B", "chr7D",
        )
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
