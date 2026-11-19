package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readByteArray
import kotlinx.io.readString
import kotlinx.io.writeString
import net.maizegenetics.tassel.filter.TaxaFilter
import net.maizegenetics.tassel.filter.filterTaxa
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.io.deflate.BGZF_EOF_BLOCK
import net.maizegenetics.tassel.io.deflate.BGZF_HEADER_LENGTH
import net.maizegenetics.tassel.io.deflate.GZIP_TRAILER
import net.maizegenetics.tassel.io.deflate.bgzfBlocks
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import net.maizegenetics.tassel.testing.SuspendLazy
import net.maizegenetics.tassel.testing.TestData
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * BCF export against htslib 1.25's encoding of the TASSEL 5 VCF-export goldens.
 *
 * Inflated bytes match the committed `.bcf`. On disk, each stored data member is that uncompressed
 * size plus the BGZF framing, and the EOF member is the fixed 28 bytes. `vcf_export_edges.blocks-l0.txt`
 * is that list. The other `blocks.txt` files are default `bcftools view -Ob` sizes, which a stored
 * writer does not reproduce.
 */
class BcfWriteTest {
    @Test
    fun writesTheTutorialExport() = runTest {
        assertMatchesHtslib("mdp_genotype", tutorial())
    }

    @Test
    fun writesTheIncludedTaxaExport() = runTest {
        val included = filterTaxa(
            tutorial(),
            TaxaFilter(taxaNames = listOf("W22", "B73", "CML247", "33-16"), includeTaxa = true),
        )
        assertMatchesHtslib("include-list", included)
    }

    @Test
    fun writesTheEdgeExport() = runTest {
        val table = HapMap.read(Buffer().also { it.write(TestData.readBytes(EDGES)) })
        assertMatchesHtslib("vcf_export_edges", table)
    }

    @Test
    fun readingBackTheTutorialMatchesVcf() = runTest {
        val table = tutorial()
        assertSameTable(vcfRoundTrip(table), bcfRoundTrip(table))
    }

    @Test
    fun readingBackTheEdgesTableMatchesVcf() = runTest {
        val table = HapMap.read(Buffer().also { it.write(TestData.readBytes(EDGES)) })
        assertSameTable(vcfRoundTrip(table), bcfRoundTrip(table))
    }

    @Test
    fun reportsProgressAtTheEndOfTheOnlyBlock() = runTest {
        // 3,093 sites fit in one block, so a block-boundary report fires once, after the write.
        val reports = mutableListOf<Progress>()
        Bcf.write(tutorial(), Buffer()) { reports += it }
        assertEquals(listOf(Progress(3093, 3093, "Writing BCF")), reports)
    }

    @Test
    fun writesASiteNameLongerThanFourteenCharacters() = runTest {
        // A typed char count of 15 is the overflow form: the next integer is the real length.
        val table = InMemoryGenotypeTable.of(
            TaxaList.of("t1"),
            PositionList(listOf(Position(Chromosome.of("1"), 10, snpId = "0123456789abcde"))),
            byteArrayOf(Nucleotide.genotype("AA")),
        )
        assertSameTable(vcfRoundTrip(table), bcfRoundTrip(table))
    }

    private suspend fun assertMatchesHtslib(name: String, table: GenotypeTable) {
        val written = Buffer().also { Bcf.write(table, it) }.readByteArray()
        val oracle = TestData.readBytes("htslib/bcf/golden/$name.bcf")
        assertSameBytes(name, inflate(oracle), inflate(written))

        val writtenBlocks = bgzfBlocks(written)
        assertEquals(
            bgzfBlocks(oracle).map { it.uncompressedSize },
            writtenBlocks.map { it.uncompressedSize },
            "$name uncompressed block sizes",
        )
        assertEquals(0, writtenBlocks.last().uncompressedSize, "$name EOF")
        assertTrue(writtenBlocks.dropLast(1).all { it.uncompressedSize > 0 }, "$name has an empty data block")
        val members = writtenBlocks.mapIndexed { index, block ->
            if (index == writtenBlocks.lastIndex) {
                BGZF_EOF_BLOCK.size
            } else {
                block.uncompressedSize + STORED_MEMBER_OVERHEAD
            }
        }
        assertEquals(members, memberSizes(written), "$name stored member sizes")
        if (name == "vcf_export_edges") {
            val expected = TestData.readText("htslib/bcf/golden/$name.blocks-l0.txt")
                .lineSequence()
                .filter { it.isNotEmpty() }
                .map { it.toInt() }
                .toList()
            assertEquals(expected, members, "$name level-0 block list")
        }
    }

    private suspend fun bcfRoundTrip(table: GenotypeTable): GenotypeTable {
        val written = Buffer().also { Bcf.write(table, it) }.readByteArray()
        return Bcf.read(Buffer().also { it.write(written) }).table
    }

    private suspend fun vcfRoundTrip(table: GenotypeTable): GenotypeTable {
        val written = Buffer().also { Vcf.write(table, it) }.readString()
        return Vcf.read(Buffer().also { it.writeString(written) })
    }

    private suspend fun assertSameTable(expected: GenotypeTable, actual: GenotypeTable) {
        assertEquals(expected.taxa.names, actual.taxa.names)
        assertEquals(expected.numberOfSites, actual.numberOfSites)
        for (site in 0 until expected.numberOfSites) {
            assertEquals(expected.positions[site], actual.positions[site], "site $site")
            assertContentEquals(expected.siteGenotypes(site), actual.siteGenotypes(site), "site $site")
        }
    }

    private fun inflate(bytes: ByteArray): ByteArray =
        Buffer().also { it.write(bytes) }.decompressing().source.readByteArray()

    private fun assertSameBytes(name: String, expected: ByteArray, actual: ByteArray) {
        val shared = minOf(expected.size, actual.size)
        var index = 0
        while (index < shared && expected[index] == actual[index]) index++
        if (index == shared && expected.size == actual.size) return
        throw AssertionError(
            "$name inflated bytes differ at offset $index " +
                "(expected ${expected.size} bytes, actual ${actual.size}): " +
                "expected ${hex(expected, index)} actual ${hex(actual, index)}",
        )
    }

    private fun hex(bytes: ByteArray, at: Int): String {
        if (at >= bytes.size) return "end"
        val end = minOf(bytes.size, at + 8)
        return (at until end).joinToString(" ") { (bytes[it].toInt() and 0xFF).toString(16).padStart(2, '0') }
    }

    /** `BSIZE + 1` of each BGZF member, which is the on-disk size. */
    private fun memberSizes(file: ByteArray): List<Int> {
        val sizes = ArrayList<Int>()
        var start = 0
        while (start < file.size) {
            val size = memberSize(file, start)
            sizes.add(size)
            start += size
        }
        return sizes
    }

    private fun memberSize(file: ByteArray, start: Int): Int {
        require(start + BGZF_HEADER_LENGTH <= file.size) { "truncated gzip header at $start" }
        val flags = file[start + 3].toInt() and 0xFF
        require(flags and 4 != 0) { "BGZF block at $start has no extra field" }
        var cursor = start + 10
        val extraLength = u16(file, cursor)
        cursor += 2
        val extraEnd = cursor + extraLength
        var blockSize = -1
        var sub = cursor
        while (sub + 4 <= extraEnd) {
            val id1 = file[sub].toInt() and 0xFF
            val id2 = file[sub + 1].toInt() and 0xFF
            val subLength = u16(file, sub + 2)
            if (id1 == 'B'.code && id2 == 'C'.code && subLength == 2) {
                blockSize = u16(file, sub + 4)
            }
            sub += 4 + subLength
        }
        require(blockSize >= 0) { "BGZF block at $start has no BC subfield" }
        return blockSize + 1
    }

    private fun u16(file: ByteArray, at: Int): Int =
        (file[at].toInt() and 0xFF) or ((file[at + 1].toInt() and 0xFF) shl 8)

    private companion object {
        const val TUTORIAL = "tutorial/mdp_genotype.hmp.txt"
        const val EDGES = "tutorial/vcf_export_edges.hmp.txt"

        /** Stored DEFLATE header in front of the uncompressed bytes: `BFINAL`, `LEN`, `NLEN`. */
        const val STORED_DEFLATE_HEADER = 5
        const val STORED_MEMBER_OVERHEAD = BGZF_HEADER_LENGTH + STORED_DEFLATE_HEADER + GZIP_TRAILER

        val tutorial = SuspendLazy {
            HapMap.read(Buffer().also { it.write(TestData.readBytes(TUTORIAL)) })
        }
    }
}
