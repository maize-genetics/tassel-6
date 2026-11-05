package net.maizegenetics.tassel.genotype

import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith

class GenotypeTableTest {
    private val taxa = TaxaList.of("t0", "t1", "t2")

    private val positions = PositionList((1..5).map { Position(Chromosome.of(1), it * 100) })

    /** Site by site, three taxa per site. */
    private val codes = listOf("AAC", "RRN", "GTG", "N-+", "CCC")

    private val genotypes = codes.joinToString("").let { text ->
        ByteArray(text.length) { Nucleotide.genotype(text[it]) }
    }

    private val table = InMemoryGenotypeTable.of(taxa, positions, genotypes)

    /** Implements only [GenotypeTable.genotype], so every other member uses the interface defaults. */
    private val minimal = object : GenotypeTable {
        override val taxa = this@GenotypeTableTest.taxa
        override val positions = this@GenotypeTableTest.positions

        override fun genotype(taxon: Int, site: Int): Byte = table.genotype(taxon, site)
    }

    @Test
    fun randomAccess() {
        assertEquals(3, table.numberOfTaxa)
        assertEquals(5, table.numberOfSites)
        assertEquals(Nucleotide.genotype('T'), table.genotype(taxon = 1, site = 2))
        assertContentEquals(genotypes.copyOfRange(3, 6), table.siteGenotypes(1))
        assertEquals("ARGNC", table.taxonGenotypes(0).map(Nucleotide::iupac).joinToString(""))
        assertEquals("CNG+C", table.taxonGenotypes(2).map(Nucleotide::iupac).joinToString(""))
        assertFailsWith<IndexOutOfBoundsException> { table.genotype(3, 0) }
        assertFailsWith<IndexOutOfBoundsException> { table.genotype(0, 5) }
        assertFailsWith<IndexOutOfBoundsException> { table.siteGenotypes(-1) }
    }

    @Test
    fun copiesItsInput() {
        val input = genotypes.copyOf()
        val copy = InMemoryGenotypeTable.of(taxa, positions, input)
        input.fill(Genotype.UNKNOWN)
        assertEquals(Nucleotide.genotype('A'), copy.genotype(0, 0))
        copy.siteGenotypes(0).fill(Genotype.UNKNOWN)
        assertEquals(Nucleotide.genotype('A'), copy.genotype(0, 0))
    }

    @Test
    fun rejectsWrongGenotypeCount() {
        assertFailsWith<IllegalArgumentException> { InMemoryGenotypeTable.of(taxa, positions, ByteArray(14)) }
    }

    @Test
    fun buildsFromAFunction() {
        val built = InMemoryGenotypeTable.of(taxa, positions) { taxon, site -> table.genotype(taxon, site) }
        for (site in 0 until 5) assertContentEquals(table.siteGenotypes(site), built.siteGenotypes(site))
    }

    @Test
    fun siteBlocksCoverEverySiteInOrder() {
        for (source in listOf(table, minimal)) {
            for (blockSize in listOf(1, 2, 3, 5, 4096)) {
                val blocks = source.siteBlocks(blockSize).toList()
                assertEquals((0 until 5).toList(), blocks.flatMap { it.sites })
                assertEquals((5 + blockSize - 1) / blockSize, blocks.size)
                for (block in blocks) {
                    assertEquals(3, block.numberOfTaxa)
                    for (site in block.sites) {
                        assertContentEquals(table.siteGenotypes(site), block.siteGenotypes(site))
                        for (taxon in 0 until 3) assertEquals(table.genotype(taxon, site), block.genotype(taxon, site))
                    }
                }
            }
        }
    }

    @Test
    fun siteBlocksForATaxaSubset() {
        for (source in listOf(table, minimal)) {
            val blocks = source.siteBlocks(blockSize = 2, taxa = intArrayOf(2, 0)).toList()
            assertEquals(listOf(0..1, 2..3, 4..4), blocks.map { it.sites })
            val columns = blocks.flatMap { block -> block.sites.map { block.siteGenotypes(it) } }
            assertEquals(
                listOf("CA", "NR", "GG", "+N", "CC"),
                columns.map {
                    it.map(Nucleotide::iupac).joinToString("")
                },
            )
        }
    }

    @Test
    fun siteBlockChecksIndices() {
        val block = table.siteBlocks(blockSize = 2).drop(1).first()
        assertEquals(2..3, block.sites)
        assertFailsWith<IndexOutOfBoundsException> { block.genotype(0, 1) }
        assertFailsWith<IndexOutOfBoundsException> { block.genotype(3, 2) }
        val destination = ByteArray(5) { Genotype.RARE }
        block.copySiteGenotypes(3, destination, destinationOffset = 2)
        assertEquals("ZZN-+", destination.map(Nucleotide::iupac).joinToString(""))
        assertFailsWith<IllegalArgumentException> { table.siteBlocks(0) }
        assertFailsWith<IllegalArgumentException> { table.siteBlocks(2, taxa = intArrayOf(3)) }
    }

    @Test
    fun statsAgreeAcrossAccessPaths() {
        for (site in 0 until 5) {
            val expected = AlleleStats.of(table.siteGenotypes(site))
            val block = table.siteBlocks(blockSize = 2).first { site in it.sites }
            for (stats in listOf(table.siteStats(site), minimal.siteStats(site), block.siteStats(site))) {
                assertEquals(expected.alleles, stats.alleles)
                assertEquals(expected.alleleCounts, stats.alleleCounts)
                assertEquals(expected.heterozygousCount, stats.heterozygousCount)
                assertEquals(expected.unknownCount, stats.unknownCount)
            }
        }
        val taxon = table.taxonStats(0)
        // A and G both have three copies; ties go to the lower allele code.
        assertEquals(listOf(Allele.A, Allele.G, Allele.C), taxon.alleles)
        assertEquals(listOf(3, 3, 2), taxon.alleleCounts)
        assertEquals(1, taxon.unknownCount)
        assertEquals(1, taxon.heterozygousCount)
    }

    @Test
    fun emptyTables() {
        val noSites = InMemoryGenotypeTable.of(taxa, PositionList(emptyList()), ByteArray(0))
        assertEquals(0, noSites.siteBlocks().count())
        assertEquals(0, noSites.taxonGenotypes(1).size)
        val noTaxa = InMemoryGenotypeTable.of(TaxaList.of(), positions, ByteArray(0))
        assertEquals(0, noTaxa.siteGenotypes(4).size)
        assertEquals(listOf(0..4), noTaxa.siteBlocks().map { it.sites }.toList())
    }
}
