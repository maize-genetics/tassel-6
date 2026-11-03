package net.maizegenetics.tassel.genotype

import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList

/**
 * A [GenotypeTable] held in memory as one byte per genotype, stored site by site (all taxa of site 0, then all taxa
 * of site 1, and so on). Site blocks over all taxa are views of that array, so reading them copies nothing.
 *
 * @see net.maizegenetics.dna.snp.genotypecall.ByteGenotypeCallTable
 * @see net.maizegenetics.util.SuperByteMatrix
 */
class InMemoryGenotypeTable private constructor(
    override val taxa: TaxaList,
    override val positions: PositionList,
    private val genotypes: ByteArray,
) : GenotypeTable {
    init {
        val expected = taxa.size.toLong() * positions.size
        require(genotypes.size.toLong() == expected) {
            "Expected $expected genotypes (${taxa.size} taxa x ${positions.size} sites), got ${genotypes.size}"
        }
    }

    override fun genotype(taxon: Int, site: Int): Byte {
        checkTaxon(taxon)
        checkSite(site)
        return genotypes[site * numberOfTaxa + taxon]
    }

    override fun siteGenotypes(site: Int): ByteArray {
        checkSite(site)
        return genotypes.copyOfRange(site * numberOfTaxa, (site + 1) * numberOfTaxa)
    }

    override fun taxonGenotypes(taxon: Int): ByteArray {
        checkTaxon(taxon)
        return ByteArray(numberOfSites) { genotypes[it * numberOfTaxa + taxon] }
    }

    override fun siteBlocks(blockSize: Int, taxa: IntArray?): Sequence<SiteBlock> {
        if (taxa != null) return super.siteBlocks(blockSize, taxa)
        require(blockSize > 0) { "blockSize must be positive: $blockSize" }
        return (0 until numberOfSites step blockSize).asSequence().map { start ->
            val sites = start until minOf(start + blockSize, numberOfSites)
            ByteSiteBlock(sites, numberOfTaxa, genotypes, offset = start * numberOfTaxa)
        }
    }

    override fun siteStats(site: Int): AlleleStats {
        checkSite(site)
        return AlleleStats.of(genotypes, site * numberOfTaxa, (site + 1) * numberOfTaxa)
    }

    override fun toString(): String = "InMemoryGenotypeTable(${taxa.size} taxa, ${positions.size} sites)"

    private fun checkSite(site: Int) {
        if (site < 0 || site >= numberOfSites) throw IndexOutOfBoundsException("Site $site of $numberOfSites")
    }

    private fun checkTaxon(taxon: Int) {
        if (taxon < 0 || taxon >= numberOfTaxa) throw IndexOutOfBoundsException("Taxon $taxon of $numberOfTaxa")
    }

    companion object {
        /**
         * A table over a copy of [genotypes], stored site by site: the genotype of taxon `t` at site `s` is at
         * `s * taxa.size + t`.
         */
        fun of(taxa: TaxaList, positions: PositionList, genotypes: ByteArray): InMemoryGenotypeTable =
            InMemoryGenotypeTable(taxa, positions, genotypes.copyOf())

        /** A table whose genotype of taxon `t` at site `s` is `genotype(t, s)`. */
        fun of(
            taxa: TaxaList,
            positions: PositionList,
            genotype: (taxon: Int, site: Int) -> Byte,
        ): InMemoryGenotypeTable {
            val genotypes = ByteArray(taxa.size * positions.size)
            for (site in positions.indices) {
                for (taxon in taxa.indices) genotypes[site * taxa.size + taxon] = genotype(taxon, site)
            }
            return InMemoryGenotypeTable(taxa, positions, genotypes)
        }

        /** A table that keeps [genotypes] without copying, for readers in this module that built the array. */
        internal fun wrap(taxa: TaxaList, positions: PositionList, genotypes: ByteArray): InMemoryGenotypeTable =
            InMemoryGenotypeTable(taxa, positions, genotypes)
    }
}
