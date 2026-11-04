package net.maizegenetics.tassel.genotype

import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList

/**
 * Diploid genotype calls for a set of taxa at a set of sites. This is a frozen contract (ADR-0012): every reader,
 * writer, and analysis works through it, and changing it needs an ADR.
 *
 * Genotypes are bytes in the TASSEL 5 nucleotide encoding (see [Genotype] and [Nucleotide]). Taxa and sites are
 * addressed by index into [taxa] and [positions].
 *
 * Analyses should read genotypes through [siteBlocks], which lets them accumulate results without the whole matrix in
 * memory and gives a natural unit for progress and cancellation. Random access ([genotype], [siteGenotypes],
 * [taxonGenotypes]) is available for code that needs it.
 *
 * Implementations are immutable and safe to read from several threads.
 *
 * @see net.maizegenetics.dna.snp.GenotypeTable
 * @see net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable
 */
interface GenotypeTable {
    /** The taxa, in table order. */
    val taxa: TaxaList

    /** The sites, in table order. */
    val positions: PositionList

    /** Number of taxa. */
    val numberOfTaxa: Int get() = taxa.size

    /** Number of sites. */
    val numberOfSites: Int get() = positions.size

    /** The genotype of [taxon] at [site]. */
    fun genotype(taxon: Int, site: Int): Byte

    /** Genotypes of all taxa at [site], in taxa order. The caller owns the returned array. */
    fun siteGenotypes(site: Int): ByteArray = ByteArray(numberOfTaxa) { genotype(it, site) }

    /** Genotypes of [taxon] at all sites, in site order. The caller owns the returned array. */
    fun taxonGenotypes(taxon: Int): ByteArray = ByteArray(numberOfSites) { genotype(taxon, it) }

    /**
     * The table as consecutive blocks of up to [blockSize] sites, covering every site once, in order.
     *
     * Blocks hold all taxa, or only [taxa] (table taxon indices, in the order given) if it is not null; inside a
     * block, taxon `i` is then `taxa[i]`. The sequence is lazy, so a table backed by storage outside memory can load
     * one block at a time.
     */
    fun siteBlocks(blockSize: Int = DEFAULT_SITE_BLOCK_SIZE, taxa: IntArray? = null): Sequence<SiteBlock> {
        require(blockSize > 0) { "blockSize must be positive: $blockSize" }
        val selected = taxa?.copyOf() ?: IntArray(numberOfTaxa) { it }
        selected.forEach { require(it in 0 until numberOfTaxa) { "Taxon index out of range: $it" } }
        return (0 until numberOfSites step blockSize).asSequence().map { start ->
            val sites = start until minOf(start + blockSize, numberOfSites)
            val genotypes = ByteArray(sites.count() * selected.size)
            for (site in sites) {
                val offset = (site - start) * selected.size
                for (taxon in selected.indices) genotypes[offset + taxon] = genotype(selected[taxon], site)
            }
            SiteBlock.of(sites, selected.size, genotypes)
        }
    }

    /** Allele statistics over all taxa at [site]. */
    fun siteStats(site: Int): AlleleStats = AlleleStats.of(siteGenotypes(site))

    /** Allele statistics over all sites of [taxon]. */
    fun taxonStats(taxon: Int): AlleleStats = AlleleStats.of(taxonGenotypes(taxon))

    companion object {
        /** Default number of sites per block. */
        const val DEFAULT_SITE_BLOCK_SIZE: Int = 4096
    }
}

/**
 * Genotypes of a run of consecutive [sites] for [numberOfTaxa] taxa, as produced by [GenotypeTable.siteBlocks].
 *
 * Sites are addressed by their table index (any value in [sites]); taxa by their index within the block.
 */
interface SiteBlock {
    /** Table indices of the sites in this block. */
    val sites: IntRange

    /** Number of taxa in this block. */
    val numberOfTaxa: Int

    /** The genotype of block taxon [taxon] at table site [site]. */
    fun genotype(taxon: Int, site: Int): Byte

    /** Copies the genotypes of every block taxon at [site] into [destination], starting at [destinationOffset]. */
    fun copySiteGenotypes(site: Int, destination: ByteArray, destinationOffset: Int = 0)

    /** Genotypes of every block taxon at [site]. The caller owns the returned array. */
    fun siteGenotypes(site: Int): ByteArray = ByteArray(numberOfTaxa).also { copySiteGenotypes(site, it) }

    /** Allele statistics over the block taxa at [site]. */
    fun siteStats(site: Int): AlleleStats = AlleleStats.of(siteGenotypes(site))

    companion object {
        /**
         * A block over [genotypes] stored site by site: the genotype of taxon `t` at site `s` is at
         * `(s - sites.first) * numberOfTaxa + t`. The block keeps [genotypes] without copying, so the caller must not
         * modify it afterward.
         */
        fun of(sites: IntRange, numberOfTaxa: Int, genotypes: ByteArray): SiteBlock =
            ByteSiteBlock(sites, numberOfTaxa, genotypes, offset = 0)
    }
}

/** A [SiteBlock] over a site-major byte array, starting at [offset]. */
internal class ByteSiteBlock(
    override val sites: IntRange,
    override val numberOfTaxa: Int,
    private val genotypes: ByteArray,
    private val offset: Int,
) : SiteBlock {
    init {
        require(numberOfTaxa >= 0) { "numberOfTaxa must not be negative: $numberOfTaxa" }
        require(offset + sites.count() * numberOfTaxa <= genotypes.size) {
            "Expected ${sites.count() * numberOfTaxa} genotypes from offset $offset, but the array has ${genotypes.size}"
        }
    }

    override fun genotype(taxon: Int, site: Int): Byte = genotypes[index(taxon, site)]

    override fun copySiteGenotypes(site: Int, destination: ByteArray, destinationOffset: Int) {
        val start = index(0, site)
        genotypes.copyInto(destination, destinationOffset, start, start + numberOfTaxa)
    }

    override fun siteStats(site: Int): AlleleStats {
        val start = index(0, site)
        return AlleleStats.of(genotypes, start, start + numberOfTaxa)
    }

    private fun index(taxon: Int, site: Int): Int {
        if (site !in sites) throw IndexOutOfBoundsException("Site $site is not in block $sites")
        if (taxon < 0 || taxon >= numberOfTaxa) throw IndexOutOfBoundsException("Taxon $taxon of $numberOfTaxa")
        return offset + (site - sites.first) * numberOfTaxa + taxon
    }
}
