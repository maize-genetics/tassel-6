package net.maizegenetics.tassel.genotype

/**
 * Haploid allele codes. An allele occupies the low four bits of a byte; nucleotide data uses codes 0 to 5,
 * plus [RARE] and [UNKNOWN].
 *
 * @see net.maizegenetics.dna.snp.NucleotideAlignmentConstants
 * @see net.maizegenetics.dna.snp.GenotypeTable
 */
object Allele {
    const val A: Byte = 0x0
    const val C: Byte = 0x1
    const val G: Byte = 0x2
    const val T: Byte = 0x3

    /** Insertion, written `+`. */
    const val INSERT: Byte = 0x4

    /** Gap (deletion), written `-`. */
    const val GAP: Byte = 0x5

    /** A code with no nucleotide meaning, written `X`. */
    const val UNDEFINED: Byte = 0x6

    /** Alleles too infrequent to be kept as their own state, written `Z`. */
    const val RARE: Byte = 0xE

    /** Missing allele, written `N`. */
    const val UNKNOWN: Byte = 0xF

    /** Number of nucleotide states counted by allele statistics: A, C, G, T, insert, and gap. */
    const val NUCLEOTIDE_STATES: Int = 6
}

/**
 * Diploid genotype bytes. The first allele is in the high four bits and the second in the low four bits, so `AC` is
 * `0x01` and a fully missing genotype is [UNKNOWN] (`0xFF`).
 *
 * Analyses keep genotypes as raw [Byte] values; these functions name the bit operations.
 *
 * @see net.maizegenetics.dna.snp.GenotypeTableUtils
 */
object Genotype {
    /** Both alleles missing (`N`). */
    const val UNKNOWN: Byte = 0xFF.toByte()

    /** Both alleles rare (`Z`). */
    const val RARE: Byte = 0xEE.toByte()

    /** Packs two allele codes into a genotype, keeping their order. */
    fun of(first: Byte, second: Byte): Byte = ((first.toInt() shl 4) or (second.toInt() and 0xF)).toByte()

    /** The allele in the high four bits. */
    fun firstAllele(genotype: Byte): Byte = ((genotype.toInt() shr 4) and 0xF).toByte()

    /** The allele in the low four bits. */
    fun secondAllele(genotype: Byte): Byte = (genotype.toInt() and 0xF).toByte()

    /**
     * True when the two allele codes differ. As in TASSEL 5, a genotype with one missing allele (such as `A/N`)
     * counts as heterozygous here; [AlleleStats.heterozygousCount] does not count it.
     */
    fun isHeterozygous(genotype: Byte): Boolean = firstAllele(genotype) != secondAllele(genotype)

    /** True when both allele codes are equal and the genotype is not [UNKNOWN]. */
    fun isHomozygous(genotype: Byte): Boolean = genotype != UNKNOWN && !isHeterozygous(genotype)
}
