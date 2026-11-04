package net.maizegenetics.tassel.genotype

/**
 * Allele and genotype counts for one site (all taxa) or one taxon (all sites).
 *
 * Alleles are the nucleotide states A, C, G, T, insert, and gap. Each genotype contributes its two allele codes, and
 * codes outside those six (unknown, rare, undefined) are not counted as alleles. [alleles] is sorted by count,
 * highest first, with ties broken by allele code (A before C), exactly as TASSEL 5 does.
 *
 * Every genotype falls into one of four groups, which sum to [numberOfGenotypes]: [unknownCount] (`N/N`),
 * [partlyUnknownCount] (one allele unknown, such as `A/N`), [heterozygousCount], and [homozygousCount].
 *
 * @see net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache
 * @see net.maizegenetics.dna.snp.genotypecall.Stats
 */
class AlleleStats private constructor(
    /** Alleles present, most frequent first. */
    val alleles: List<Byte>,
    /** Number of allele copies (gametes) for each entry of [alleles]. */
    val alleleCounts: List<Int>,
    /** Number of genotypes counted. */
    val numberOfGenotypes: Int,
    /** Genotypes with both alleles unknown. */
    val unknownCount: Int,
    /** Genotypes with exactly one allele unknown. */
    val partlyUnknownCount: Int,
    /**
     * Genotypes with two known, different allele codes. TASSEL 5's `GenotypeTable.heterozygousCount(site)` also
     * counts partly unknown genotypes; that value is `heterozygousCount + partlyUnknownCount`.
     */
    val heterozygousCount: Int,
    /** Genotypes with two known, equal allele codes. */
    val homozygousCount: Int,
) {
    /** Most frequent allele, or [Allele.UNKNOWN] if there are none. */
    val majorAllele: Byte get() = alleles.getOrElse(0) { Allele.UNKNOWN }

    /** Second most frequent allele, or [Allele.UNKNOWN] if there are fewer than two. */
    val minorAllele: Byte get() = alleles.getOrElse(1) { Allele.UNKNOWN }

    /** Copies of [majorAllele], or 0. */
    val majorAlleleCount: Int get() = alleleCounts.getOrElse(0) { 0 }

    /** Copies of [minorAllele], or 0. */
    val minorAlleleCount: Int get() = alleleCounts.getOrElse(1) { 0 }

    /** Copies of all alleles except the major one. */
    val allMinorAlleleCount: Int get() = alleleGameteCount - majorAlleleCount

    /** Total copies of the six counted alleles. TASSEL 5: `AlleleFreqCache.totalGametesNonMissingForSite`. */
    val alleleGameteCount: Int get() = alleleCounts.sum()

    /** [majorAlleleCount] over [alleleGameteCount], or 0.0 if there are no alleles. */
    val majorAlleleFrequency: Double
        get() = if (alleles.isEmpty()) 0.0 else majorAlleleCount.toDouble() / alleleGameteCount.toDouble()

    /** [minorAlleleCount] over [alleleGameteCount], or 0.0 if there are fewer than two alleles. */
    val minorAlleleFrequency: Double
        get() = if (alleles.size < 2) 0.0 else minorAlleleCount.toDouble() / alleleGameteCount.toDouble()

    /** Genotypes that are not `N/N`. TASSEL 5: `totalNonMissingForSite` and `totalNonMissingForTaxon`. */
    val nonMissingCount: Int get() = numberOfGenotypes - unknownCount

    /**
     * Allele codes that are not unknown, including rare and undefined codes. TASSEL 5:
     * `GenotypeTable.totalGametesNonMissingForSite` and `totalGametesNonMissingForTaxon`.
     */
    val nonMissingGameteCount: Int get() = 2 * numberOfGenotypes - 2 * unknownCount - partlyUnknownCount

    /** [heterozygousCount] over [numberOfGenotypes]; NaN when there are no genotypes, as in TASSEL 5. */
    val proportionHeterozygous: Double get() = heterozygousCount.toDouble() / numberOfGenotypes.toDouble()

    /** [nonMissingCount] over [numberOfGenotypes]; NaN when there are no genotypes, as in TASSEL 5. */
    val proportionNotMissing: Double get() = nonMissingCount.toDouble() / numberOfGenotypes.toDouble()

    /** True if an insert or gap allele is present. */
    val hasIndel: Boolean get() = alleles.any { it == Allele.INSERT || it == Allele.GAP }

    companion object {
        /** Counts the genotypes in [genotypes] from [fromIndex] (inclusive) to [toIndex] (exclusive). */
        fun of(genotypes: ByteArray, fromIndex: Int = 0, toIndex: Int = genotypes.size): AlleleStats {
            val counts = IntArray(Allele.NUCLEOTIDE_STATES)
            val groups = IntArray(GROUPS)
            for (i in fromIndex until toIndex) {
                val genotype = genotypes[i].toInt() and 0xFF
                groups[GROUP_OF[genotype].toInt()]++
                val first = genotype shr 4
                val second = genotype and 0xF
                if (first < Allele.NUCLEOTIDE_STATES) counts[first]++
                if (second < Allele.NUCLEOTIDE_STATES) counts[second]++
            }
            val present = (0 until Allele.NUCLEOTIDE_STATES)
                .filter { counts[it] > 0 }
                .sortedWith(compareByDescending<Int> { counts[it] }.thenBy { it })
            return AlleleStats(
                alleles = present.map { it.toByte() },
                alleleCounts = present.map { counts[it] },
                numberOfGenotypes = toIndex - fromIndex,
                unknownCount = groups[UNKNOWN],
                partlyUnknownCount = groups[PARTLY_UNKNOWN],
                heterozygousCount = groups[HETEROZYGOUS],
                homozygousCount = groups[HOMOZYGOUS],
            )
        }

        private const val UNKNOWN = 0
        private const val PARTLY_UNKNOWN = 1
        private const val HETEROZYGOUS = 2
        private const val HOMOZYGOUS = 3
        private const val GROUPS = 4

        private val GROUP_OF = ByteArray(256) {
            val genotype = it.toByte()
            val first = Genotype.firstAllele(genotype)
            val second = Genotype.secondAllele(genotype)
            when {
                genotype == Genotype.UNKNOWN -> UNKNOWN
                first == Allele.UNKNOWN || second == Allele.UNKNOWN -> PARTLY_UNKNOWN
                first == second -> HOMOZYGOUS
                else -> HETEROZYGOUS
            }.toByte()
        }
    }
}
