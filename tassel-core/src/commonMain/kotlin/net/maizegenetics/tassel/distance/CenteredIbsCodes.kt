package net.maizegenetics.tassel.distance

private const val ALLELE_CODE_COUNT = 8
private const val GENOTYPE_BYTE_COUNT = 256

/**
 * How many copies of an allele each diploid genotype carries, indexed by
 * `allele shl 6 or firstAllele shl 3 or secondAllele` with the three low bits of each allele code.
 * The values are the codes packed into an allele-count word: 7 when the genotype or the allele is
 * unknown, 1 for no copies, 2 for one, and 4 for two.
 */
private val COUNT_CODES = ByteArray(512) { index ->
    val allele = index shr 6
    val firstAllele = (index shr 3) and 0x7
    val secondAllele = index and 0x7
    val unknown = 0x7
    when {
        allele == unknown || (firstAllele == unknown && secondAllele == unknown) -> 7
        firstAllele == allele && secondAllele == allele -> 4
        firstAllele == allele || secondAllele == allele -> 2
        else -> 1
    }.toByte()
}

/**
 * Copy-count codes for one allele, indexed by the raw genotype byte. Built once from [COUNT_CODES]:
 * 1 for no copies, 2 for one, 4 for two, and 7 when the genotype is not comparable with that allele.
 */
internal val CODE_BY_ALLELE_AND_GENOTYPE: Array<ByteArray> = Array(ALLELE_CODE_COUNT) { allele ->
    val alleleBits = allele shl 6
    ByteArray(GENOTYPE_BYTE_COUNT) { genotype ->
        COUNT_CODES[alleleBits or ((genotype and 0x70) shr 1) or (genotype and 0x7)]
    }
}

/** Centered dosage `copies - twiceP`, or 0 when the genotype is not comparable at [allele]. */
internal fun centeredDosage(genotype: Int, allele: Int, twiceP: Double): Double {
    val code = CODE_BY_ALLELE_AND_GENOTYPE[allele and 0x7][genotype and 0xFF].toInt()
    return when (code) {
        1 -> 0.0 - twiceP
        2 -> 1.0 - twiceP
        4 -> 2.0 - twiceP
        else -> 0.0
    }
}
