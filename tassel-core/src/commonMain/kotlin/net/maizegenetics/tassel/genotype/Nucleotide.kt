package net.maizegenetics.tassel.genotype

/**
 * Text encoding of nucleotide alleles and genotypes, with the same tables as TASSEL 5.
 *
 * Genotypes are written as one IUPAC character (`A`, `R` for A/G, `N` for missing, `0` for an indel heterozygote) or
 * as two allele characters (`AG`, `A-`, `NN`). Allele characters are `A C G T + -`, `N` or `X` for missing, and `Z` for
 * rare.
 *
 * @see net.maizegenetics.dna.snp.NucleotideAlignmentConstants
 */
object Nucleotide {
    /** Text for each allele code 0 to 15: `A C G T + -`, `X` for codes 6 to 13, `Z` for rare, and `N` for unknown. */
    val alleleStrings: List<String> = List(16) { code ->
        when (code) {
            in 0..5 -> ALLELE_CHARS[code].toString()
            Allele.RARE.toInt() -> "Z"
            Allele.UNKNOWN.toInt() -> "N"
            else -> UNDEFINED_STRING
        }
    }

    /**
     * Genotype for a one-character code, such as `A`, `R`, or `N`.
     *
     * Like TASSEL 5, characters without a genotype meaning are not rejected: they read whatever the shared lookup
     * table holds at that index, usually `0x66` (undefined).
     *
     * @throws IllegalArgumentException if [code] is outside 0 to 255
     */
    fun genotype(code: Char): Byte {
        require(code.code < 256) { "Not a nucleotide genotype character: $code" }
        return DIPLOID_TABLE[code.code]
    }

    /**
     * Genotype for a one- or two-character code, such as `R` or `AG`. Two-character codes keep allele order.
     *
     * Two-character codes go through TASSEL 5's hashed lookup table, so an invalid pair can return a valid-looking
     * genotype instead of failing; this matches TASSEL 5.
     *
     * @throws IllegalArgumentException if [code] is not one or two characters long
     */
    fun genotype(code: String): Byte = when (code.length) {
        1 -> genotype(code[0])
        2 -> DIPLOID_TABLE[pairIndex(code[0], code[1])]
        else -> throw IllegalArgumentException("Nucleotide genotype must be one or two characters: $code")
    }

    /**
     * Allele code for one character: `A C G T` (either case), `+`, `-`, or `N`/`n` for [Allele.UNKNOWN]. Any other
     * character gives [Allele.UNDEFINED].
     *
     * @throws IllegalArgumentException if [code] is outside 0 to 255
     */
    fun allele(code: Char): Byte {
        require(code.code < 256) { "Not a nucleotide allele character: $code" }
        return ALLELE_TABLE[code.code]
    }

    /** One-character IUPAC code for a genotype, or `X` if the genotype has none. */
    fun iupac(genotype: Byte): Char = IUPAC_TABLE[genotype.toInt() and 0xFF]

    /** Text for an allele code, such as `A` or `+`. */
    fun alleleString(allele: Byte): String = alleleStrings[allele.toInt() and 0xF]

    /** Genotype as two allele strings separated by a colon, such as `A:G`. */
    fun genotypeString(genotype: Byte): String =
        alleleString(Genotype.firstAllele(genotype)) + ":" + alleleString(Genotype.secondAllele(genotype))

    /** Complement of an allele code: A and T swap, C and G swap, and other codes are unchanged. */
    fun complement(allele: Byte): Byte = when (allele) {
        Allele.A -> Allele.T
        Allele.T -> Allele.A
        Allele.C -> Allele.G
        Allele.G -> Allele.C
        else -> allele
    }

    /** Complement of both alleles of a genotype, keeping their order. */
    fun complementGenotype(genotype: Byte): Byte =
        Genotype.of(complement(Genotype.firstAllele(genotype)), complement(Genotype.secondAllele(genotype)))

    /** Complement of an IUPAC genotype character, such as `R` to `Y`; `X` for characters without one. */
    fun complementIupac(code: Char): Char = if (code.code < 256) IUPAC_COMPLEMENT_TABLE[code.code] else 'X'

    private const val ALLELE_CHARS = "ACGT+-"
    private const val UNDEFINED_STRING = "X"
    private const val UNDEFINED_GENOTYPE: Byte = 0x66

    /** Characters accepted in two-character genotype codes, and the allele code each one stands for. */
    private val PAIR_ALLELES = mapOf(
        'A' to Allele.A, 'C' to Allele.C, 'G' to Allele.G, 'T' to Allele.T, '+' to Allele.INSERT, '-' to Allele.GAP,
        'N' to Allele.UNKNOWN, 'X' to Allele.UNKNOWN, 'Z' to Allele.RARE,
    )

    private val SINGLE_GENOTYPES = mapOf(
        'A' to 0x00, 'C' to 0x11, 'G' to 0x22, 'T' to 0x33, '+' to 0x44, '-' to 0x55, 'Z' to 0xEE, 'N' to 0xFF,
        'X' to 0xFF, 'R' to 0x02, 'Y' to 0x13, 'S' to 0x21, 'W' to 0x03, 'K' to 0x23, 'M' to 0x01, '0' to 0x54,
    )

    /** TASSEL 5's hash from a two-character code to a slot in the same table that one-character codes index. */
    private fun pairIndex(first: Char, second: Char): Int =
        (((second.code shl 2) xor PAIR_MASK_SECOND) xor (first.code and PAIR_MASK_FIRST)) and 0xFF

    private const val PAIR_MASK_SECOND = -127
    private const val PAIR_MASK_FIRST = 0x2F

    private val DIPLOID_TABLE = ByteArray(256) { UNDEFINED_GENOTYPE }.also { table ->
        for ((first, firstAllele) in PAIR_ALLELES) {
            for ((second, secondAllele) in PAIR_ALLELES) {
                table[pairIndex(first, second)] = Genotype.of(firstAllele, secondAllele)
            }
        }
        for ((code, genotype) in SINGLE_GENOTYPES) table[code.code] = genotype.toByte()
    }

    private val ALLELE_TABLE = ByteArray(256) { Allele.UNDEFINED }.also { table ->
        ALLELE_CHARS.forEachIndexed { allele, code ->
            table[code.code] = allele.toByte()
            table[code.lowercaseChar().code] = allele.toByte()
        }
        table['N'.code] = Allele.UNKNOWN
        table['n'.code] = Allele.UNKNOWN
    }

    /** IUPAC codes for heterozygotes of two different A, C, G, T alleles, keyed by `smaller * 4 + larger`. */
    private val IUPAC_HETEROZYGOTES = mapOf(1 to 'M', 2 to 'R', 3 to 'W', 6 to 'S', 7 to 'Y', 11 to 'K')

    private val IUPAC_TABLE = CharArray(256) { iupacOf(it shr 4, it and 0xF) }

    private val IUPAC_COMPLEMENT_TABLE = CharArray(256) { 'X' }.also { table ->
        for (genotype in 0..255) {
            val code = IUPAC_TABLE[genotype]
            if (code != 'X') table[code.code] = iupac(complementGenotype(genotype.toByte()))
        }
    }

    private fun iupacOf(first: Int, second: Int): Char {
        val rare = Allele.RARE.toInt()
        fun valid(allele: Int) = allele < Allele.NUCLEOTIDE_STATES || allele >= rare
        fun text(allele: Int) = alleleString(allele.toByte())[0]
        return when {
            !valid(first) || !valid(second) -> 'X'
            first == second -> text(first)
            first >= rare && second >= rare -> 'N'
            first >= rare -> text(second)
            second >= rare -> text(first)
            first >= Allele.INSERT || second >= Allele.INSERT -> '0'
            else -> IUPAC_HETEROZYGOTES.getValue(minOf(first, second) * 4 + maxOf(first, second))
        }
    }
}
