package net.maizegenetics.tassel.io

import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.Nucleotide

/**
 * The variants a site's calls index into, from its `REF` and `ALT` columns as TASSEL 5 builds them.
 *
 * Symbolic alleles become one character (`<INS>` to `+`, `<DEL>` to `-`, `*` to `N`), and a missing `ALT` leaves
 * the reference as the only variant. [splitVariants], [allelesOf], and [call] are the rest of the same rules,
 * moved out of [VcfRecordParser] so the BCF reader can apply them. [call] still reads the ASCII field.
 *
 * @see net.maizegenetics.dna.snp.io.BuilderFromVCF
 */
internal fun variantsOf(reference: String, alternate: String): String {
    if (alternate == ".") return reference
    // TASSEL 5's substitutions, in its order: <INS> and any I become +, <DEL> and any D become -, * becomes N.
    return (reference + "/" + alternate).replace(',', '/')
        .replace("<INS>", "+").replace('I', '+')
        .replace("<DEL>", "-").replace('D', '-')
        .replace("*", "N")
}

/** Splits the variants on `/` as Java's `String.split` does, which drops trailing empty parts. */
internal fun splitVariants(variants: String, site: String): List<String> {
    val parts = variants.split('/')
    var size = parts.size
    while (size > 0 && parts[size - 1].isEmpty()) size--
    require(size > 0) { "VCF site $site has no REF or ALT allele: [$variants]" }
    return if (size == parts.size) parts else parts.subList(0, size)
}

/**
 * One allele code per variant. A site is an indel when the reference is longer than one character or some
 * alternate is longer than the reference; TASSEL 5 then drops the shared first base and calls what is left of a
 * one-character variant a gap, so `AT`/`A` becomes `T` and `-`.
 */
internal fun allelesOf(variants: List<String>, site: String): ByteArray {
    val indel = variants[0].length > 1 || variants.drop(1).any { it.length > variants[0].length }
    return ByteArray(variants.size) { index ->
        val variant = variants[index]
        val text = when {
            !indel -> variant
            variant.length > 1 -> variant.substring(1)
            else -> "-"
        }
        require(text.isNotEmpty()) { "VCF site $site has an empty allele in [${variants.joinToString("/")}]" }
        Nucleotide.allele(text[0])
    }
}

/**
 * One taxon's call, from the `GT` field of `window[from until to]`.
 *
 * The separator between the two allele indexes is ignored, so `0|1` and `0/1` give the same genotype. TASSEL 5
 * also reads only the characters at offsets 0 and 2, which is why an allele index above 9 is not supported.
 */
internal fun call(window: ByteArray, from: Int, to: Int, alleles: ByteArray, site: String): Byte {
    require(from < to) { "VCF site $site has an empty genotype column" }
    // TASSEL 5 treats any column starting with a dot as missing, which covers ./. and .|. and a bare dot.
    if (window[from] == DOT) return Genotype.UNKNOWN
    var end = from
    while (end < to && window[end] != COLON) end++
    val first: Int
    val second: Int
    if (end - from == 1) {
        first = digit(window[from])
        second = first
    } else {
        require(end - from >= 3) {
            // TASSEL 5 indexes character 2 unguarded here and dies with a string-index error instead.
            "VCF site $site has a genotype that is not diploid: [${window.decodeToString(from, end)}]"
        }
        first = digit(window[from])
        second = digit(window[from + 2])
    }
    require(first <= alleles.size - 1 && second <= alleles.size - 1) {
        "VCF site $site has a genotype naming an allele it does not list: " +
            "[${window.decodeToString(from, end)}] with ${alleles.size} variants"
    }
    // A half-missing call such as 0/. gives a negative index, which TASSEL 5 checks only after the bound above.
    return if (first < 0 || second < 0) Genotype.UNKNOWN else Genotype.of(alleles[first], alleles[second])
}

private fun digit(byte: Byte): Int = (byte.toInt() and 0xFF) - '0'.code

private val DOT = '.'.code.toByte()
private val COLON = ':'.code.toByte()
