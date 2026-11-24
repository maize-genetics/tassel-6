package net.maizegenetics.tassel.io

import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.AlleleStats
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.position.Position

/**
 * One site as [Vcf.write] prints it, and as the BCF writer will encode it.
 *
 * [alleles] is the VCF index order: index 0 is [reference], then [alternates]. It is empty when the site has no
 * called allele, and [alleleCount] is then 1, which is the BCF `n_allele` for a site written as `.` / `.`.
 * [calls] holds two indexes per taxon, the high nibble then the low nibble. [MISSING] is the `.` in a `GT` field.
 *
 * A gap that is one of [alleles] keeps its index. TASSEL 5 has a later branch that rewrites a gap call as `.` and
 * decrements later indexes. That branch does not run once the variants have been resolved, and the `n_prefix`
 * golden is `1/1` for the gap, so it is not ported.
 *
 * @see net.maizegenetics.dna.snp.ExportUtils.writeToVCF
 */
internal class VcfSiteRecord(
    val reference: String,
    val alternates: List<String>,
    val alleles: ByteArray,
    val calls: ByteArray,
) {
    /** BCF `n_allele`. A site with no alleles is still one allele, the missing one. */
    val alleleCount: Int get() = if (alleles.isEmpty()) 1 else alleles.size

    init {
        require(calls.size % 2 == 0) { "VCF export calls must be two per taxon: ${calls.size}" }
    }

    companion object {
        /** A call whose allele is not in [alleles]. [Vcf.write] prints it as `.`. */
        const val MISSING: Byte = -1
    }
}

/**
 * The `##` lines [Vcf.write] emits. TASSEL 5 writes these even when the record carries only `GT`.
 * There is no `##SAMPLE` line and no `##contig` line.
 */
internal val VCF_EXPORT_HEADER: List<String> = listOf(
    "##fileformat=VCFv4.0",
    "##Tassel=<ID=GenotypeTable,Version=5,Description=\"Reference allele is not known. " +
        "The major allele was used as reference allele\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference " +
        "and alternate alleles in the order listed\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">",
    "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">",
    "##FORMAT=<ID=PL,Number=.,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for " +
        "AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">",
    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
    "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">",
)

/**
 * REF, ALT, and GT indexes for one site, using TASSEL 5's unknown-reference rules.
 *
 * [stats] supplies the alleles in frequency order. [knownVariants] is the site's alleles column, or the REF/ALT
 * text a reader stored. [Position] has no reference allele, so `resolveRefSorted` is not applied: with an unknown
 * reference that function copies this frequency order through unchanged.
 *
 * @see net.maizegenetics.dna.snp.ExportUtils.writeToVCF
 * @see net.maizegenetics.dna.snp.io.VCFUtil.indelInKnownVariant
 */
internal fun vcfSiteRecord(knownVariants: List<String>, stats: AlleleStats, genotypes: ByteArray): VcfSiteRecord {
    val resolved = resolveAlleles(knownVariants, stats.alleles)
    if (resolved.alleles.isEmpty()) {
        return VcfSiteRecord(
            reference = ".",
            alternates = emptyList(),
            alleles = ByteArray(0),
            calls = ByteArray(genotypes.size * 2) { VcfSiteRecord.MISSING },
        )
    }
    val alleles = ByteArray(resolved.alleles.size) { resolved.alleles[it].toByte() }
    val redirect = IntArray(ALLELE_CODES) { -1 }
    for (index in alleles.indices) redirect[alleles[index].toInt() and 0xF] = index
    val calls = ByteArray(genotypes.size * 2)
    for (taxon in genotypes.indices) {
        val genotype = genotypes[taxon]
        calls[taxon * 2] = callIndex(redirect, Genotype.firstAllele(genotype))
        calls[taxon * 2 + 1] = callIndex(redirect, Genotype.secondAllele(genotype))
    }
    return VcfSiteRecord(
        reference = resolved.texts[0],
        alternates = resolved.texts.drop(1),
        alleles = alleles,
        calls = calls,
    )
}

/** Appends one VCF record, including its trailing `\n`. */
internal fun appendVcfRecord(line: StringBuilder, position: Position, record: VcfSiteRecord) {
    line.append(position.chromosome.name).append('\t')
    line.append(position.position).append('\t')
    line.append(position.name).append('\t')
    line.append(record.reference).append('\t')
    if (record.alternates.isEmpty()) {
        line.append('.')
    } else {
        record.alternates.forEachIndexed { index, alternate ->
            if (index > 0) line.append(',')
            line.append(alternate)
        }
    }
    line.append("\t.\tPASS\t.\tGT")
    for (call in record.calls.indices step 2) {
        line.append('\t').appendCall(record.calls[call]).append('/').appendCall(record.calls[call + 1])
    }
    line.append('\n')
}

private fun StringBuilder.appendCall(index: Byte): StringBuilder =
    if (index == VcfSiteRecord.MISSING) append('.') else append(index.toInt())

/**
 * Allele codes and the REF/ALT text in the same order.
 *
 * Codes are 0 to 15, the low four bits of a nucleotide allele.
 */
private class ResolvedAlleles(val alleles: List<Int>, val texts: List<String>)

/**
 * Frequency order when the site lists no variants. A gap is written `N`, and every other allele is prefixed
 * with `N`, which is the no-known-variant branch of `writeToVCF`.
 */
private fun allelesFromFrequency(present: List<Byte>): ResolvedAlleles {
    val gapAt = present.indexOf(Allele.GAP)
    val texts = present.mapIndexed { index, allele ->
        val text = Nucleotide.alleleString(allele)
        when {
            gapAt < 0 -> text
            index == gapAt -> "N"
            else -> "N$text"
        }
    }
    return ResolvedAlleles(present.map { it.toInt() and 0xF }, texts)
}

/**
 * Merges [present] with [original], which is what `writeToVCF` does when `knownVariants` is not empty.
 *
 * A variant longer than one character is an indel: its allele is the second character, and a one-character
 * variant in that set is a gap. The written text stays the original string (`AT` stays `AT`). Otherwise the
 * allele is the first character, and only the last variant that starts with `N` is removed. Present alleles
 * stay in frequency order. A known allele that no genotype carries is appended in known-variant order. When a
 * gap is present and the known variants are not indels, that appended text is prefixed with `N`, and a gap
 * that was not in the known list is written `N`.
 */
private fun allelesFromKnownVariants(original: List<String>, present: List<Byte>): ResolvedAlleles {
    val tempAlleles = ArrayList<Int>()
    val variants: List<String>
    if (original.any { it.length > 1 }) {
        for (variant in original) {
            tempAlleles += if (variant.length > 1) alleleCode(variant[1]) else alleleCode('-')
        }
        variants = original
    } else {
        var lastN = -1
        for ((index, variant) in original.withIndex()) {
            require(variant.isNotEmpty()) { "VCF export variant is empty" }
            if (variant[0] == 'N') {
                lastN = index
            } else {
                tempAlleles += alleleCode(variant[0])
            }
        }
        variants = if (lastN < 0) original else original.filterIndexed { index, _ -> index != lastN }
    }

    val indelsInKnownVariants = variants.any { it.length > 1 }
    val indelsExist = indelsInKnownVariants || present.any { it == Allele.GAP }
    val prefixed = if (indelsExist && !indelsInKnownVariants) variants.map { "N$it" } else variants

    val ordered = ArrayList<Int>()
    val textByAllele = mutableMapOf<Int, String>()
    for (allele in present) {
        val code = allele.toInt() and 0xF
        ordered += code
        val knownAt = tempAlleles.indexOf(code)
        textByAllele[code] = if (knownAt < 0) {
            if (indelsExist && allele == Allele.GAP) "N" else Nucleotide.alleleString(allele)
        } else {
            variants[knownAt]
        }
    }
    for (index in tempAlleles.indices) {
        val code = tempAlleles[index]
        if (code !in ordered) {
            ordered += code
            textByAllele[code] = prefixed[index]
        }
    }
    return ResolvedAlleles(ordered, ordered.map { textByAllele.getValue(it) })
}

private fun resolveAlleles(knownVariants: List<String>, present: List<Byte>): ResolvedAlleles =
    if (knownVariants.isEmpty()) allelesFromFrequency(present) else allelesFromKnownVariants(knownVariants, present)

private fun alleleCode(base: Char): Int = Nucleotide.allele(base).toInt() and 0xF

private fun callIndex(redirect: IntArray, allele: Byte): Byte {
    val index = redirect[allele.toInt() and 0xF]
    return if (index < 0) VcfSiteRecord.MISSING else index.toByte()
}

private const val ALLELE_CODES = 16
