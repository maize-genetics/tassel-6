package net.maizegenetics.tassel.app.theme

import androidx.compose.ui.graphics.Color
import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.Nucleotide

/**
 * The colors the genotype grid paints calls with. Part of the app theme, so no drawing code names a color, and a
 * pure function, so the mapping can be tested without a window.
 *
 * A call is colored by its two alleles:
 * - both alleles the same nucleotide (`A/A`, `+/+`): that nucleotide's color.
 * - two different known alleles (`A/C`): [heterozygous].
 * - anything with an allele that is unknown, rare, or undefined, including a half call such as `A/N`: [missing].
 *
 * The last rule is deliberate. `Genotype.isHeterozygous` calls `A/N` heterozygous, as TASSEL 5 does, but a grid that
 * painted it as a real heterozygote would hide missing data, so the grid shows it as missing.
 */
data class AlleleColors(
    val a: Color,
    val c: Color,
    val g: Color,
    val t: Color,
    val insert: Color,
    val gap: Color,
    val heterozygous: Color,
    val missing: Color,
) {
    /** The color for one genotype byte. */
    fun of(genotype: Byte): Color {
        val first = Genotype.firstAllele(genotype)
        val second = Genotype.secondAllele(genotype)
        if (!first.isNucleotide() || !second.isNucleotide()) return missing
        return if (first == second) nucleotide(first) else heterozygous
    }

    /** The color for a homozygous call of nucleotide [allele]; [missing] for any other allele code. */
    fun nucleotide(allele: Byte): Color = when (allele) {
        Allele.A -> a
        Allele.C -> c
        Allele.G -> g
        Allele.T -> t
        Allele.INSERT -> insert
        Allele.GAP -> gap
        else -> missing
    }

    /** Label and color for each entry of the grid legend, in the order the legend shows them. */
    fun legend(): List<Pair<String, Color>> = buildList {
        for (allele in NUCLEOTIDES) add(Nucleotide.alleleString(allele) to nucleotide(allele))
        add("het" to heterozygous)
        add("N" to missing)
    }

    private fun Byte.isNucleotide(): Boolean = this in NUCLEOTIDES

    companion object {
        /** The allele codes that get a color of their own, in legend order. */
        private val NUCLEOTIDES = listOf(Allele.A, Allele.C, Allele.G, Allele.T, Allele.INSERT, Allele.GAP)
    }
}
