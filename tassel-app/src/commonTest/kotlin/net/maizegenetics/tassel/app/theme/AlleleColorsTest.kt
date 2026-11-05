package net.maizegenetics.tassel.app.theme

import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.Nucleotide
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotEquals
import kotlin.test.assertTrue

/** The allele-to-color mapping the grid paints with. A pure function, so it needs no window to test. */
class AlleleColorsTest {
    private val colors = LightAlleleColors

    /** Every color the mapping can return. */
    private val allColors = setOf(
        colors.a,
        colors.c,
        colors.g,
        colors.t,
        colors.insert,
        colors.gap,
        colors.heterozygous,
        colors.missing,
    )

    @Test
    fun homozygousCallsTakeTheirNucleotideColor() {
        assertEquals(colors.a, colors.of(Nucleotide.genotype('A')))
        assertEquals(colors.c, colors.of(Nucleotide.genotype('C')))
        assertEquals(colors.g, colors.of(Nucleotide.genotype('G')))
        assertEquals(colors.t, colors.of(Nucleotide.genotype('T')))
        assertEquals(colors.insert, colors.of(Nucleotide.genotype('+')))
        assertEquals(colors.gap, colors.of(Nucleotide.genotype('-')))
    }

    @Test
    fun theFourBasesAreToldApart() {
        val distinct = setOf(colors.a, colors.c, colors.g, colors.t, colors.insert, colors.gap, colors.heterozygous)
        assertEquals(7, distinct.size, "two allele colors are the same, so the grid cannot be read")
        assertTrue(colors.missing !in distinct)
    }

    @Test
    fun twoKnownAllelesAreHeterozygous() {
        // M is the IUPAC code for A/C; R is A/G.
        assertEquals(colors.heterozygous, colors.of(Nucleotide.genotype('M')))
        assertEquals(colors.heterozygous, colors.of(Nucleotide.genotype('R')))
        assertEquals(colors.heterozygous, colors.of(Genotype.of(Allele.A, Allele.GAP)))
    }

    @Test
    fun aMissingCallIsMissing() {
        assertEquals(colors.missing, colors.of(Genotype.UNKNOWN))
        assertEquals(colors.missing, colors.of(Nucleotide.genotype('N')))
    }

    @Test
    fun aHalfMissingCallIsMissingRatherThanHeterozygous() {
        val halfMissing = Genotype.of(Allele.A, Allele.UNKNOWN)
        // TASSEL 5, and Genotype.isHeterozygous with it, calls this heterozygous. Painting it as a heterozygote
        // would hide missing data, so the grid does not.
        assertTrue(Genotype.isHeterozygous(halfMissing))
        assertEquals(colors.missing, colors.of(halfMissing))
        assertNotEquals(colors.heterozygous, colors.of(halfMissing))
    }

    @Test
    fun rareAndUndefinedAllelesAreMissing() {
        assertEquals(colors.missing, colors.of(Genotype.RARE))
        assertEquals(colors.missing, colors.of(Genotype.of(Allele.UNDEFINED, Allele.UNDEFINED)))
        assertEquals(colors.missing, colors.of(Genotype.of(Allele.A, Allele.RARE)))
    }

    @Test
    fun everyGenotypeByteGetsAColor() {
        // The grid paints whatever a reader produced, so no byte may fall through the mapping.
        for (value in Byte.MIN_VALUE..Byte.MAX_VALUE) {
            val color = colors.of(value.toByte())
            assertTrue(color in allColors, "no color for genotype byte $value")
        }
    }

    @Test
    fun nonNucleotideAllelesFallBackToMissing() {
        assertEquals(colors.missing, colors.nucleotide(Allele.UNKNOWN))
        assertEquals(colors.missing, colors.nucleotide(Allele.RARE))
        assertEquals(colors.missing, colors.nucleotide(Allele.UNDEFINED))
    }

    @Test
    fun theLegendNamesEveryColorTheGridUses() {
        val legend = LightAlleleColors.legend()

        assertEquals(listOf("A", "C", "G", "T", "+", "-", "het", "N"), legend.map { it.first })
        assertEquals(allColors, legend.map { it.second }.toSet())
    }

    @Test
    fun theDarkPaletteCoversTheSameCases() {
        assertEquals(LightAlleleColors.legend().map { it.first }, DarkAlleleColors.legend().map { it.first })
        assertEquals(allColors.size, DarkAlleleColors.legend().map { it.second }.toSet().size)
    }
}
