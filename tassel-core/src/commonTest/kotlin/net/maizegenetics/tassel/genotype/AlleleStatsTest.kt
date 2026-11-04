package net.maizegenetics.tassel.genotype

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Allele statistics compared with TASSEL 5 `v5.2.98`. Expected values were printed by a small Java program calling
 * `AlleleFreqCache.allelesSortedByFrequencyNucleotide`, `allelesSortedByFrequencyAndCountsNucleotide` (whose
 * `Stats` gives the frequencies and proportions), and `AlleleFreqCache.proportionHeterozygous(byte[])`.
 */
class AlleleStatsTest {
    private class Case(
        val genotypes: ByteArray,
        val alleles: List<Int>,
        val counts: List<Int>,
        val majorFrequency: Double,
        val minorFrequency: Double,
        val allMinor: Int,
        val proportionHeterozygous: Double,
        val proportionNotMissing: Double,
        val hasIndel: Boolean,
        /** TASSEL 5 `AlleleFreqCache.proportionHeterozygous(byte[])`: differing allele codes among non-`N/N`. */
        val heterozygousAmongNonMissing: Double,
    )

    private fun iupac(codes: String) = ByteArray(codes.length) { Nucleotide.genotype(codes[it]) }

    private fun raw(vararg genotypes: Int) = ByteArray(genotypes.size) { genotypes[it].toByte() }

    private val cases = listOf(
        Case(iupac("N"), listOf(), listOf(), 0.0, 0.0, 0, 0.0, 0.0, false, 0.0),
        Case(iupac("NNNN"), listOf(), listOf(), 0.0, 0.0, 0, 0.0, 0.0, false, 0.0),
        Case(iupac("A"), listOf(0), listOf(2), 1.0, 0.0, 0, 0.0, 1.0, false, 0.0),
        Case(iupac("AAAA"), listOf(0), listOf(8), 1.0, 0.0, 0, 0.0, 1.0, false, 0.0),
        Case(iupac("ACAC"), listOf(0, 1), listOf(4, 4), 0.5, 0.5, 4, 0.0, 1.0, false, 0.0),
        Case(iupac("CACA"), listOf(0, 1), listOf(4, 4), 0.5, 0.5, 4, 0.0, 1.0, false, 0.0),
        Case(iupac("AMCN"), listOf(0, 1), listOf(3, 3), 0.5, 0.5, 3, 0.25, 0.75, false, 0.3333333333333333),
        Case(iupac("RRRR"), listOf(0, 2), listOf(4, 4), 0.5, 0.5, 4, 1.0, 1.0, false, 1.0),
        Case(
            iupac("AAGGTT"), listOf(0, 2, 3), listOf(4, 4, 4), 0.3333333333333333, 0.3333333333333333, 8, 0.0, 1.0,
            false, 0.0,
        ),
        Case(iupac("AGTGGTTA"), listOf(2, 3, 0), listOf(6, 6, 4), 0.375, 0.375, 10, 0.0, 1.0, false, 0.0),
        Case(iupac("+-0-"), listOf(5, 4), listOf(5, 3), 0.625, 0.375, 3, 0.25, 1.0, true, 0.25),
        Case(iupac("AACCNNMM"), listOf(0, 1), listOf(6, 6), 0.5, 0.5, 6, 0.25, 0.75, false, 0.3333333333333333),
        Case(iupac("GGGCCCAAATTT"), listOf(0, 1, 2, 3), listOf(6, 6, 6, 6), 0.25, 0.25, 18, 0.0, 1.0, false, 0.0),
        Case(
            iupac("ACGT+-"), listOf(0, 1, 2, 3, 4, 5), listOf(2, 2, 2, 2, 2, 2), 0.16666666666666666,
            0.16666666666666666, 10, 0.0, 1.0, true, 0.0,
        ),
        Case(
            iupac("TTTTTA"), listOf(3, 0), listOf(10, 2), 0.8333333333333334, 0.16666666666666666, 2, 0.0, 1.0, false,
            0.0,
        ),
        Case(iupac("AZXN"), listOf(0), listOf(2), 1.0, 0.0, 0, 0.0, 0.5, false, 0.0),
        Case(
            raw(0x0F, 0xF0, 0x00, 0x01), listOf(0, 1), listOf(5, 1), 0.8333333333333334, 0.16666666666666666, 1,
            0.25, 1.0, false, 0.75,
        ),
        Case(raw(0x0E, 0xE0, 0xEE, 0x66, 0x06, 0x11), listOf(0, 1), listOf(3, 2), 0.6, 0.4, 2, 0.5, 1.0, false, 0.5),
        Case(
            raw(0x12, 0x21, 0x22, 0x2F, 0xFF, 0x77), listOf(2, 1), listOf(5, 2), 0.7142857142857143,
            0.2857142857142857, 2, 0.3333333333333333, 0.8333333333333334, false, 0.6,
        ),
    )

    @Test
    fun matchesTassel5() {
        for (case in cases) {
            val name = case.genotypes.joinToString(",") { (it.toInt() and 0xFF).toString(16) }
            val stats = AlleleStats.of(case.genotypes)
            assertEquals(case.alleles.map { it.toByte() }, stats.alleles, name)
            assertEquals(case.counts, stats.alleleCounts, name)
            assertEquals(case.alleles.getOrElse(0) { 15 }.toByte(), stats.majorAllele, name)
            assertEquals(case.alleles.getOrElse(1) { 15 }.toByte(), stats.minorAllele, name)
            assertEquals(case.majorFrequency, stats.majorAlleleFrequency, name)
            assertEquals(case.minorFrequency, stats.minorAlleleFrequency, name)
            assertEquals(case.counts.sum(), stats.alleleGameteCount, name)
            assertEquals(case.allMinor, stats.allMinorAlleleCount, name)
            assertEquals(case.proportionHeterozygous, stats.proportionHeterozygous, name)
            assertEquals(case.proportionNotMissing, stats.proportionNotMissing, name)
            assertEquals(case.hasIndel, stats.hasIndel, name)
            val heterozygousAmongNonMissing =
                if (stats.nonMissingCount == 0) 0.0 else differingCount(stats).toDouble() / stats.nonMissingCount
            assertEquals(case.heterozygousAmongNonMissing, heterozygousAmongNonMissing, name)
        }
    }

    @Test
    fun genotypeGroupsAddUp() {
        for (case in cases) {
            val stats = AlleleStats.of(case.genotypes)
            val groups = stats.unknownCount + stats.partlyUnknownCount + stats.heterozygousCount + stats.homozygousCount
            assertEquals(case.genotypes.size, groups)
            assertEquals(case.genotypes.count { it != Genotype.UNKNOWN }, stats.nonMissingCount)
            val knownGametes = case.genotypes.sumOf { genotype ->
                listOf(Genotype.firstAllele(genotype), Genotype.secondAllele(genotype)).count { it != Allele.UNKNOWN }
            }
            assertEquals(knownGametes, stats.nonMissingGameteCount)
            assertEquals(case.genotypes.count(Genotype::isHeterozygous), differingCount(stats))
        }
    }

    @Test
    fun emptyInputGivesNaNProportions() {
        val stats = AlleleStats.of(ByteArray(0))
        assertEquals(Allele.UNKNOWN, stats.majorAllele)
        assertEquals(0.0, stats.minorAlleleFrequency)
        assertTrue(stats.proportionHeterozygous.isNaN())
        assertTrue(stats.proportionNotMissing.isNaN())
    }

    @Test
    fun countsOnlyTheGivenRange() {
        val genotypes = iupac("TTAACCGG")
        val stats = AlleleStats.of(genotypes, fromIndex = 2, toIndex = 6)
        assertEquals(listOf(Allele.A, Allele.C), stats.alleles)
        assertEquals(4, stats.numberOfGenotypes)
    }

    private fun differingCount(stats: AlleleStats) = stats.heterozygousCount + stats.partlyUnknownCount
}
