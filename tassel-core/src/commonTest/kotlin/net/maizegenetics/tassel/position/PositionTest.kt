package net.maizegenetics.tassel.position

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNotEquals
import kotlin.test.assertNull
import kotlin.test.assertTrue

class PositionTest {
    /** Names and numbers that TASSEL 5 `v5.2.98` `net.maizegenetics.dna.map.Chromosome` gives each input name. */
    private val tassel5Names = listOf(
        "1" to "1", "chr2" to "2", "Chr10" to "10", "CHROMOSOME3" to "3", "chr1 extra" to "1", " 5 " to "5",
        "10" to "10", "2" to "2", "X" to "X", "x" to "X", "Y" to "Y", "scaffold_1" to "SCAFFOLD_1", "Mt" to "MT",
        "Pt" to "PT", "01" to "01", "1A" to "1A", "10A" to "10A", "UNKNOWN" to "UNKNOWN", "0" to "0", "-1" to "-1",
        "+3" to "+3",
    )

    /** Sign of TASSEL 5's `compareTo` for every pair of [tassel5Names], one row per name. */
    private val tassel5Comparisons = (
        "0---0---------0---++- +0--+--0------++--++- ++0+++0+------++--+++ ++-0+--+------++--++0 " +
            "0---0---------0---++- ++-++0-+------++--+++ ++0+++0+------++--+++ +0--+--0------++--++- " +
            "++++++++00-++++++++++ ++++++++00-++++++++++ ++++++++++0++++++++++ ++++++++---0+++++-+++ " +
            "++++++++----0-+++-+++ ++++++++----+0+++-+++ 0---0---------0---++- +---+---------+0--+++ " +
            "++++++++------++0-+++ ++++++++---++++++0+++ ------------------0+- -------------------0- " +
            "++-0+--+------+---++0"
        ).split(' ')

    @Test
    fun chromosomeNamesMatchTassel5() {
        for ((input, name) in tassel5Names) assertEquals(name, Chromosome.of(input).name, input)
        assertEquals(1, Chromosome.of("01").number)
        assertEquals(3, Chromosome.of("+3").number)
        assertNull(Chromosome.of("X").number)
        assertFailsWith<IllegalArgumentException> { Chromosome.of("") }
    }

    @Test
    fun chromosomeOrderMatchesTassel5() {
        val chromosomes = tassel5Names.map { Chromosome.of(it.first) }
        for (i in chromosomes.indices) {
            for (j in chromosomes.indices) {
                val sign = chromosomes[i].compareTo(chromosomes[j]).let {
                    if (it <
                        0
                    ) {
                        '-'
                    } else if (it > 0) {
                        '+'
                    } else {
                        '0'
                    }
                }
                assertEquals(tassel5Comparisons[i][j], sign, "${tassel5Names[i].first} vs ${tassel5Names[j].first}")
            }
        }
    }

    @Test
    fun chromosomesAreEqualByName() {
        assertEquals(Chromosome.of("chr1"), Chromosome.of("1"))
        assertNotEquals(Chromosome.of("01"), Chromosome.of("1"))
    }

    @Test
    fun defaultSiteNameIsDropped() {
        val chr1 = Chromosome.of(1)
        assertNull(Position(chr1, 100, snpId = "S1_100").snpId)
        assertEquals("S1_100", Position(chr1, 100).name)
        assertEquals("PZB00859.1", Position(chr1, 157104, snpId = "PZB00859.1").name)
        assertEquals(Position(chr1, 100), Position(chr1, 100, snpId = "S1_100"))
        assertEquals(Strand.PLUS, Position(chr1, 100).strand)
    }

    @Test
    fun positionsSortByChromosomePositionStrandAndName() {
        val chr1 = Chromosome.of(1)
        val sorted = listOf(
            Position(chr1, 5, snpId = "b"),
            Position(chr1, 5, snpId = "a"),
            Position(Chromosome.of(2), 1),
            Position(chr1, 5, strand = Strand.MINUS),
            Position(chr1, 4),
        ).sorted()
        assertEquals(listOf("S1_4", "S1_5", "a", "b", "S2_1"), sorted.map { it.name })
        assertEquals(Strand.MINUS, sorted[1].strand)
        assertEquals(Strand.MINUS, Strand.of("-"))
        assertFailsWith<IllegalArgumentException> { Strand.of("?") }
    }

    @Test
    fun positionListFindsChromosomesAndSites() {
        val chr1 = Chromosome.of(1)
        val chr2 = Chromosome.of(2)
        val list = PositionList(
            listOf(
                Position(chr1, 10),
                Position(chr1, 20, snpId = "a"),
                Position(chr1, 20, snpId = "b"),
                Position(chr1, 30),
                Position(chr2, 5),
            ),
        )
        assertEquals(listOf(chr1, chr2), list.chromosomes)
        assertEquals(0..3, list.siteRange(chr1))
        assertEquals(4..4, list.siteRange(chr2))
        assertTrue(list.siteRange(Chromosome.of(3)).isEmpty())
        assertEquals(1, list.siteOfPhysicalPosition(chr1, 20))
        assertEquals(0, list.siteOfPhysicalPosition(chr1, 10))
        assertEquals(-1, list.siteOfPhysicalPosition(chr1, 9))
        assertEquals(-4, list.siteOfPhysicalPosition(chr1, 25))
        assertEquals(-5, list.siteOfPhysicalPosition(chr1, 31))
        assertEquals(-6, list.siteOfPhysicalPosition(chr2, 6))
        assertEquals(Int.MIN_VALUE, list.siteOfPhysicalPosition(Chromosome.of(3), 1))
    }

    @Test
    fun positionListRejectsUnsortedSites() {
        val chr1 = Chromosome.of(1)
        assertFailsWith<IllegalArgumentException> { PositionList(listOf(Position(chr1, 2), Position(chr1, 1))) }
        assertTrue(PositionList(emptyList()).chromosomes.isEmpty())
    }
}
