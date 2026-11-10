package net.maizegenetics.tassel.spike.gram

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.distance.DistanceMatrix
import net.maizegenetics.tassel.distance.centeredIbs
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class GramSpikeTest {
    @Test
    fun parkMillerFirstDraw() {
        assertEquals(3_909_951, ParkMiller(PARK_MILLER_SEED).next())
    }

    @Test
    fun dosageCodesStayInRangeAndRepeat() {
        val first = generateDosages(6, 12)
        val second = generateDosages(6, 12)
        assertTrue(first.contentEquals(second))
        assertTrue(first.all { it.toInt() in 0..MISSING_CODE })
    }

    @Test
    fun triangleIndexPacksTheUpperTriangle() {
        assertEquals(0, triangleIndex(0, 0, 3))
        assertEquals(1, triangleIndex(0, 1, 3))
        assertEquals(2, triangleIndex(0, 2, 3))
        assertEquals(3, triangleIndex(1, 1, 3))
        assertEquals(4, triangleIndex(1, 2, 3))
        assertEquals(5, triangleIndex(2, 2, 3))
    }

    @Test
    fun handDosageGram() {
        val gram = dosageGram(byteArrayOf(2, 0, 1, 2), taxa = 2, sites = 2, panelColumns = 2)
        assertEquals(5.0, gram[triangleIndex(0, 0, 2)])
        assertEquals(2.0, gram[triangleIndex(0, 1, 2)])
        assertEquals(4.0, gram[triangleIndex(1, 1, 2)])
    }

    @Test
    fun dosageGramDoesNotDependOnPanelWidth() {
        val codes = generateDosages(taxa = 5, sites = 17)
        val narrow = dosageGram(codes, 5, 17, panelColumns = 1)
        val wide = dosageGram(codes, 5, 17, panelColumns = 8)
        for (index in narrow.indices) assertEquals(narrow[index], wide[index])
    }

    @Test
    fun publishedChecksums() {
        for (shape in GRAM_SHAPES) {
            val codes = generateDosages(shape.taxa, shape.sites)
            assertEquals(PUBLISHED_CHECKSUMS.getValue(shape), dosageChecksum(codes).toString(), shape.label)
        }
        val tutorial = generateDosages(281, 3_093)
        assertEquals(TUTORIAL_PREVIEW, tutorial.take(16).map { it.toInt() })
    }

    @Test
    fun centeredAgreesWithCenteredIbsOnASmallShape() = runTest {
        assertCenteredWithinTolerance(taxa = 8, sites = 40)
    }

    @Test
    fun centeredAgreesWithCenteredIbsOnTheTutorialShape() = runTest {
        assertCenteredWithinTolerance(taxa = 281, sites = 3_093)
    }

    private suspend fun assertCenteredWithinTolerance(taxa: Int, sites: Int) {
        val codes = generateDosages(taxa, sites)
        val ours = centeredGram(codes, taxa, sites)
        val matrix = centeredIbs(genotypeTable(codes, taxa, sites))
        var maxAbs = 0.0
        var index = 0
        for (first in 0 until taxa) {
            for (second in first until taxa) {
                val scaled = (ours.values[index] / ours.sumPk).toFloat()
                val diff = abs(scaled - matrix.distance(first, second)).toDouble()
                if (diff > maxAbs) maxAbs = diff
                index++
            }
        }
        val sumPk = matrix.annotations.getValue(DistanceMatrix.CENTERED_IBS_SUMPK).toDouble()
        assertTrue(abs(ours.sumPk - sumPk) <= 1.0e-4, "SumPk ${ours.sumPk} vs $sumPk")
        assertTrue(maxAbs <= SCALED_TOLERANCE, "maxAbs=$maxAbs")
    }
}
