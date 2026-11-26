package net.maizegenetics.tassel.dataset

import net.maizegenetics.tassel.distance.DistanceMatrix
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertSame

class DataSetKinshipTest {
    private val kinship = DistanceMatrix(TaxaList.of("A"), floatArrayOf(1f))

    @Test
    fun exposesDistanceMatricesInOrderThroughNestedDataSets() {
        val second = DistanceMatrix(TaxaList.of("B"), floatArrayOf(2f))
        val data = DataSet.of(kinship, "other", DataSet.of(second))
        assertEquals(2, data.distanceMatrices().size)
        assertSame(kinship, data.distanceMatrices()[0])
        assertSame(second, data.distanceMatrices()[1])
    }

    @Test
    fun hasNoDistanceMatricesOrPairingsByDefault() {
        assertEquals(emptyList(), DataSet.of("x").distanceMatrices())
        assertEquals(emptyList(), DataSet.of(kinship).genotypePhenotypes())
    }
}
