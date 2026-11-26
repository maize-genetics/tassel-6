package net.maizegenetics.tassel.distance

import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith

/** [DistanceMatrix.subset]. The value at row `i`, column `j` of the fixture is `10 * i + j` for `j <= i`. */
class DistanceMatrixSubsetTest {
    private val matrix = DistanceMatrix(
        TaxaList.of("A", "B", "C", "D"),
        FloatArray(DistanceMatrix.triangleSize(4)).also { values ->
            var next = 0
            for (row in 0 until 4) for (column in 0..row) values[next++] = (10 * row + column).toFloat()
        },
        mapOf("Matrix_Type" to "Test"),
    )

    @Test
    fun reordersToTheRequestedNames() {
        val subset = matrix.subset(listOf("D", "B", "A")).matrix
        assertEquals(listOf("D", "B", "A"), subset.taxa.names)
        assertEquals(33f, matrix.distance(3, 3))
        assertEquals(33f, subset.distance(0, 0))
        assertEquals(31f, subset.distance(0, 1))
        assertEquals(31f, subset.distance(1, 0))
        assertEquals(30f, subset.distance(2, 0))
        assertEquals(10f, subset.distance(1, 2))
        assertEquals(mapOf("Matrix_Type" to "Test"), subset.annotations)
    }

    @Test
    fun keepsTheFullMatrixInGivenOrderWhenAskedForEveryTaxon() {
        val subset = matrix.subset(matrix.taxa.names).matrix
        assertEquals(matrix.toDelimitedString(), subset.toDelimitedString())
    }

    @Test
    fun partlyOverlappingNamesAreAnErrorByDefault() {
        val error = assertFailsWith<IllegalArgumentException> { matrix.subset(listOf("C", "X", "A", "Y", "X")) }
        assertContains(error.message.orEmpty(), "2 taxa are not in the matrix")
        assertContains(error.message.orEmpty(), "'X', 'Y'")
    }

    @Test
    fun dropReportsWhatWasMissingAndKeepsTheRest() {
        val result = matrix.subset(listOf("C", "X", "A", "Y", "X"), MissingTaxa.Drop)
        assertEquals(listOf("C", "A"), result.matrix.taxa.names)
        assertEquals(listOf("X", "Y"), result.missing)
        assertEquals(20f, result.matrix.distance(0, 1))
        assertEquals(22f, result.matrix.distance(0, 0))
    }

    @Test
    fun duplicatedNamesKeepTheirFirstPositionAndAreReported() {
        val result = matrix.subset(listOf("B", "A", "B", "C", "A", "B"))
        assertEquals(listOf("B", "A", "C"), result.matrix.taxa.names)
        assertEquals(listOf("B", "A"), result.repeated)
        assertEquals(emptyList(), result.missing)
    }

    @Test
    fun namesMatchExactlyIncludingCase() {
        assertFailsWith<IllegalArgumentException> { matrix.subset(listOf("a")) }
    }

    @Test
    fun noNamesGivesAnEmptyMatrix() {
        assertEquals(0, matrix.subset(emptyList()).matrix.numberOfTaxa)
        assertEquals(0, matrix.subset(listOf("X"), MissingTaxa.Drop).matrix.numberOfTaxa)
    }

    @Test
    fun doesNotChangeTheSource() {
        matrix.subset(listOf("D", "A"))
        assertEquals(listOf("A", "B", "C", "D"), matrix.taxa.names)
        assertEquals(31f, matrix.distance(3, 1))
    }
}
