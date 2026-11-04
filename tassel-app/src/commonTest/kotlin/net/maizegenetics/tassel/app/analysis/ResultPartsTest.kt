package net.maizegenetics.tassel.app.analysis

import net.maizegenetics.tassel.analysis.genotypeSummary
import net.maizegenetics.tassel.analysis.principalComponents
import net.maizegenetics.tassel.app.testing.testTable
import net.maizegenetics.tassel.app.testing.written
import net.maizegenetics.tassel.distance.DistanceMatrix
import net.maizegenetics.tassel.table.PhenotypeTable
import net.maizegenetics.tassel.table.resultTable
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertIs
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * Turning an operation's output into panes the UI can show and export.
 *
 * Every assertion about exported text compares against the core writer's own `toDelimitedString`, so the app can
 * never drift from the TASSEL 5 format by reformatting something on the way out.
 */
class ResultPartsTest {
    private val matrix = DistanceMatrix(
        taxa = TaxaList.of(listOf("B73", "Mo17", "CML247")),
        lowerTriangle = floatArrayOf(2f, 0.5f, 2f, -0.25f, 0.75f, 2f),
        annotations = mapOf(DistanceMatrix.MATRIX_TYPE to "Centered_IBS"),
    )

    @Test
    fun aResultTableIsOnePaneAndExportsExactlyWhatTheCoreWriterWrites() {
        val table = resultTable("Overall Summary", listOf("Stat Type", "Value")) {
            row("Number of Taxa", 3.0)
        }

        val parts = resultParts(table, "mdp.summary")

        assertEquals(1, parts.size)
        assertEquals("Overall Summary", parts[0].title)
        assertEquals("mdp.summary.txt", parts[0].fileName)
        assertNull(parts[0].chart)
        assertEquals(table.toDelimitedString(), written(parts[0].write))
    }

    @Test
    fun aGenotypeSummaryBecomesOnePanePerTableWithFileNamesFromTheirTitles() {
        val summary = genotypeSummary(testTable(numberOfTaxa = 3, numberOfSites = 4))

        val parts = resultParts(summary, "mdp.summary")

        assertEquals(summary.tables.map { it.title }, parts.map { it.title })
        assertEquals(
            listOf(
                "mdp.summary.overall_summary.txt",
                "mdp.summary.allele_summary.txt",
                "mdp.summary.site_summary.txt",
                "mdp.summary.taxa_summary.txt",
            ),
            parts.map { it.fileName },
        )
        summary.tables.forEachIndexed { index, table ->
            assertEquals(table.toDelimitedString(), written(parts[index].write), table.title)
        }
    }

    @Test
    fun aDistanceMatrixGetsAHeatmapAndTheSquareMatrixAsItsTable() {
        val parts = resultParts(matrix, "mdp.kinship")

        assertEquals(1, parts.size)
        assertEquals("Centered_IBS", parts[0].title)
        assertEquals("mdp.kinship.txt", parts[0].fileName)
        assertEquals(matrix.toDelimitedString(), written(parts[0].write))

        val heatmap = assertIs<Chart.Heatmap>(parts[0].chart)
        assertEquals(3, heatmap.size)
        assertEquals(listOf("B73", "Mo17", "CML247"), heatmap.labels)
        // The chart reads the matrix itself, in either order, so nothing is copied for it.
        assertEquals(0.5f, heatmap.value(1, 0))
        assertEquals(0.5f, heatmap.value(0, 1))
    }

    @Test
    fun theMatrixTableIsTheNameThenTheWholeRow() {
        val table = resultParts(matrix, "mdp.kinship")[0].table

        assertEquals(listOf("Taxa", "B73", "Mo17", "CML247"), table.columnNames)
        assertEquals(3, table.rowCount)
        assertEquals("Mo17", table.cell(1, 0))
        assertEquals("0.5", table.cell(1, 1))
        assertEquals("-0.25", table.cell(2, 1))
    }

    @Test
    fun principalComponentsGiveTheComponentsFirstWithAScatterOfTheFirstTwo() {
        val result = principalComponents(testTable(numberOfTaxa = 4, numberOfSites = 6), numberOfComponents = 3)

        val parts = resultParts(result, "mdp.pca")

        // TASSEL 5 exports the components as "pcs", then the eigenvalues, then the eigenvectors.
        assertEquals(
            listOf("mdp.pca.pcs.txt", "mdp.pca.proportion_of_variance_explained.txt", "mdp.pca.eigenvectors.txt"),
            parts.map { it.fileName },
        )
        assertEquals(result.principalComponents.toDelimitedString(), written(parts[0].write))

        val scatter = assertIs<Chart.Scatter>(parts[0].chart)
        assertEquals("PC1", scatter.xLabel)
        assertEquals("PC2", scatter.yLabel)
        assertEquals(4, scatter.x.size)
        assertEquals(4, scatter.y.size)
    }

    @Test
    fun theComponentsTableShowsTheTaxonNameAndTheFloatsTasselWrites() {
        val result = principalComponents(testTable(numberOfTaxa = 4, numberOfSites = 6), numberOfComponents = 2)

        val table = resultParts(result, "mdp.pca")[0].table

        assertEquals(listOf("Taxa", "PC1", "PC2"), table.columnNames)
        assertEquals("taxon0", table.cell(0, 0))
        assertEquals(
            PhenotypeTable.formatCell(result.principalComponents.attributes[1].values[0]),
            table.cell(0, 1),
        )
    }

    @Test
    fun aPhenotypeWithOneCovariateHasNothingToPlotAgainst() {
        val phenotype = PhenotypeTable(
            listOf(PhenotypeTable.taxa(listOf("B73", "Mo17")), PhenotypeTable.Attribute("PC1", floatArrayOf(1f, 2f))),
        )

        val parts = resultParts(phenotype, "mdp.pca")

        assertEquals(1, parts.size)
        assertNull(parts[0].chart)
        assertEquals(phenotype.toDelimitedString(), written(parts[0].write))
    }

    @Test
    fun anOutputTypeNothingKnowsHowToShowSaysSoRatherThanDrawingNothing() {
        val failure = assertFailsWith<IllegalArgumentException> { resultParts("a bare string", "mdp.mystery") }

        assertTrue(failure.message!!.contains("resultParts"), failure.message!!)
    }

    @Test
    fun aGenotypeExportKeepsTheHapMapExtensionAndDoesNotDoubleIt() {
        assertEquals("mdp_genotype.hmp.txt", genotypeFileName("mdp_genotype.hmp.txt"))
        assertEquals("filtered.hmp.txt", genotypeFileName("filtered"))
    }
}
