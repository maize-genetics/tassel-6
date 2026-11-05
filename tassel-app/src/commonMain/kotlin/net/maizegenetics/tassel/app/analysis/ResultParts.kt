package net.maizegenetics.tassel.app.analysis

import kotlinx.io.Sink
import net.maizegenetics.tassel.analysis.GenotypeSummary
import net.maizegenetics.tassel.analysis.PrincipalComponentsResult
import net.maizegenetics.tassel.app.export.appendText
import net.maizegenetics.tassel.distance.DistanceMatrix
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.table.PhenotypeTable
import net.maizegenetics.tassel.table.ResultTable
import net.maizegenetics.tassel.table.floatToString

/** TASSEL 5's extension for a table, matrix, or phenotype export, as `tassel-cli` uses it. */
private const val TEXT_EXTENSION = ".txt"

/** TASSEL 5's extension for a HapMap export. */
private const val HAPMAP_EXTENSION = ".hmp.txt"

/**
 * Read-only access to any table the app shows, whatever core type it came from, with cells already spelled the way
 * the exported file spells them. One renderer then draws a genotype summary, a set of eigenvalues, and a kinship
 * matrix without knowing which is which.
 */
interface TableView {
    val columnNames: List<String>

    val rowCount: Int

    /** The cell at [row] and [column], formatted as TASSEL 5 writes it. */
    fun cell(row: Int, column: Int): String
}

/** A picture drawn beside a table. */
sealed interface Chart {
    /** A square symmetric matrix, drawn as a heatmap of [size] by [size] cells. */
    data class Heatmap(val size: Int, val labels: List<String>, val value: (Int, Int) -> Float) : Chart

    /** Two columns of one table, drawn as points. */
    data class Scatter(val xLabel: String, val yLabel: String, val x: FloatArray, val y: FloatArray) : Chart {
        init {
            require(x.size == y.size) { "scatter has ${x.size} x values and ${y.size} y values" }
        }

        override fun equals(other: Any?): Boolean = other is Scatter &&
            xLabel == other.xLabel &&
            yLabel == other.yLabel &&
            x.contentEquals(other.x) &&
            y.contentEquals(other.y)

        override fun hashCode(): Int = 31 * (31 * xLabel.hashCode() + yLabel.hashCode()) + x.contentHashCode()
    }
}

/**
 * One pane of a result: a table to show, an optional picture, and the bytes to export.
 *
 * @property title what the tab is called, which is the core type's own title where it has one.
 * @property fileName the name the save dialog offers, following `tassel-cli`'s extensions.
 * @property table the rows and columns to render.
 * @property chart a picture of the same data, or null.
 * @property write fills a sink with the TASSEL 5 text of this part.
 */
class ResultPart(
    val title: String,
    val fileName: String,
    val table: TableView,
    val chart: Chart? = null,
    val write: (Sink) -> Unit,
)

/**
 * The panes [result] becomes.
 *
 * Like `tassel-cli`'s `outputParts`, this switches on the *output type*, never on the operation, so an operation
 * added to the registry whose output is one of these types needs no change here. [GenotypeTable] is missing on
 * purpose: a filtered genotype table goes back into the data tree as a dataset rather than into a result pane, and
 * [genotypeTablePart] exports it from there.
 *
 * @param baseName file-name stem for the exported files, normally `<dataset>.<operation>`.
 * @throws IllegalArgumentException for a type nothing knows how to show.
 */
fun resultParts(result: Any, baseName: String): List<ResultPart> = when (result) {
    is ResultTable -> listOf(tablePart(result, baseName, null))
    is PhenotypeTable -> listOf(phenotypePart(result, baseName, null))
    is DistanceMatrix -> listOf(matrixPart(result, baseName))
    is GenotypeSummary -> result.tables.map { tablePart(it, baseName, partName(it.title)) }
    is PrincipalComponentsResult -> listOf(
        // TASSEL 5 exports the components themselves as "pcs"; see testdata/golden/MANIFEST [pca].
        phenotypePart(result.principalComponents, baseName, "pcs", scatterOf(result.principalComponents)),
    ) + result.tables.map { tablePart(it, baseName, partName(it.title)) }
    else -> throw IllegalArgumentException(
        "no view for a ${result::class.simpleName} result; add one to resultParts()",
    )
}

/**
 * The file name a genotype dataset is exported under. TASSEL 5 writes genotypes as HapMap rather than as a
 * delimited table, and `HapMap.write` takes the sink straight, which is why this is not a [ResultPart].
 */
fun genotypeFileName(baseName: String): String = baseName.removeSuffix(HAPMAP_EXTENSION) + HAPMAP_EXTENSION

private fun tablePart(table: ResultTable, baseName: String, part: String?) = ResultPart(
    title = table.title,
    fileName = fileName(baseName, part),
    table = ResultTableView(table),
    write = { sink -> sink.appendText { table.writeDelimited(it) } },
)

private fun phenotypePart(table: PhenotypeTable, baseName: String, part: String?, chart: Chart? = null) = ResultPart(
    title = if (part == null) "Phenotype" else "Principal Components",
    fileName = fileName(baseName, part),
    table = PhenotypeTableView(table),
    chart = chart,
    write = { sink -> sink.appendText { table.writeDelimited(it) } },
)

private fun matrixPart(matrix: DistanceMatrix, baseName: String) = ResultPart(
    title = matrix.annotations[DistanceMatrix.MATRIX_TYPE] ?: "Distance Matrix",
    fileName = fileName(baseName, null),
    table = DistanceMatrixView(matrix),
    chart = if (matrix.numberOfTaxa == 0) {
        null
    } else {
        Chart.Heatmap(matrix.numberOfTaxa, matrix.taxa.names) { row, column -> matrix.distance(row, column) }
    },
    write = { sink -> sink.appendText { matrix.writeDelimited(it) } },
)

/** `run.sites.txt`, the way `tassel-cli`'s `outputPath` names one part of a multi-part result. */
private fun fileName(baseName: String, part: String?): String =
    if (part == null) baseName + TEXT_EXTENSION else "$baseName.$part$TEXT_EXTENSION"

/** `Site Summary` becomes `site_summary`, so a title the app never saw still makes a usable file name. */
private fun partName(title: String): String =
    title.lowercase().map { if (it.isLetterOrDigit()) it else '_' }.joinToString("").trim('_')

/** PC1 against PC2, when the components table has both. */
private fun scatterOf(components: PhenotypeTable): Chart.Scatter? {
    val numeric = components.attributes.filter { it.type == PhenotypeTable.AttributeType.covariate }
    if (numeric.size < 2) return null
    return Chart.Scatter(
        xLabel = numeric[0].name,
        yLabel = numeric[1].name,
        x = FloatArray(components.numberOfObservations) { (numeric[0].values[it] as? Float) ?: 0f },
        y = FloatArray(components.numberOfObservations) { (numeric[1].values[it] as? Float) ?: 0f },
    )
}

private class ResultTableView(private val table: ResultTable) : TableView {
    override val columnNames: List<String> get() = table.columnNames
    override val rowCount: Int get() = table.rowCount
    override fun cell(row: Int, column: Int): String = ResultTable.formatCell(table[row, column])
}

private class PhenotypeTableView(private val table: PhenotypeTable) : TableView {
    override val columnNames: List<String> = table.attributes.map { it.name }
    override val rowCount: Int get() = table.numberOfObservations
    override fun cell(row: Int, column: Int): String = PhenotypeTable.formatCell(table.attributes[column].values[row])
}

/** The square matrix as TASSEL 5 writes it: the taxon name, then its whole row. */
private class DistanceMatrixView(private val matrix: DistanceMatrix) : TableView {
    override val columnNames: List<String> = listOf("Taxa") + matrix.taxa.names
    override val rowCount: Int get() = matrix.numberOfTaxa
    override fun cell(row: Int, column: Int): String = if (column == 0) {
        matrix.taxa[row].name
    } else {
        floatToString(matrix.distance(row, column - 1))
    }
}
