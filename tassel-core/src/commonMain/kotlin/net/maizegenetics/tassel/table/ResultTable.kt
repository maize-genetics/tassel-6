package net.maizegenetics.tassel.table

/**
 * A finished table of results: a title, column names, and rows of cells. Analyses return one of
 * these, and the CLI, the app, and the exporters render it without knowing which analysis made it.
 *
 * A cell is a `Double`, an `Int`, a `Long`, a `String`, or `null`. `Double` cells go through
 * [DoubleFormat] when written, so an analysis that holds a count as a `Double` prints it the way
 * TASSEL 5 does. Nothing else is allowed, because a `Float` would silently print its own, different
 * digits.
 *
 * Build one with [resultTable] or [ResultTableBuilder]; the type itself is read only.
 *
 * @see net.maizegenetics.util.TableReport
 * @see net.maizegenetics.util.SimpleTableReport
 */
class ResultTable(val title: String, val columnNames: List<String>, val rows: List<List<Any?>>) {
    init {
        require(columnNames.isNotEmpty()) { "result table '$title' has no columns" }
        rows.forEachIndexed { row, cells ->
            require(cells.size == columnNames.size) {
                "result table '$title' row $row has ${cells.size} cells, expected ${columnNames.size}"
            }
            cells.forEachIndexed { column, cell ->
                require(cell == null || cell is Double || cell is Int || cell is Long || cell is String) {
                    "result table '$title' row $row column '${columnNames[column]}' holds a " +
                        "${cell!!::class.simpleName}; cells must be Double, Int, Long, String, or null"
                }
            }
        }
    }

    val rowCount: Int get() = rows.size

    val columnCount: Int get() = columnNames.size

    /** The cell at [row] and [column], both zero based. */
    operator fun get(row: Int, column: Int): Any? = rows[row][column]

    /** The index of the column named [name], or `-1`. */
    fun columnIndex(name: String): Int = columnNames.indexOf(name)

    /**
     * Appends the table as delimited text, exactly as TASSEL 5's `saveDelimitedTableReport` writes
     * it: one header line, `\n` endings including after the last row, a `null` cell written as
     * nothing, and a `Double` cell written by [DoubleFormat]. The title is not written.
     *
     * @see net.maizegenetics.util.TableReportUtils
     */
    fun writeDelimited(out: Appendable, delimiter: String = "\t", includeHeader: Boolean = true) {
        if (includeHeader) {
            columnNames.forEachIndexed { column, name ->
                if (column != 0) out.append(delimiter)
                out.append(name)
            }
            out.append('\n')
        }
        for (cells in rows) {
            cells.forEachIndexed { column, cell ->
                if (column != 0) out.append(delimiter)
                out.append(formatCell(cell))
            }
            out.append('\n')
        }
    }

    /** The whole table as delimited text; see [writeDelimited]. */
    fun toDelimitedString(delimiter: String = "\t", includeHeader: Boolean = true): String =
        StringBuilder().also { writeDelimited(it, delimiter, includeHeader) }.toString()

    override fun toString(): String = "ResultTable('$title', $rowCount rows, $columnCount columns)"

    companion object {
        /** One cell as TASSEL 5 writes it: nothing for `null`, [DoubleFormat] for a `Double`. */
        fun formatCell(cell: Any?): String = when (cell) {
            null -> ""
            is Double -> DoubleFormat.format(cell)
            else -> cell.toString()
        }
    }
}

/** Collects the rows of a [ResultTable] one at a time. Not thread safe. */
class ResultTableBuilder(val title: String, val columnNames: List<String>) {
    private val rows = mutableListOf<List<Any?>>()

    /** How many rows have been added. */
    val rowCount: Int get() = rows.size

    /** Adds one row; it must have one cell per column name. */
    fun row(cells: List<Any?>) {
        require(cells.size == columnNames.size) {
            "result table '$title' row $rowCount has ${cells.size} cells, expected ${columnNames.size}"
        }
        rows += cells.toList()
    }

    /** Adds one row; it must have one cell per column name. */
    fun row(vararg cells: Any?) = row(cells.asList())

    fun build(): ResultTable = ResultTable(title, columnNames, rows.toList())
}

/**
 * Builds a [ResultTable] row at a time:
 *
 * ```kotlin
 * resultTable("Overview", listOf("Stat Type", "Value")) {
 *     row("Number of Taxa", taxaCount.toDouble())
 * }
 * ```
 */
fun resultTable(title: String, columnNames: List<String>, build: ResultTableBuilder.() -> Unit): ResultTable =
    ResultTableBuilder(title, columnNames).apply(build).build()
