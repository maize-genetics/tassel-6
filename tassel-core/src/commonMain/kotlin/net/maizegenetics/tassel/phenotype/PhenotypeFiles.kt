package net.maizegenetics.tassel.phenotype

import kotlinx.io.Source
import kotlinx.io.readLine
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report

/**
 * Reads TASSEL phenotype files.
 *
 * This reads the `<Phenotype>` format TASSEL 5 writes and
 * [PhenotypeTable.writeDelimited] writes back:
 *
 * ```
 * <Phenotype>
 * taxa     factor  data   covariate
 * Taxon    env     EarHT  Q1
 * B73      IA      64.75  0.014
 * B73      NY      NA     0.014
 * ```
 *
 * The second line types each column by prefix, ignoring case: `tax`, `fac`, `cov`, or `dat`, as
 * TASSEL 5 matches them. The third line names the columns. Then each line is one observation.
 * Cells are separated by runs of whitespace, as TASSEL 5 splits them. Numbers are parsed exactly
 * as Java's `Float.parseFloat` parses them.
 *
 * TASSEL 5's legacy `<Trait>`, `<Covariate>`, `<Factor>`, and `<Header>` files are recognized but
 * not yet read.
 *
 * Where this departs from TASSEL 5 (ADR-0023):
 * - It is more forgiving. Lines are trimmed, a byte-order mark is dropped, and blank lines and
 *   lines starting with `#` are skipped. TASSEL 5 rejects all of these.
 * - Factor cells can be missing too, through the same tokens as numbers.
 * - The taxa column keeps the name the file gives it.
 * - An infinite number is stored as missing, with a warning.
 * - A repeated attribute name, or a type other than the four above, is an error.
 *
 * @see net.maizegenetics.phenotype.PhenotypeBuilder.importPhenotypeFile
 */
object PhenotypeFiles {
    /**
     * Reads one phenotype from [source].
     *
     * @throws IllegalArgumentException if [source] is empty, is not a `<Phenotype>` file, or has a
     *   malformed header, a line of the wrong width, or a cell that is not a number in a numeric
     *   column. The message names the line.
     */
    suspend fun read(
        source: Source,
        options: PhenotypeReadOptions = PhenotypeReadOptions(),
        progress: ProgressReporter = ProgressReporter.None,
    ): PhenotypeRead {
        val lines = Lines(source)
        val tag = lines.next() ?: throw IllegalArgumentException("Phenotype input is empty")
        val format = tag.text.substringBefore(' ').substringBefore('\t').lowercase()
        when {
            format.startsWith("<phenotype") -> Unit
            LEGACY_TAGS.any { format.startsWith(it) } -> throw IllegalArgumentException(
                "Line ${tag.number}: TASSEL 5's legacy <Trait>, <Covariate>, <Factor>, and <Header> phenotype " +
                    "formats are not read yet; only <Phenotype> files are",
            )
            else -> throw IllegalArgumentException(
                "Line ${tag.number}: not a TASSEL phenotype file; it must start with <Phenotype>, " +
                    "found '${excerpt(tag.text)}'",
            )
        }
        return readPhenotype(lines, options, progress)
    }

    /**
     * Whether [token] is one of the missing-value spellings a `<Phenotype>` file is read with by
     * default: `NA` or `NaN` in any case, `.`, or `?`.
     */
    fun isDefaultMissingToken(token: String): Boolean = token.equals("NA", ignoreCase = true) ||
        token.equals("NaN", ignoreCase = true) ||
        token == "." ||
        token == FactorAttribute.MISSING

    private suspend fun readPhenotype(
        lines: Lines,
        options: PhenotypeReadOptions,
        progress: ProgressReporter,
    ): PhenotypeRead {
        val typeLine = lines.next() ?: throw IllegalArgumentException("Phenotype input ends before its types line")
        val nameLine = lines.next() ?: throw IllegalArgumentException("Phenotype input ends before its names line")
        val typeTokens = splitOnWhitespace(typeLine.text)
        val names = splitOnWhitespace(nameLine.text)
        require(typeTokens.size == names.size) {
            "Line ${nameLine.number} names ${names.size} attributes, but line ${typeLine.number} types " +
                "${typeTokens.size}; the counts must match"
        }
        val types = typeTokens.map { token ->
            attributeType(token) ?: throw IllegalArgumentException(
                "Line ${typeLine.number}: unknown attribute type '$token'; expected taxa, factor, covariate, or data",
            )
        }
        val taxaColumns = types.count { it == AttributeType.taxa }
        require(taxaColumns == 1) {
            "Line ${typeLine.number} has $taxaColumns taxa columns; a phenotype needs exactly one"
        }
        names.groupingBy { it }.eachCount().filterValues { it > 1 }.keys.firstOrNull()?.let { repeated ->
            throw IllegalArgumentException("Line ${nameLine.number} names the attribute '$repeated' more than once")
        }

        val isMissing: (String) -> Boolean = options.missingTokens?.let { tokens -> { it in tokens } }
            ?: ::isDefaultMissingToken
        val columns = types.map { type ->
            when (type) {
                AttributeType.taxa -> TaxaColumn()
                AttributeType.factor -> FactorColumn()
                AttributeType.covariate, AttributeType.data -> NumericColumn()
            }
        }
        var infinite = 0
        var observations = 0
        while (true) {
            val line = lines.next() ?: break
            val cells = splitOnWhitespace(line.text)
            require(cells.size == names.size) {
                "Line ${line.number} has ${cells.size} values, but the header names ${names.size} attributes: " +
                    "'${excerpt(line.text)}'"
            }
            cells.forEachIndexed { column, cell ->
                when (val target = columns[column]) {
                    is TaxaColumn -> target.taxa += cell
                    is FactorColumn -> target.labels += if (isMissing(cell)) null else cell
                    is NumericColumn -> {
                        val value = if (isMissing(cell)) {
                            Float.NaN
                        } else {
                            JavaFloat.parse(cell) ?: throw IllegalArgumentException(
                                "Line ${line.number}, attribute '${names[column]}': '$cell' is not a number",
                            )
                        }
                        if (value.isInfinite()) infinite++
                        target.add(value)
                    }
                }
            }
            if (++observations % PROGRESS_OBSERVATIONS == 0) progress.report(observations, 0, PROGRESS_MESSAGE)
        }

        val attributes = columns.mapIndexed { column, values ->
            when (values) {
                is TaxaColumn -> TaxaAttribute(names[column], values.taxa)
                is FactorColumn -> FactorAttribute.of(names[column], values.labels)
                is NumericColumn -> NumericAttribute(names[column], types[column], values.toArray())
            }
        }
        val warnings = if (infinite == 0) {
            emptyList()
        } else {
            listOf("$infinite infinite values were read as missing")
        }
        return PhenotypeRead(PhenotypeTable(attributes, options.name), warnings)
    }

    /** The type a `<Phenotype>` types-line token names, by TASSEL 5's case-insensitive prefixes. */
    private fun attributeType(token: String): AttributeType? {
        val lower = token.lowercase()
        return when {
            lower.startsWith("tax") -> AttributeType.taxa
            lower.startsWith("fac") -> AttributeType.factor
            lower.startsWith("cov") -> AttributeType.covariate
            lower.startsWith("dat") -> AttributeType.data
            else -> null
        }
    }

    /** [line] split on runs of Java `\s` whitespace; [line] is already trimmed. */
    private fun splitOnWhitespace(line: String): List<String> {
        val parts = ArrayList<String>()
        var start = -1
        for (index in line.indices) {
            if (line[index].isFieldSeparator()) {
                if (start >= 0) parts += line.substring(start, index)
                start = -1
            } else if (start < 0) {
                start = index
            }
        }
        if (start >= 0) parts += line.substring(start)
        return parts
    }

    private fun excerpt(line: String): String =
        line.take(EXCERPT_LENGTH) + if (line.length > EXCERPT_LENGTH) "..." else ""

    /** The lines of [source] that carry content, trimmed, with their 1-based line numbers. */
    private class Lines(private val source: Source) {
        private var number = 0

        fun next(): Line? {
            while (true) {
                var text = source.readLine() ?: return null
                number++
                if (number == 1) text = text.removePrefix(BYTE_ORDER_MARK)
                text = text.trim { it.isFieldSeparator() }
                if (text.isNotEmpty() && !text.startsWith("#")) return Line(number, text)
            }
        }
    }

    private class Line(val number: Int, val text: String)

    private sealed interface Column

    private class TaxaColumn : Column {
        val taxa = ArrayList<String>()
    }

    private class FactorColumn : Column {
        val labels = ArrayList<String?>()
    }

    private class NumericColumn : Column {
        private var values = FloatArray(INITIAL_CAPACITY)
        private var size = 0

        fun add(value: Float) {
            if (size == values.size) values = values.copyOf(values.size * 2)
            values[size++] = value
        }

        fun toArray(): FloatArray = values.copyOf(size)
    }

    private val LEGACY_TAGS = listOf("<trait", "<cov", "<fac", "<head")
    private const val BYTE_ORDER_MARK = "\uFEFF"
    private const val EXCERPT_LENGTH = 40
    private const val INITIAL_CAPACITY = 64
    private const val PROGRESS_OBSERVATIONS = 256
    private const val PROGRESS_MESSAGE = "Reading phenotype"
}

/**
 * How [PhenotypeFiles.read] reads a file.
 *
 * @property name the table's display name, which TASSEL 5 takes from the file name.
 * @property missingTokens the exact cell spellings read as missing, in numeric and factor columns.
 *   `null` means the defaults: `NA` and `NaN` in any case, `.`, and `?`.
 */
class PhenotypeReadOptions(val name: String = PhenotypeTable.DEFAULT_NAME, val missingTokens: Set<String>? = null)

/** A phenotype read from a file, and anything worth telling the user about how it was read. */
class PhenotypeRead(val table: PhenotypeTable, val warnings: List<String>)
