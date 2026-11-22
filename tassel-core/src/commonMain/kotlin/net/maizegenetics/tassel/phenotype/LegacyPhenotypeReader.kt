package net.maizegenetics.tassel.phenotype

import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report

/**
 * Reads TASSEL 5's legacy phenotype files: `<Trait>`, `<Covariate>`, `<Factor>`, and `<Header>`.
 *
 * Tag lines come first, then one taxon per line, as TASSEL 5 reads them:
 *
 * ```
 * <Covariate>                     optional: every trait is a covariate rather than data
 * <Trait>   Q1     Q2     Q3      names the columns after the taxon
 * 33-16     0.014  0.972  0.014
 * ```
 *
 * - `<Trait>` makes each column numeric: data, or covariates after a `<Covariate>` line.
 * - `<Factor>` after `<Trait>` makes each column a factor.
 * - `<Header name=x>` lines, one per factor `x`, give each column a level of that factor. The file
 *   is reshaped to long format: the same trait can head several columns, one per level
 *   combination, and each data line becomes one observation per distinct combination, in
 *   first-appearance order. The traits come out in first-appearance order too, after the taxa and
 *   factor columns. A (trait, combination) pair the file does not give is missing.
 *
 * Where this departs from TASSEL 5 (ADR-0023):
 * - An unparseable cell is an error naming its line and column. TASSEL 5 treats it as missing.
 * - `-999` and `-99` are missing in `<Header>` files too; TASSEL 5 reads them there as numbers.
 * - The default missing spellings also apply to factor cells.
 * - A repeated (trait, level combination), a header of the wrong width, a tag after the first
 *   data line, a `<Covariate>` with `<Factor>`, and a `<Factor>` with `<Header>` are errors.
 *   TASSEL 5 overwrites, or fails with an index error, or ignores the tag.
 * - `<Header>` rows and columns keep first-appearance order; TASSEL 5 sorts both.
 *
 * @see net.maizegenetics.phenotype.PhenotypeBuilder.importTraits
 */
internal object LegacyPhenotypeReader {
    /** Whether [line], a trimmed line, is one of the four legacy tags, by TASSEL 5's prefixes. */
    fun recognizes(line: String): Boolean = tagOf(line) != null

    suspend fun read(
        first: PhenotypeFiles.Line,
        lines: PhenotypeFiles.Lines,
        options: PhenotypeReadOptions,
        progress: ProgressReporter,
    ): PhenotypeRead {
        val layout = readTags(first, lines)
        val isMissing: (String) -> Boolean = options.missingTokens?.let { tokens -> { it in tokens } }
            ?: PhenotypeFiles::isLegacyDefaultMissingToken
        val reader = Cells(layout, isMissing)
        val target = if (layout.headers.isEmpty()) Wide(layout, reader) else Reshaped(layout, reader)
        var observations = 0
        var line = layout.firstData
        while (line != null) {
            if (line.text.startsWith("<")) {
                throw IllegalArgumentException(
                    "Line ${line.number}: the tag '${PhenotypeFiles.excerpt(line.text)}' comes after data; " +
                        "every tag must come before the first data line",
                )
            }
            val cells = PhenotypeFiles.splitOnWhitespace(line.text)
            require(cells.size == layout.traits.size + 1) {
                "Line ${line.number} has ${cells.size} values, but the <Trait> line names ${layout.traits.size} " +
                    "traits, so ${layout.traits.size + 1} are expected: '${PhenotypeFiles.excerpt(line.text)}'"
            }
            target.add(line.number, cells)
            if (++observations % PhenotypeFiles.PROGRESS_OBSERVATIONS == 0) {
                progress.report(observations, 0, PhenotypeFiles.PROGRESS_MESSAGE)
            }
            line = lines.next()
        }
        val warnings = if (reader.infinite == 0) {
            emptyList()
        } else {
            listOf("${reader.infinite} infinite values were read as missing")
        }
        return PhenotypeRead(PhenotypeTable(target.attributes(), options.name), warnings)
    }

    private enum class Tag(val display: String) {
        Trait("<Trait>"),
        Covariate("<Covariate>"),
        Factor("<Factor>"),
        Header("<Header name=x>"),
    }

    /** What the tag lines say, and the first line after them. */
    private class Layout(
        val traits: List<String>,
        val type: AttributeType,
        val headers: List<HeaderLine>,
        val firstData: PhenotypeFiles.Line?,
    )

    private class HeaderLine(val number: Int, val factor: String, val levels: List<String>)

    private fun tagOf(line: String): Tag? {
        if (!line.startsWith("<")) return null
        val word = line.drop(1).takeWhile { it != '>' && it != '=' && !it.isFieldSeparator() }.lowercase()
        return when {
            word.startsWith("trai") -> Tag.Trait
            word.startsWith("cov") -> Tag.Covariate
            word.startsWith("fac") -> Tag.Factor
            word.startsWith("head") -> Tag.Header
            else -> null
        }
    }

    private fun readTags(first: PhenotypeFiles.Line, lines: PhenotypeFiles.Lines): Layout {
        var traits: List<String>? = null
        var traitLine = 0
        var covariate = false
        var factor = false
        val headers = ArrayList<HeaderLine>()
        var line: PhenotypeFiles.Line? = first
        while (line != null && line.text.startsWith("<")) {
            when (tagOf(line.text)) {
                Tag.Trait -> {
                    require(traits == null) { "Line ${line.number}: a second <Trait> line; line $traitLine has one" }
                    traitLine = line.number
                    traits = line.text.split(TRAIT_SEPARATORS).filter { it.isNotEmpty() }.drop(1)
                }
                Tag.Covariate -> covariate = true
                Tag.Factor -> factor = true
                Tag.Header -> headers += parseHeader(line)
                null -> throw IllegalArgumentException(
                    "Line ${line.number}: unknown tag '${PhenotypeFiles.excerpt(line.text)}'; the accepted tags are " +
                        "<Phenotype>, <Trait>, <Covariate>, <Factor>, and <Header>",
                )
            }
            line = lines.next()
        }
        val names = traits ?: throw IllegalArgumentException(
            "Line ${first.number}: ${tagOf(first.text)?.display} needs a <Trait> line that names the columns",
        )
        require(names.isNotEmpty()) { "Line $traitLine: <Trait> names no traits" }
        require(!(factor && covariate)) { "<Covariate> and <Factor> both appear; a file's columns have one type" }
        require(!(factor && headers.isNotEmpty())) {
            "<Factor> and <Header> both appear; a <Header> file's columns are numeric"
        }
        headers.groupingBy { it.factor }.eachCount().filterValues { it > 1 }.keys.firstOrNull()?.let { repeated ->
            throw IllegalArgumentException("<Header name=$repeated> appears more than once")
        }
        headers.forEach { header ->
            require(header.levels.size == names.size) {
                "Line ${header.number}: <Header name=${header.factor}> gives ${header.levels.size} levels, but the " +
                    "<Trait> line names ${names.size} traits"
            }
        }
        val type = when {
            factor -> AttributeType.factor
            covariate -> AttributeType.covariate
            else -> AttributeType.data
        }
        return Layout(names, type, headers, line)
    }

    /** `<Header name=x> a b c`, split as TASSEL 5 splits it, on `<`, `>`, `=`, and whitespace. */
    private fun parseHeader(line: PhenotypeFiles.Line): HeaderLine {
        val tokens = line.text.split(HEADER_SEPARATORS).filter { it.isNotEmpty() }
        require(tokens.size >= 3 && tokens[1].equals("name", ignoreCase = true)) {
            "Line ${line.number}: malformed header '${PhenotypeFiles.excerpt(line.text)}'; expected " +
                "<Header name=factor> followed by one level per trait"
        }
        return HeaderLine(line.number, tokens[2], tokens.drop(3))
    }

    private val TRAIT_SEPARATORS = Regex("[<>\\s]+")
    private val HEADER_SEPARATORS = Regex("[<>=\\s]+")
    private const val TAXA_NAME = "Taxa"

    /** Parses cells, and counts the infinite numbers it turns into missing ones. */
    private class Cells(private val layout: Layout, private val isMissing: (String) -> Boolean) {
        var infinite = 0
            private set

        fun number(cell: String, line: Int, trait: Int): Float {
            if (isMissing(cell)) return Float.NaN
            val value = JavaFloat.parse(cell) ?: throw IllegalArgumentException(
                "Line $line, column ${trait + 2}, trait '${layout.traits[trait]}': '$cell' is not a number",
            )
            if (value.isInfinite()) infinite++
            return value
        }

        fun level(cell: String): String? = if (isMissing(cell)) null else cell
    }

    private sealed interface Target {
        fun add(line: Int, cells: List<String>)

        fun attributes(): List<PhenotypeAttribute>
    }

    /** `<Trait>` and `<Factor>`: one observation per line, one column per trait. */
    private class Wide(private val layout: Layout, private val parser: Cells) : Target {
        private val taxa = ArrayList<String>()
        private val numbers = if (layout.type == AttributeType.factor) {
            emptyList()
        } else {
            layout.traits.map { PhenotypeFiles.NumericColumn() }
        }
        private val factors = if (layout.type == AttributeType.factor) {
            layout.traits.map { ArrayList<String?>() }
        } else {
            emptyList()
        }

        override fun add(line: Int, cells: List<String>) {
            taxa += cells[0]
            layout.traits.indices.forEach { trait ->
                val cell = cells[trait + 1]
                if (factors.isEmpty()) {
                    numbers[trait].add(parser.number(cell, line, trait))
                } else {
                    factors[trait] += parser.level(cell)
                }
            }
        }

        override fun attributes(): List<PhenotypeAttribute> = listOf<PhenotypeAttribute>(
            TaxaAttribute(TAXA_NAME, taxa),
        ) + layout.traits.mapIndexed { trait, name ->
            if (factors.isEmpty()) {
                NumericAttribute(name, layout.type, numbers[trait].toArray())
            } else {
                FactorAttribute.of(name, factors[trait])
            }
        }
    }

    /**
     * `<Header>`: each trait name heads columns for several level combinations. Every data line
     * becomes one observation per combination, and a trait with no column for a combination is
     * missing there.
     */
    private class Reshaped(private val layout: Layout, private val parser: Cells) : Target {
        /** Each distinct level combination, in first-appearance order, one level per header line. */
        private val combinations = LinkedHashMap<List<String>, Int>()
        private val traitNames = LinkedHashMap<String, Int>()

        /** For each file column, its trait's and its combination's index. */
        private val traitOf: IntArray
        private val combinationOf: IntArray

        private val taxa = ArrayList<String>()
        private val factors: List<ArrayList<String>>
        private val numbers: List<PhenotypeFiles.NumericColumn>

        init {
            val seen = HashSet<Pair<Int, Int>>()
            traitOf = IntArray(layout.traits.size)
            combinationOf = IntArray(layout.traits.size)
            layout.traits.forEachIndexed { column, name ->
                val levels = layout.headers.map { it.levels[column] }
                traitOf[column] = traitNames.getOrPut(name) { traitNames.size }
                combinationOf[column] = combinations.getOrPut(levels) { combinations.size }
                require(seen.add(traitOf[column] to combinationOf[column])) {
                    "Column ${column + 2} repeats trait '$name' for " +
                        layout.headers.indices.joinToString(", ") { "${layout.headers[it].factor}=${levels[it]}" }
                }
            }
            factors = layout.headers.map { ArrayList() }
            numbers = traitNames.keys.map { PhenotypeFiles.NumericColumn() }
        }

        override fun add(line: Int, cells: List<String>) {
            val block = FloatArray(traitNames.size * combinations.size) { Float.NaN }
            layout.traits.indices.forEach { column ->
                block[traitOf[column] * combinations.size + combinationOf[column]] =
                    parser.number(cells[column + 1], line, column)
            }
            combinations.keys.forEachIndexed { combination, levels ->
                taxa += cells[0]
                levels.forEachIndexed { header, level -> factors[header] += level }
                numbers.forEachIndexed { trait, values -> values.add(block[trait * combinations.size + combination]) }
            }
        }

        override fun attributes(): List<PhenotypeAttribute> = listOf<PhenotypeAttribute>(
            TaxaAttribute(TAXA_NAME, taxa),
        ) + layout.headers.mapIndexed { header, line ->
            FactorAttribute.of(line.factor, factors[header])
        } + traitNames.keys.mapIndexed { trait, name ->
            NumericAttribute(name, layout.type, numbers[trait].toArray())
        }
    }
}
