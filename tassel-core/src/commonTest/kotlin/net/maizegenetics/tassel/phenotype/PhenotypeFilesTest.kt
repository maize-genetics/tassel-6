package net.maizegenetics.tassel.phenotype

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.writeString
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.Tolerance
import net.maizegenetics.tassel.testing.assertMatchesGolden
import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.goldenTolerance
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * Reading `<Phenotype>` files: TASSEL 5's own test inputs against what TASSEL 5 exports from them
 * (`testdata/golden/phenotype-format/`, ADR-0009), and the reading rules of ADR-0023 on small
 * hand-written files.
 */
class PhenotypeFilesTest {
    @Test
    fun theGoldensAreComparedAsPhenotypes() {
        GOLDENS.keys.forEach { golden -> assertEquals(Tolerance.Phenotype(), goldenTolerance(golden)) }
    }

    @Test
    fun readsTassel5sTestInputsToTheValuesTassel5Exports() = runTest {
        GOLDENS.forEach { (golden, input) ->
            assertMatchesGolden(golden, readFile(input).table.toDelimitedString())
        }
    }

    /** Writing a golden back out after reading it gives the golden's exact bytes, floats and all. */
    @Test
    fun writesEachGoldenBackByteForByte() = runTest {
        GOLDENS.keys.forEach { golden ->
            val text = goldenText(golden)
            assertEquals(normalizeLineEndings(text), read(text).table.toDelimitedString(), golden)
        }
    }

    @Test
    fun readsTheBasicTestFileIntoTypedAttributes() = runTest {
        val table = readFile("phenotype/basic_phenotype_test.txt").table
        assertEquals(listOf("Taxon", "Factor1", "Trait1", "Trait2", "Trait3"), table.attributes.map { it.name })
        assertEquals(
            listOf(
                AttributeType.taxa,
                AttributeType.factor,
                AttributeType.data,
                AttributeType.data,
                AttributeType.covariate,
            ),
            table.attributes.map { it.type },
        )
        assertEquals(7, table.numberOfObservations)
        assertEquals(listOf("t1", "t2", "t3", "t4"), table.taxa.names)
        assertTrue(table.hasReplicatedTaxa)
        assertEquals(listOf("A", "B"), (table.attribute("Factor1") as FactorAttribute).levels)
        val trait1 = table.attribute("Trait1") as NumericAttribute
        assertTrue(trait1.isMissing(3))
        assertEquals(1.1f, trait1.floatValue(0))
        assertEquals(1, trait1.missingCount)
    }

    /** `earht_with_rep.txt` starts with a lowercase `<phenotype>` tag; TASSEL 5 accepts that too. */
    @Test
    fun acceptsTheTagInAnyCase() = runTest {
        val table = readFile("phenotype/earht_with_rep.txt").table
        assertEquals(listOf("taxon", "rep", "earht"), table.attributes.map { it.name })
        assertEquals(listOf("1", "2"), (table.attribute("rep") as FactorAttribute).levels)
        assertTrue(table.hasReplicatedTaxa)
    }

    @Test
    fun readsWhatItWritesBack() = runTest {
        val table = PhenotypeTable(
            listOf(
                TaxaAttribute("Taxon", listOf("B73", "B73", "Mo17", "W22")),
                FactorAttribute.of("env", listOf("IA", "NY", null, "IA")),
                FactorAttribute.of("rep", listOf("10", "2", "2", null)),
                NumericAttribute("EarHT", AttributeType.data, floatArrayOf(64.75f, Float.NaN, 1.0e-7f, -3.4028235e38f)),
                NumericAttribute("Q1", AttributeType.covariate, floatArrayOf(0.014f, 0.014f, Float.NaN, 1.4e-45f)),
            ),
        )
        assertEquals(table, read(table.toDelimitedString()).table)
    }

    @Test
    fun forgivesWhatTassel5Rejects() = runTest {
        val text = "\uFEFF# a comment before the tag\r\n" +
            "  <Phenotype>  \r\n" +
            "\r\n" +
            "taxa\t data \tfactor\r\n" +
            "# a comment between header lines\r\n" +
            "Taxa  EarHT   env\r\n" +
            "   B73   64.75   IA\r\n" +
            "\r\n" +
            "Mo17\t.\tNA\r\n"
        val table = read(text).table
        assertEquals(listOf("Taxa", "EarHT", "env"), table.attributes.map { it.name })
        assertEquals(listOf("B73", "Mo17"), table.taxa.names)
        assertTrue(table.attribute("EarHT").isMissing(1))
        assertTrue(table.attribute("env").isMissing(1))
    }

    @Test
    fun readsTheDefaultMissingSpellingsInEveryCase() = runTest {
        val table = read(header("data", "factor") + rows("NA ?", "na NA", "NaN nan", "NAN .", "? x")).table
        assertEquals(5, table.attribute("value").missingCount)
        assertEquals(4, table.attribute("level").missingCount)
    }

    /** -999 is a number in a `<Phenotype>` file, as it is in TASSEL 5; only `<Trait>` files treat it as missing. */
    @Test
    fun readsMinus999AsANumber() = runTest {
        val value = read(header("data", "factor") + rows("-999 a")).table.attribute("value") as NumericAttribute
        assertEquals(-999f, value.floatValue(0))
    }

    @Test
    fun replacesTheMissingSpellingsWhenTold() = runTest {
        val options = PhenotypeReadOptions(missingTokens = setOf("-9", "none"))
        val table = read(header("data", "factor") + rows("-9 none", "1.5 NA"), options).table
        val value = table.attribute("value") as NumericAttribute
        val level = table.attribute("level") as FactorAttribute
        assertTrue(value.isMissing(0))
        assertEquals(1.5f, value.floatValue(1))
        assertNull(level.level(0))
        assertEquals("NA", level.level(1))
        val failure =
            assertFailsWith<IllegalArgumentException> { read(header("data", "factor") + rows("NA a"), options) }
        assertEquals("Line 4, attribute 'value': 'NA' is not a number", failure.message)
    }

    @Test
    fun readsJavaFloatSpellingsAndWarnsAboutInfinity() = runTest {
        val read = read(header("data", "factor") + rows("1.5f a", "2D b", "Infinity c", "1e39 d"))
        val value = read.table.attribute("value") as NumericAttribute
        assertEquals(listOf(1.5f, 2f), listOf(value.floatValue(0), value.floatValue(1)))
        assertTrue(value.isMissing(2) && value.isMissing(3))
        assertEquals(listOf("2 infinite values were read as missing"), read.warnings)
    }

    @Test
    fun typesColumnsByTassel5sPrefixesInAnyCase() = runTest {
        val text = "<Phenotype>\nTAXON\tFactors\tcov\tDataPoint\nT\tF\tC\tD\nb73\tx\t1\t2\n"
        assertEquals(
            listOf(AttributeType.taxa, AttributeType.factor, AttributeType.covariate, AttributeType.data),
            read(text).table.attributes.map { it.type },
        )
    }

    @Test
    fun namesTheTableAsTold() = runTest {
        assertEquals(
            "mdp_traits",
            read(header("data", "factor") + rows("1 a"), PhenotypeReadOptions(name = "mdp_traits")).table.name,
        )
        assertEquals(PhenotypeTable.DEFAULT_NAME, read(header("data", "factor") + rows("1 a")).table.name)
    }

    @Test
    fun readsAHeaderWithNoObservations() = runTest {
        val table = read(header("data", "factor")).table
        assertEquals(0, table.numberOfObservations)
        assertEquals(3, table.numberOfAttributes)
    }

    @Test
    fun reportsProgressWhileReading() = runTest {
        val reports = mutableListOf<Progress>()
        val text = header("data", "factor") + rows(*Array(600) { "$it a" })
        PhenotypeFiles.read(source(text), progress = { reports += it })
        assertEquals(listOf(256L, 512L), reports.map { it.completed })
    }

    @Test
    fun explainsEveryMalformedFile() = runTest {
        val cases = mapOf(
            "" to "Phenotype input is empty",
            "# only a comment\n\n" to "Phenotype input is empty",
            "<Phenotype>\n" to "Phenotype input ends before its types line",
            "<Phenotype>\ntaxa data\n" to "Phenotype input ends before its names line",
            "Taxa\tEarHT\nB73\t1\n" to
                "Line 1: not a TASSEL phenotype file; it must start with <Phenotype>, found 'Taxa\tEarHT'",
            "<Phenotype>\ntaxa data\nTaxa\n" to "Line 3 names 1 attributes, but line 2 types 2; the counts must match",
            "<Phenotype>\ntaxa trait\nTaxa EarHT\n" to
                "Line 2: unknown attribute type 'trait'; expected taxa, factor, covariate, or data",
            "<Phenotype>\ndata data\nA B\n" to "Line 2 has 0 taxa columns; a phenotype needs exactly one",
            "<Phenotype>\ntaxa taxa\nA B\n" to "Line 2 has 2 taxa columns; a phenotype needs exactly one",
            "<Phenotype>\ntaxa data data\nTaxa x x\n" to "Line 3 names the attribute 'x' more than once",
            "<Phenotype>\ntaxa data\nTaxa EarHT\nB73 1 2\n" to
                "Line 4 has 3 values, but the header names 2 attributes: 'B73 1 2'",
            "<Phenotype>\ntaxa data\nTaxa EarHT\nB73 tall\n" to "Line 4, attribute 'EarHT': 'tall' is not a number",
            "<Phenotype>\ntaxa data\nTaxa EarHT\nB73 0x1p3\n" to "Line 4, attribute 'EarHT': '0x1p3' is not a number",
        )
        cases.forEach { (text, expected) ->
            val failure = assertFailsWith<IllegalArgumentException>(text) { read(text) }
            assertEquals(expected, failure.message, text)
        }
    }

    /** The tutorial traits are a `<Trait>` file, which issue #157 teaches the reader. */
    @Test
    fun saysLegacyFormatsAreNotReadYet() = runTest {
        listOf("tutorial/mdp_traits.txt", "tutorial/mdp_population_structure.txt").forEach { input ->
            val failure = assertFailsWith<IllegalArgumentException> { readFile(input) }
            assertEquals(
                "Line 1: TASSEL 5's legacy <Trait>, <Covariate>, <Factor>, and <Header> phenotype formats are " +
                    "not read yet; only <Phenotype> files are",
                failure.message,
            )
        }
    }

    private companion object {
        /** Each golden and the input TASSEL 5 exported it from. */
        val GOLDENS = mapOf(
            "golden/phenotype-format/basic_phenotype_test.txt" to "phenotype/basic_phenotype_test.txt",
            "golden/phenotype-format/earht_with_rep.txt" to "phenotype/earht_with_rep.txt",
        )

        fun header(valueType: String, levelType: String) =
            "<Phenotype>\ntaxa\t$valueType\t$levelType\nTaxa\tvalue\tlevel\n"

        /** One `t<index>` observation per entry, each entry holding the other cells. */
        fun rows(vararg cells: String) = cells.mapIndexed { index, rest -> "t$index $rest\n" }.joinToString("")

        fun source(text: String) = Buffer().apply { writeString(text) }

        suspend fun read(text: String, options: PhenotypeReadOptions = PhenotypeReadOptions()) =
            PhenotypeFiles.read(source(text), options)

        suspend fun readFile(path: String) = PhenotypeFiles.read(Buffer().apply { write(TestData.readBytes(path)) })
    }
}
