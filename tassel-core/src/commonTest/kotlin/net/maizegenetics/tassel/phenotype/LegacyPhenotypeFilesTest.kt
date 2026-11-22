package net.maizegenetics.tassel.phenotype

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.writeString
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.Tolerance
import net.maizegenetics.tassel.testing.assertMatchesGolden
import net.maizegenetics.tassel.testing.goldenTolerance
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * Reading TASSEL 5's legacy `<Trait>`, `<Covariate>`, `<Factor>`, and `<Header>` files: TASSEL 5's
 * own inputs against what TASSEL 5 exports from them (`testdata/golden/phenotype-legacy/`), and
 * the reading rules of ADR-0023 on small hand-written files.
 */
class LegacyPhenotypeFilesTest {
    @Test
    fun theGoldensAreComparedAsPhenotypesWithRowsUnordered() {
        GOLDENS.keys.forEach { golden ->
            assertEquals(Tolerance.Phenotype(rowsUnordered = true), goldenTolerance(golden))
        }
    }

    @Test
    fun readsTassel5sInputsToTheValuesTassel5Exports() = runTest {
        GOLDENS.forEach { (golden, input) -> assertMatchesGolden(golden, readFile(input).table.toDelimitedString()) }
    }

    /** `mdp_traits.txt` has CRLF line endings, a trailing tab on the `<Trait>` line, and `-999`. */
    @Test
    fun readsTheTutorialTraits() = runTest {
        val table = readFile("tutorial/mdp_traits.txt").table
        assertEquals(listOf("Taxa", "EarHT", "dpoll", "EarDia"), table.attributes.map { it.name })
        assertEquals(
            listOf(AttributeType.taxa, AttributeType.data, AttributeType.data, AttributeType.data),
            table.attributes.map { it.type },
        )
        assertEquals(listOf("811", "33-16", "38-11"), table.taxa.names.take(3))
        assertEquals(59.5f, (table.attribute("EarHT") as NumericAttribute).floatValue(0))
        assertTrue(table.attribute("dpoll").isMissing(0))
        assertTrue(table.attribute("EarDia").isMissing(1))
        assertEquals(37.897f, (table.attribute("EarDia") as NumericAttribute).floatValue(2))
    }

    @Test
    fun aCovariateLineMakesEveryColumnACovariate() = runTest {
        val table = readFile("tutorial/mdp_population_structure.txt").table
        assertEquals(listOf("Taxa", "Q1", "Q2", "Q3"), table.attributes.map { it.name })
        assertEquals(
            List(3) { AttributeType.covariate },
            table.attributes.drop(1).map { it.type },
        )
    }

    @Test
    fun readsMinus999AndMinus99AsMissingButMinus999Point0AsANumber() = runTest {
        val table = read("<Trait>\tA\tB\nt1\t-999\t-99\nt2\t-999.0\t-99.0\nt3\t-9990\t-100\n").table
        val a = table.attribute("A") as NumericAttribute
        val b = table.attribute("B") as NumericAttribute
        assertEquals(listOf(true, false, false), List(3) { a.isMissing(it) })
        assertEquals(listOf(true, false, false), List(3) { b.isMissing(it) })
        assertEquals(-999f, a.floatValue(1))
        assertEquals(-99f, b.floatValue(1))
    }

    @Test
    fun readsTheOtherMissingSpellingsInAnyCase() = runTest {
        val table = read("<Trait>\tA\nt1\tNA\nt2\tnan\nt3\t.\nt4\t?\nt5\tNaN\nt6\t1\n").table
        assertEquals(5, table.attribute("A").missingCount)
    }

    @Test
    fun anUnparseableCellIsAnErrorNamingItsLineAndColumn() = runTest {
        val failure = assertFailsWith<IllegalArgumentException> { read("<Trait>\tA\tB\nt1\t1\t2\nt2\t3\ttall\n") }
        assertEquals("Line 3, column 3, trait 'B': 'tall' is not a number", failure.message)
        assertFailsWith<IllegalArgumentException> { read("<Trait>\tA\nt1\t0x1p3\n") }
    }

    @Test
    fun replacesTheMissingSpellingsWhenTold() = runTest {
        val options = PhenotypeReadOptions(missingTokens = setOf("-1", "none"))
        val table = read("<Trait>\tA\nt1\t-1\nt2\t-999\nt3\tnone\n", options).table
        val a = table.attribute("A") as NumericAttribute
        assertEquals(listOf(true, false, true), List(3) { a.isMissing(it) })
        assertEquals(-999f, a.floatValue(1))
        assertFailsWith<IllegalArgumentException> { read("<Trait>\tA\nt1\tNA\n", options) }
    }

    @Test
    fun aFactorLineMakesEveryColumnAFactor() = runTest {
        val table = read("<Trait>\tenv\trep\n<Factor>\nt1\tIA\t10\nt2\tNY\t2\nt3\t?\t10\n").table
        assertEquals(List(2) { AttributeType.factor }, table.attributes.drop(1).map { it.type })
        assertEquals(listOf("IA", "NY"), (table.attribute("env") as FactorAttribute).levels)
        assertEquals(listOf("10", "2"), (table.attribute("rep") as FactorAttribute).levels)
        assertTrue(table.attribute("env").isMissing(2))
    }

    @Test
    fun theTagsMayComeInAnyOrderAndCase() = runTest {
        val table = read("# comment\n<trait>\tQ1\n<COVARIATE>\nt1\t0.5\n").table
        assertEquals(AttributeType.covariate, table.attribute("Q1").type)
    }

    @Test
    fun detectsTheFormatAfterABomBlankAndCommentLines() = runTest {
        val table = read("\r\n# note\r\n\r\n  <Trait>  Q1  \r\n t1   0.5  \r\n").table
        assertEquals(listOf("Taxa", "Q1"), table.attributes.map { it.name })
        assertEquals(listOf("t1"), table.taxa.names)
    }

    @Test
    fun aCovariateWithoutATraitLineIsAnError() = runTest {
        listOf("<Covariate>\nt1\t1\n", "<Factor>\nt1\tA\n", "<Header name=loc>\tA\nt1\t1\n").forEach { text ->
            val failure = assertFailsWith<IllegalArgumentException>(text) { read(text) }
            assertTrue(failure.message!!.endsWith("needs a <Trait> line that names the columns"), failure.message)
        }
    }

    @Test
    fun readsAHeaderFileAsLongFormat() = runTest {
        val table = readFile("phenotype/version4format_test.txt").table
        assertEquals(listOf("Taxa", "loc", "trait1", "trait2"), table.attributes.map { it.name })
        assertEquals(
            listOf(AttributeType.taxa, AttributeType.factor, AttributeType.data, AttributeType.data),
            table.attributes.map { it.type },
        )
        assertEquals(6, table.numberOfObservations)
        assertEquals(
            listOf("taxon1", "taxon1", "taxon2", "taxon2", "taxon3", "taxon3"),
            List(6) { (table.attribute("Taxa") as TaxaAttribute).taxon(it) },
        )
        val loc = table.attribute("loc") as FactorAttribute
        assertEquals(listOf("A", "B", "A", "B", "A", "B"), List(6) { loc.level(it) })
        val trait2 = table.attribute("trait2") as NumericAttribute
        assertEquals(listOf(21f, 24f, 22f, 25f, 23f, 26f), List(6) { trait2.floatValue(it) })
    }

    @Test
    fun aHeaderCombinationTheFileDoesNotGiveForATraitIsMissing() = runTest {
        val table = read(
            "<Trait>\tx\ty\tx\n<Header name=loc>\tA\tA\tB\nt1\t1\t2\t3\nt2\t4\t5\t6\n",
        ).table
        assertEquals(listOf("Taxa", "loc", "x", "y"), table.attributes.map { it.name })
        val x = table.attribute("x") as NumericAttribute
        val y = table.attribute("y") as NumericAttribute
        assertEquals(listOf(1f, 3f, 4f, 6f), List(4) { x.floatValue(it) })
        assertEquals(listOf(false, true, false, true), List(4) { y.isMissing(it) })
        assertEquals(2, y.missingCount)
    }

    @Test
    fun twoHeaderLinesGiveTwoFactorsAndTheirDistinctCombinations() = runTest {
        val table = read(
            "<Covariate>\n<Trait>\ty\ty\ty\n<Header name=env>\tIA\tIA\tNY\n<Header name=rep>\t1\t2\t1\nt1\t1\t2\t3\n",
        ).table
        assertEquals(listOf("Taxa", "env", "rep", "y"), table.attributes.map { it.name })
        assertEquals(AttributeType.covariate, table.attribute("y").type)
        val env = table.attribute("env") as FactorAttribute
        val rep = table.attribute("rep") as FactorAttribute
        val y = table.attribute("y") as NumericAttribute
        assertEquals(listOf("IA", "IA", "NY"), List(3) { env.level(it) })
        assertEquals(listOf("1", "2", "1"), List(3) { rep.level(it) })
        assertEquals(listOf(1f, 2f, 3f), List(3) { y.floatValue(it) })
    }

    @Test
    fun keepsTraitsAndLevelsInFirstAppearanceOrder() = runTest {
        val table = read("<Trait>\tz\ta\tz\ta\n<Header name=loc>\tB\tB\tA\tA\nt1\t1\t2\t3\t4\n").table
        assertEquals(listOf("Taxa", "loc", "z", "a"), table.attributes.map { it.name })
        val loc = table.attribute("loc") as FactorAttribute
        assertEquals(listOf("B", "A"), List(2) { loc.level(it) })
        assertEquals(listOf("A", "B"), loc.levels)
    }

    @Test
    fun readsMinus999AsMissingInAHeaderFileToo() = runTest {
        val table = read("<Trait>\tx\tx\n<Header name=loc>\tA\tB\nt1\t-999\t5\n").table
        assertTrue(table.attribute("x").isMissing(0))
    }

    @Test
    fun explainsEveryMalformedFile() = runTest {
        val cases = mapOf(
            "<Trait>\tA\tB\nt1\t1\n" to
                "Line 2 has 2 values, but the <Trait> line names 2 traits, so 3 are expected: 't1\t1'",
            "<Trait>\tA\nt1\t1\t2\n" to
                "Line 2 has 3 values, but the <Trait> line names 1 traits, so 2 are expected: 't1\t1\t2'",
            "<Trait>\nt1\n" to "Line 1: <Trait> names no traits",
            "<Trait>\tA\n<Trait>\tB\n" to "Line 2: a second <Trait> line; line 1 has one",
            "<Trait>\tA\nt1\t1\n<Factor>\n" to
                "Line 3: the tag '<Factor>' comes after data; every tag must come before the first data line",
            "<Covariate>\n<Factor>\n<Trait>\tA\nt1\tx\n" to
                "<Covariate> and <Factor> both appear; a file's columns have one type",
            "<Trait>\tA\n<Factor>\n<Header name=x>\tA\n" to
                "<Factor> and <Header> both appear; a <Header> file's columns are numeric",
            "<Trait>\tA\tB\n<Header name=loc>\tA\nt1\t1\t2\n" to
                "Line 2: <Header name=loc> gives 1 levels, but the <Trait> line names 2 traits",
            "<Trait>\tA\tB\n<Header name=loc>\tA\tB\tC\nt1\t1\t2\n" to
                "Line 2: <Header name=loc> gives 3 levels, but the <Trait> line names 2 traits",
            "<Trait>\tA\tB\n<Header name=loc>\tA\tB\n<Header name=env>\tx\nt1\t1\t2\n" to
                "Line 3: <Header name=env> gives 1 levels, but the <Trait> line names 2 traits",
            "<Trait>\tA\tA\n<Header name=loc>\tx\tx\nt1\t1\t2\n" to "Column 3 repeats trait 'A' for loc=x",
            "<Trait>\tA\n<Header name=loc>\tx\n<Header name=loc>\ty\n" to "<Header name=loc> appears more than once",
            "<Trait>\tA\n<Header loc>\tx\n" to
                "Line 2: malformed header '<Header loc>\tx'; expected <Header name=factor> followed by one " +
                "level per trait",
            "<Marker>\nt1\n" to "Line 1: not a TASSEL phenotype file; it must start with <Phenotype>, <Trait>, " +
                "<Covariate>, <Factor>, or <Header>, found '<Marker>'",
            "<Trait>\tA\n<Numeric>\nt1\t1\n" to "Line 2: unknown tag '<Numeric>'; the accepted tags are " +
                "<Phenotype>, <Trait>, <Covariate>, <Factor>, and <Header>",
        )
        cases.forEach { (text, expected) ->
            val failure = assertFailsWith<IllegalArgumentException>(text) { read(text) }
            assertEquals(expected, failure.message, text)
        }
    }

    @Test
    fun anInfiniteNumberIsMissingWithAWarning() = runTest {
        val read = read("<Trait>\tA\nt1\tInfinity\nt2\t1\n")
        assertTrue(read.table.attribute("A").isMissing(0))
        assertEquals(listOf("1 infinite values were read as missing"), read.warnings)
    }

    @Test
    fun readsWhatItWritesBackAsAPhenotypeFile() = runTest {
        val table = readFile("phenotype/version4format_test.txt").table
        assertEquals(table, read(table.toDelimitedString()).table)
    }

    private companion object {
        /** Each golden and the input TASSEL 5 exported it from. */
        val GOLDENS = mapOf(
            "golden/phenotype-legacy/mdp_traits.txt" to "tutorial/mdp_traits.txt",
            "golden/phenotype-legacy/mdp_population_structure.txt" to "tutorial/mdp_population_structure.txt",
            "golden/phenotype-legacy/version4format_test.txt" to "phenotype/version4format_test.txt",
            "golden/phenotype-legacy/version4format_covariate.txt" to "phenotype/version4format_covariate.txt",
            "golden/phenotype-legacy/mdp_traits_nomissing.txt" to "phenotype/mdp_traits_nomissing.txt",
        )

        suspend fun read(text: String, options: PhenotypeReadOptions = PhenotypeReadOptions()) =
            PhenotypeFiles.read(Buffer().apply { writeString(text) }, options)

        suspend fun readFile(path: String) = PhenotypeFiles.read(Buffer().apply { write(TestData.readBytes(path)) })
    }
}
