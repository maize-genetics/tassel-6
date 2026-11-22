package net.maizegenetics.tassel.table

import net.maizegenetics.tassel.testing.goldenText
import net.maizegenetics.tassel.testing.normalizeLineEndings
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * The `Phenotype` export format, against the shape of `testdata/golden/pca/pcs.txt`.
 *
 * @see net.maizegenetics.phenotype.PhenotypeUtils.write
 */
class PhenotypeTableTest {
    @Test
    fun writesTheThreeHeaderLinesAndOneRowPerObservation() {
        assertEquals(
            "<Phenotype>\n" +
                "taxa\tcovariate\tcovariate\n" +
                "Taxa\tPC1\tPC2\n" +
                "one\t0.5\t-1.25\n" +
                "two\t2.0\t3.5\n",
            table.toDelimitedString(),
        )
    }

    @Test
    fun writesFloatsInJavasSpellingOnEveryTarget() {
        val tiny = PhenotypeTable(
            listOf(PhenotypeTable.taxa(listOf("one")), PhenotypeTable.Attribute("PC1", floatArrayOf(3.1861675e-4f))),
        )
        assertTrue(tiny.toDelimitedString().endsWith("one\t3.1861675E-4\n"), tiny.toDelimitedString())
    }

    /** A `Double` still goes through `DoubleFormat`, as `saveDelimitedTableReport` does. */
    @Test
    fun formatsCellsTheWayTassel5Does() {
        assertEquals("", PhenotypeTable.formatCell(null))
        assertEquals("0.12346", PhenotypeTable.formatCell(0.123456789))
        assertEquals("0.12345679", PhenotypeTable.formatCell(0.123456789f))
        assertEquals("7", PhenotypeTable.formatCell(7))
        assertEquals("name", PhenotypeTable.formatCell("name"))
    }

    @Test
    fun rejectsAttributesOfDifferentLengthsAndUnusableCells() {
        assertFailsWith<IllegalArgumentException> {
            PhenotypeTable(
                listOf(PhenotypeTable.taxa(listOf("one", "two")), PhenotypeTable.Attribute("PC1", floatArrayOf(1.0f))),
            )
        }
        assertFailsWith<IllegalArgumentException> { PhenotypeTable(emptyList()) }
        assertFailsWith<IllegalArgumentException> {
            PhenotypeTable.Attribute("PC1", PhenotypeTable.AttributeType.covariate, listOf(true))
        }
    }

    /** The format the PCA golden is written in, read back off the golden itself. */
    @Test
    fun matchesTheShapeOfTheGoldenPhenotype() {
        val lines = normalizeLineEndings(goldenText("golden/pca/pcs.txt")).lines()
        assertEquals(PhenotypeTable.HEADER, lines[0])
        assertEquals(listOf(PhenotypeTable.AttributeType.taxa.name) + List(5) { "covariate" }, lines[1].split("\t"))
        assertEquals(listOf(PhenotypeTable.TAXA) + (1..5).map { "PC$it" }, lines[2].split("\t"))
    }

    private val table = PhenotypeTable(
        listOf(
            PhenotypeTable.taxa(listOf("one", "two")),
            PhenotypeTable.Attribute("PC1", floatArrayOf(0.5f, 2.0f)),
            PhenotypeTable.Attribute("PC2", floatArrayOf(-1.25f, 3.5f)),
        ),
    )
}
