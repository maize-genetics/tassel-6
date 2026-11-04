package net.maizegenetics.tassel.table

import net.maizegenetics.tassel.testing.assertMatchesGolden
import net.maizegenetics.tassel.testing.goldenText
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class ResultTableTest {
    @Test
    fun writesTheAlleleSummaryGoldenFromItsValues() {
        val written = alleleSummary().toDelimitedString()
        assertMatchesGolden("golden/summary/alleles.txt", written)
        // The golden tolerance is exact, which normalizes line endings; TASSEL writes \n only, and
        // a line after the last row.
        assertEquals(goldenText("golden/summary/alleles.txt"), written)
        assertTrue(written.startsWith("Alleles\tNumber\tProportion\tFrequency\n"), "header line")
        assertTrue(written.endsWith("T:T\t19\t0.00614\t\n"), "last row, empty cell, trailing newline")
        assertTrue('\r' !in written, "TASSEL writes \\n endings")
    }

    @Test
    fun writesTassel5CellText() {
        assertEquals("", ResultTable.formatCell(null))
        assertEquals("0.28342", ResultTable.formatCell(0.2834169223812696))
        assertEquals("281", ResultTable.formatCell(281.0))
        assertEquals("281", ResultTable.formatCell(281))
        assertEquals("869133", ResultTable.formatCell(869133L))
        assertEquals("PZB00859.1", ResultTable.formatCell("PZB00859.1"))
    }

    @Test
    fun writesHeaderAndDelimiterAsAsked() {
        val table = resultTable("Overview", listOf("Stat Type", "Value")) {
            row("Number of Taxa", 281.0)
            row("Missing", null)
        }
        assertEquals("Stat Type\tValue\nNumber of Taxa\t281\nMissing\t\n", table.toDelimitedString())
        assertEquals("Number of Taxa,281\nMissing,\n", table.toDelimitedString(",", includeHeader = false))

        val out = StringBuilder("header: ")
        table.writeDelimited(out, delimiter = "|", includeHeader = false)
        assertEquals("header: Number of Taxa|281\nMissing|\n", out.toString())
    }

    @Test
    fun readsItsOwnShape() {
        val table = alleleSummary()
        assertEquals("Allele Summary", table.title)
        assertEquals(27, table.rowCount)
        assertEquals(4, table.columnCount)
        assertEquals("C", table[0, 0])
        assertEquals(246327, table[0, 1])
        assertEquals(null, table[11, 3])
        assertEquals(2, table.columnIndex("Proportion"))
        assertEquals(-1, table.columnIndex("Nope"))
    }

    @Test
    fun rejectsRaggedRows() {
        val problem = assertFailsWith<IllegalArgumentException> {
            resultTable("Overview", listOf("Stat Type", "Value")) { row("Number of Taxa") }
        }
        assertTrue("expected 2" in problem.message.orEmpty(), problem.message.orEmpty())
    }

    @Test
    fun rejectsACellTypeThatWouldPrintItsOwnDigits() {
        // A Float would print Float.toString, not DoubleFormat, and quietly break parity.
        val problem = assertFailsWith<IllegalArgumentException> {
            ResultTable("PCs", listOf("Taxa", "PC1"), listOf(listOf("33-16", 0.9156939f)))
        }
        assertTrue("Float" in problem.message.orEmpty(), problem.message.orEmpty())
    }

    /** `testdata/golden/summary/alleles.txt`, with the values from its `.full.txt` twin. */
    private fun alleleSummary(): ResultTable =
        resultTable("Allele Summary", listOf("Alleles", "Number", "Proportion", "Frequency")) {
            row("C", 246327, 0.2834169223812696, 0.29404384748162277)
            row("G", 235079, 0.27047528974276663, 0.2806169588479233)
            row("T", 178046, 0.2048547230400871, 0.2125359009313352)
            row("A", 168648, 0.1940416484013379, 0.20131738213870473)
            row("N", 31411, 0.036140613692035625, 0.037495732474496314)
            row("Y", 3938, 0.004530952109746149, 0.0047008434779079455)
            row("R", 2698, 0.003104242963965239, 0.003220638827677917)
            row("S", 994, 0.0011436684604082458, 0.0011865511470392327)
            row("K", 839, 9.653298171856321E-4, 0.0010015255657604791)
            row("W", 596, 6.857408474882441E-4, 7.114532028524976E-4)
            row("M", 557, 6.408685437096509E-4, 6.64898379175908E-4)
            row("C:T", 547, 0.17685095376656967, null)
            row("G:A", 486, 0.1571290009699321, null)
            row("T:C", 408, 0.13191076624636275, null)
            row("A:G", 373, 0.12059489169091497, null)
            row("G:C", 219, 0.07080504364694472, null)
            row("C:A", 184, 0.05948916909149693, null)
            row("G:T", 165, 0.0533462657613967, null)
            row("C:G", 155, 0.05011315874555448, null)
            row("A:T", 121, 0.03912059489169092, null)
            row("T:A", 110, 0.035564177174264466, null)
            row("T:G", 108, 0.03491755577109602, null)
            row("A:C", 78, 0.02521823472356935, null)
            row("C:C", 50, 0.016165535079211122, null)
            row("G:G", 46, 0.014872292272874232, null)
            row("A:A", 24, 0.007759456838021339, null)
            row("T:T", 19, 0.006142903330100227, null)
        }
}
