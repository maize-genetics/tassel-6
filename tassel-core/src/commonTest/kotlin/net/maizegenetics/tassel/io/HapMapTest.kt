package net.maizegenetics.tassel.io

import kotlinx.io.Buffer
import kotlinx.io.readString
import kotlinx.io.writeString
import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.testing.GoldenManifest
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.Tolerance
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * The tutorial HapMap file read and written back, against the TASSEL 5 export in `testdata/golden/hapmap/`
 * (ADR-0009), plus the parsing rules on small hand-written inputs.
 */
class HapMapTest {
    @Test
    fun readsTutorialDimensions() {
        val table = readTutorial()
        assertEquals(281, table.numberOfTaxa)
        assertEquals(3093, table.numberOfSites)
        assertEquals("33-16", table.taxa[0].name)
        assertEquals("38-11", table.taxa[1].name)
        assertEquals("YU796NS", table.taxa[280].name)
        assertEquals((1..10).map { it.toString() }, table.positions.chromosomes.map { it.name })
    }

    @Test
    fun readsTutorialPositions() {
        val table = readTutorial()
        val first = table.positions[0]
        assertEquals("PZB00859.1", first.name)
        assertEquals("1", first.chromosome.name)
        assertEquals(157104, first.position)
        assertEquals(listOf("A", "C"), first.knownVariants)
        val last = table.positions[3092]
        assertEquals("PZA03429.1", last.name)
        assertEquals("10", last.chromosome.name)
        assertEquals(148907116, last.position)
    }

    @Test
    fun readsTutorialGenotypes() {
        val table = readTutorial()
        // Site 0 of the input reads CC for 33-16, AA for A188, NN for CML328, and AC for A654.
        assertEquals(Nucleotide.genotype("CC"), table.genotype(0, 0))
        assertEquals(Nucleotide.genotype("AA"), table.genotype(4, 0))
        assertEquals(Genotype.UNKNOWN, table.genotype(17, 0))
        // TASSEL 5 reads the two characters in reverse, so AC becomes C/A; both write as the IUPAC code M.
        assertEquals(Genotype.of(Allele.C, Allele.A), table.genotype(91, 0))
        assertEquals('M', Nucleotide.iupac(table.genotype(91, 0)))
    }

    @Test
    fun writesTassel5Export() {
        assertEquals(Tolerance.Exact(), GoldenManifest.load().tolerance(GOLDEN))
        val written = Buffer().also { HapMap.write(readTutorial(), it) }.readString().lines()
        val golden = TestData.readText(GOLDEN).replace("\r\n", "\n").lines()
        assertEquals(golden.size, written.size, "number of lines")
        written.forEachIndexed { index, line ->
            if (line != golden[index]) reportLineDifference(index, golden[index], line)
        }
    }

    @Test
    fun readingBackWhatWasWrittenGivesTheSameFile() {
        // One IUPAC character per genotype loses allele order, so A/C and C/A both come back as A/C; writing the
        // re-read table must still give the same file.
        val table = readTutorial()
        val written = Buffer().also { HapMap.write(table, it) }.readString()
        val reread = HapMap.read(Buffer().also { it.writeString(written) })
        assertEquals(table.taxa.names, reread.taxa.names)
        assertEquals(table.positions.map { it.name }, reread.positions.map { it.name })
        assertEquals(written, Buffer().also { HapMap.write(reread, it) }.readString())
    }

    @Test
    fun readsOneLetterGenotypes() {
        val table = read(header("T1", "T2"), site("s1", "A/C", "1", "10", "A", "M"))
        assertEquals(Nucleotide.genotype("AA"), table.genotype(0, 0))
        assertEquals(Nucleotide.genotype("AC"), table.genotype(1, 0))
    }

    @Test
    fun readsTwoLetterGenotypes() {
        val table = read(header("T1", "T2"), site("s1", "A/C", "1", "10", "AA", "CA"))
        assertEquals(Nucleotide.genotype("AA"), table.genotype(0, 0))
        // Reversed, as in TASSEL 5: the second character becomes the first allele.
        assertEquals(Nucleotide.genotype("AC"), table.genotype(1, 0))
    }

    @Test
    fun readsIndelsAndMissingGenotypes() {
        val table = read(header("T1", "T2", "T3"), site("s1", "+/-", "1", "10", "+", "0", "N"))
        assertEquals(Nucleotide.genotype("++"), table.genotype(0, 0))
        assertEquals(Nucleotide.genotype("-+"), table.genotype(1, 0))
        assertEquals(Genotype.UNKNOWN, table.genotype(2, 0))
    }

    @Test
    fun skipsAnnotationLines() {
        val table = read(
            "##fileformat=HapMap",
            "##SAMPLE=<ID=T1,Project=test>",
            header("T1", "T2"),
            site("s1", "A/C", "1", "10", "A", "C"),
        )
        assertEquals(listOf("T1", "T2"), table.taxa.names)
    }

    @Test
    fun sortsSitesThatAreOutOfOrder() {
        val table = read(
            header("T1"),
            site("s3", "A/C", "2", "5", "A"),
            site("s2", "A/C", "1", "20", "C"),
            site("s1", "A/C", "1", "10", "G"),
        )
        assertEquals(listOf("s1", "s2", "s3"), table.positions.map { it.name })
        assertEquals(listOf("G", "C", "A").map { Nucleotide.genotype(it[0]) }, (0..2).map { table.genotype(0, it) })
    }

    @Test
    fun writesTwoLetterGenotypesWhenAsked() {
        val table = read(header("T1", "T2"), site("s1", "A/C", "1", "10", "A", "M"))
        val written = Buffer().also { HapMap.write(table, it, diploid = true) }.readString()
        assertEquals("s1\tA/C\t1\t10\t+\tNA\tNA\tNA\tNA\tNA\tNA\tAA\tAC", written.lines()[1])
    }

    @Test
    fun writesNaWhenNoAlleleIsCalled() {
        val table = read(header("T1", "T2"), site("s1", "A/C", "1", "10", "N", "N"))
        val written = Buffer().also { HapMap.write(table, it) }.readString()
        assertEquals("s1\tNA\t1\t10\t+\tNA\tNA\tNA\tNA\tNA\tNA\tN\tN", written.lines()[1])
    }

    @Test
    fun writesTheHeaderTassel5Writes() {
        val table = read(header("T1", "T2"), site("s1", "A/C", "1", "10", "A", "C"))
        val written = Buffer().also { HapMap.write(table, it) }.readString()
        assertEquals(HapMap.FIXED_COLUMNS + listOf("T1", "T2"), written.lines()[0].split('\t'))
        assertContains(written, "panelLSID")
        assertTrue(written.endsWith("\n") && !written.contains("\r"))
    }

    @Test
    fun rejectsMalformedInput() {
        assertFailsWith<IllegalArgumentException> { read() }
        assertFailsWith<IllegalArgumentException> { read("##only annotations") }
        assertFailsWith<IllegalArgumentException> { read(HapMap.FIXED_COLUMNS.joinToString("\t")) }
        // Three characters per genotype is neither one nor two letters.
        assertFailsWith<IllegalArgumentException> {
            read(header("T1", "T2"), site("s1", "A/C", "1", "10", "AAA", "CCC"))
        }
        assertFailsWith<IllegalArgumentException> {
            read(header("T1", "T2"), site("s1", "A/C", "1", "ten", "A", "C"))
        }
        assertFailsWith<IllegalArgumentException> {
            read(header("T1", "T2"), site("s1", "A/C", "1", "10", "A", "C"), "s2\tA/C\t1\t20\tA")
        }
    }

    @Test
    fun rejectsGenotypesWithoutANucleotideMeaning() {
        val failure = assertFailsWith<IllegalArgumentException> {
            read(header("T1", "T2"), site("s1", "A/C", "1", "10", "A", "!"))
        }
        assertContains(failure.message ?: "", "s1")
    }

    @Test
    fun rejectsDelimitersTassel5Rejects() {
        val table = read(header("T1"), site("s1", "A/C", "1", "10", "A"))
        assertFailsWith<IllegalArgumentException> { HapMap.write(table, Buffer(), delimiter = ',') }
    }

    /** The tutorial file is read once and shared: parsing it per test is slow on Node and Native. */
    private fun readTutorial(): GenotypeTable = tutorial

    private fun read(vararg lines: String): GenotypeTable =
        HapMap.read(Buffer().also { buffer -> lines.forEach { buffer.writeString(it + "\n") } })

    private fun header(vararg taxa: String): String = (HapMap.FIXED_COLUMNS + taxa.toList()).joinToString("\t")

    private fun site(snpId: String, alleles: String, chromosome: String, position: String, vararg calls: String) =
        (listOf(snpId, alleles, chromosome, position, "+", "NA", "NA", "NA", "NA", "NA", "NA") + calls.toList())
            .joinToString("\t")

    private fun reportLineDifference(line: Int, expected: String, actual: String): Nothing {
        val column = expected.commonPrefixWith(actual).length
        throw AssertionError(
            "Line ${line + 1} differs at column ${column + 1}:\n" +
                "  expected: ...${expected.drop(column).take(40)}\n" +
                "  actual:   ...${actual.drop(column).take(40)}",
        )
    }

    private companion object {
        const val TUTORIAL = "tutorial/mdp_genotype.hmp.txt"
        const val GOLDEN = "golden/hapmap/mdp_genotype.hmp.txt"

        val tutorial: GenotypeTable by lazy {
            HapMap.read(Buffer().also { it.write(TestData.readBytes(TUTORIAL)) })
        }
    }
}
