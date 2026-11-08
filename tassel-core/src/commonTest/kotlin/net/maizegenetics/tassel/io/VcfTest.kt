package net.maizegenetics.tassel.io

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import kotlinx.io.readLine
import kotlinx.io.readString
import kotlinx.io.writeString
import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.testing.GoldenManifest
import net.maizegenetics.tassel.testing.SuspendLazy
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.Tolerance
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

/**
 * The tutorial VCF read and written back as HapMap, against the TASSEL 5 export in `testdata/golden/vcf/`
 * (ADR-0009), plus the parsing rules on small hand-written inputs.
 */
class VcfTest {
    @Test
    fun readsTutorialDimensions() = runTest {
        val table = tutorial()
        assertEquals(281, table.numberOfTaxa)
        assertEquals(3093, table.numberOfSites)
        assertEquals("33-16", table.taxa[0].name)
        assertEquals("YU796NS", table.taxa[280].name)
        assertEquals((1..10).map { it.toString() }, table.positions.chromosomes.map { it.name })
    }

    @Test
    fun readsTutorialPositions() = runTest {
        val table = tutorial()
        val first = table.positions[0]
        assertEquals("PZA01271.1", first.name)
        assertEquals("1", first.chromosome.name)
        assertEquals(157104, first.position)
        // REF=C, ALT=G. TASSEL 5 also stores C as the position's reference allele and G as its alternate; Position
        // has no such slots (#10), so only the variant text survives (see docs/porting-ledger.md).
        assertEquals(listOf("C", "G"), first.knownVariants)
        val last = table.positions[3092]
        assertEquals("PZA03429.1", last.name)
        assertEquals("10", last.chromosome.name)
        assertEquals(148907116, last.position)
    }

    @Test
    fun readsTutorialGenotypes() = runTest {
        val table = tutorial()
        // Site 0 (REF=C, ALT=G) reads 0/0 for 33-16, 1/1 for 38-11, and ./. for CI3A.
        assertEquals(Genotype.of(Allele.C, Allele.C), table.genotype(0, 0))
        assertEquals(Genotype.of(Allele.G, Allele.G), table.genotype(1, 0))
        assertEquals("CI3A", table.taxa[TAXON_CI3A].name)
        assertEquals(Genotype.UNKNOWN, table.genotype(TAXON_CI3A, 0))
        assertEquals('N', Nucleotide.iupac(table.genotype(TAXON_CI3A, 0)))
    }

    @Test
    fun writesTassel5HapMapExport() = runTest {
        assertEquals(Tolerance.Exact(), GoldenManifest.load().tolerance(GOLDEN))
        val written = Buffer().also { HapMap.write(tutorial(), it) }.readString().lines()
        val goldenLines = TestData.readText(GOLDEN).replace("\r\n", "\n").lines()

        // Maintainer decision (tassel-6#14, 2026-11-04): compare from the rs# header row down. The golden's first
        // 281 lines are ##SAMPLE= taxa annotations, which TASSEL 5 copies out of the VCF into Taxon annotations and
        // re-emits on export. Taxon is frozen as a name (#10) and HapMap.write emits no ## lines, so they cannot
        // round trip; the gap is recorded in docs/porting-ledger.md. This is a decision, not an oversight: do not
        // "fix" it by editing the golden or the MANIFEST, whose tolerance for this output stays exact. The shape of
        // the skipped block is asserted so the test still fails if the golden or the taxa count changes.
        val skipped = goldenLines.takeWhile { it.startsWith("##") }
        assertEquals(tutorial().numberOfTaxa, skipped.size, "skipped ##SAMPLE= lines")
        assertTrue(skipped.all { it.startsWith("##SAMPLE=<ID=") }, "every skipped line is a ##SAMPLE= annotation")
        val golden = goldenLines.drop(skipped.size)
        assertTrue(golden.first().startsWith("rs#\talleles\tchrom\tpos\t"), "first compared line: ${golden.first()}")

        assertEquals(golden.size, written.size, "number of lines")
        written.forEachIndexed { index, line ->
            if (line != golden[index]) reportLineDifference(skipped.size + index, golden[index], line)
        }
    }

    @Test
    fun readsUnphasedAndPhasedAndMissingCalls() = runTest {
        // Phase is not kept: the genotype byte has no phase bit, so 0|1 and 0/1 are the same afterwards.
        val table = read(header("T1", "T2", "T3", "T4", "T5"), site("s1", "A", "C", "0/0", "0/1", "0|1", "1|1", "./."))
        assertEquals(Nucleotide.genotype("AA"), table.genotype(0, 0))
        assertEquals(Genotype.of(Allele.A, Allele.C), table.genotype(1, 0))
        assertEquals(table.genotype(1, 0), table.genotype(2, 0))
        assertEquals(Nucleotide.genotype("CC"), table.genotype(3, 0))
        assertEquals(Genotype.UNKNOWN, table.genotype(4, 0))
    }

    @Test
    fun readsTheMissingFormsTassel5Reads() = runTest {
        // TASSEL 5 treats any taxon column starting with a dot as missing, and a negative allele index as missing.
        val table = read(header("T1", "T2", "T3", "T4"), site("s1", "A", "C", ".", ".|.", ".:5", "0/."))
        (0..3).forEach { assertEquals(Genotype.UNKNOWN, table.genotype(it, 0), "taxon $it") }
    }

    @Test
    fun readsASingleCharacterCallAsHomozygous() = runTest {
        val table = read(header("T1", "T2"), site("s1", "A", "C", "0", "1"))
        assertEquals(Nucleotide.genotype("AA"), table.genotype(0, 0))
        assertEquals(Nucleotide.genotype("CC"), table.genotype(1, 0))
    }

    @Test
    fun readsMultiallelicSites() = runTest {
        val table = read(header("T1", "T2"), site("s1", "A", "C,G", "0/2", "1/2"))
        assertEquals(listOf("A", "C", "G"), table.positions[0].knownVariants)
        assertEquals(Genotype.of(Allele.A, Allele.G), table.genotype(0, 0))
        assertEquals(Genotype.of(Allele.C, Allele.G), table.genotype(1, 0))
    }

    @Test
    fun readsASiteWithNoAlternateAllele() = runTest {
        val table = read(header("T1"), site("s1", "A", ".", "0/0"))
        assertEquals(listOf("A"), table.positions[0].knownVariants)
        assertEquals(Nucleotide.genotype("AA"), table.genotype(0, 0))
    }

    @Test
    fun readsIndelsTheWayTassel5Does() = runTest {
        // REF longer than ALT: TASSEL 5 drops the shared first base, so AT/A becomes T and a gap.
        val deletion = read(header("T1", "T2"), site("s1", "AT", "A", "0/0", "1/1"))
        assertEquals(listOf("AT", "A"), deletion.positions[0].knownVariants)
        assertEquals(Nucleotide.genotype("TT"), deletion.genotype(0, 0))
        assertEquals(Nucleotide.genotype("--"), deletion.genotype(1, 0))

        // ALT longer than REF: the same rule, so A/AT becomes a gap and T.
        val insertion = read(header("T1", "T2"), site("s1", "A", "AT", "0/0", "0/1"))
        assertEquals(Nucleotide.genotype("--"), insertion.genotype(0, 0))
        assertEquals(Genotype.of(Allele.GAP, Allele.T), insertion.genotype(1, 0))
    }

    @Test
    fun readsSymbolicAlleles() = runTest {
        // TASSEL 5 rewrites <INS> to +, <DEL> to -, and * to N before taking the first character of each variant.
        val table = read(header("T1", "T2", "T3"), site("s1", "A", "<INS>,<DEL>", "0/0", "1/1", "2/2"))
        assertEquals(Nucleotide.genotype("AA"), table.genotype(0, 0))
        assertEquals(Nucleotide.genotype("++"), table.genotype(1, 0))
        assertEquals(Nucleotide.genotype("--"), table.genotype(2, 0))

        val star = read(header("T1", "T2"), site("s1", "A", "*", "0/0", "1/1"))
        assertEquals(Nucleotide.genotype("AA"), star.genotype(0, 0))
        assertEquals(Genotype.UNKNOWN, star.genotype(1, 0))
    }

    @Test
    fun namesASiteWithoutAnIdTheWayTassel5Does() = runTest {
        val table = read(header("T1"), site(".", "A", "C", "0/0"))
        assertEquals("S1_10", table.positions[0].name)
        assertEquals(null, table.positions[0].snpId)
    }

    @Test
    fun skipsAnnotationLinesAndCommentsInTheBody() = runTest {
        val table = read(
            "##fileformat=VCFv4.0",
            "##SAMPLE=<ID=T1,FATHER=T1,INBREEDF=0.97,MOTHER=T1>",
            header("T1", "T2"),
            site("s1", "A", "C", "0/0", "1/1"),
            "#a comment in the body",
            site("s2", "A", "C", "0/1", "0/0"),
        )
        assertEquals(listOf("T1", "T2"), table.taxa.names)
        assertEquals(listOf("s1", "s2"), table.positions.map { it.name })
    }

    @Test
    fun sortsSitesThatAreOutOfOrder() = runTest {
        val table = read(
            header("T1"),
            site("s3", "A", "C", "0/0", chromosome = "2", position = "5"),
            site("s2", "A", "C", "1/1", chromosome = "1", position = "20"),
            site("s1", "A", "C", "0/1", chromosome = "1", position = "10"),
        )
        assertEquals(listOf("s1", "s2", "s3"), table.positions.map { it.name })
        assertEquals(
            listOf(Genotype.of(Allele.A, Allele.C), Nucleotide.genotype("CC"), Nucleotide.genotype("AA")),
            (0..2).map { table.genotype(0, it) },
        )
    }

    @Test
    fun readsAHeaderWithoutAnIdColumn() = runTest {
        // HeaderPositions finds columns by name, so a file may leave any of them out.
        val table = read(
            "#CHROM\tPOS\tREF\tALT\tINFO\tFORMAT\tT1",
            "1\t10\tA\tC\tNS=1\tGT\t0/1",
        )
        assertEquals("S1_10", table.positions[0].name)
        assertEquals(Genotype.of(Allele.A, Allele.C), table.genotype(0, 0))
    }

    @Test
    fun rejectsMalformedInput() = runTest {
        assertFailsWith<IllegalArgumentException> { read() }
        assertFailsWith<IllegalArgumentException> { read("##only annotations") }
        // A #CHROM header with no taxon column.
        assertFailsWith<IllegalArgumentException> { read(FIXED_HEADER) }
        // A header missing REF.
        assertFailsWith<IllegalArgumentException> {
            read("#CHROM\tPOS\tID\tALT\tINFO\tFORMAT\tT1", "1\t10\t.\tC\tNS=1\tGT\t0/0")
        }
        assertFailsWith<IllegalArgumentException> { read(header("T1"), site("s1", "A", "C", "0/0", position = "ten")) }
        // Fewer genotype columns than taxa, and more.
        assertFailsWith<IllegalArgumentException> { read(header("T1", "T2"), site("s1", "A", "C", "0/0")) }
        assertFailsWith<IllegalArgumentException> { read(header("T1"), site("s1", "A", "C", "0/0", "1/1")) }
        // A site line that stops before the genotypes start.
        assertFailsWith<IllegalArgumentException> { read(header("T1"), "1\t10\ts1\tA\tC\t.\tPASS\tNS=1") }
    }

    @Test
    fun rejectsAFormatWithoutGtFirst() = runTest {
        val misplaced = assertFailsWith<IllegalArgumentException> {
            read(header("T1"), "1\t10\ts1\tA\tC\t.\tPASS\tNS=1\tAD:GT\t2,0:0/0")
        }
        assertContains(misplaced.message ?: "", "first FORMAT field")
        val absent = assertFailsWith<IllegalArgumentException> {
            read(header("T1"), "1\t10\ts1\tA\tC\t.\tPASS\tNS=1\tAD:DP\t2,0:2")
        }
        assertContains(absent.message ?: "", "no GT field")
    }

    @Test
    fun rejectsACallNamingAnAlleleTheSiteDoesNotList() = runTest {
        val failure = assertFailsWith<IllegalArgumentException> {
            read(header("T1"), site("s1", "A", "C", "0/2"))
        }
        assertContains(failure.message ?: "", "s1")
    }

    @Test
    fun doesNotReadAlleleIndicesOfTenOrMore() = runTest {
        // TASSEL 5 reads only the characters at offsets 0 and 2 of the GT field, so 10/10 is read as allele 1 and
        // then '/' - '0', which is -1 and therefore missing. Matched deliberately rather than fixed (tassel-6#14).
        val alts = (1..10).joinToString(",") { listOf("C", "G", "T")[it % 3] }
        val table = read(header("T1"), site("s1", "A", alts, "10/10"))
        assertEquals(11, table.positions[0].knownVariants.size)
        assertEquals(Genotype.UNKNOWN, table.genotype(0, 0))
    }

    @Test
    fun keepsContigIdsAndOptionalLengthsInOrder() = runTest {
        val sample = "##SAMPLE=<ID=T1,FATHER=T1>"
        val lines = arrayOf(
            "##contig=<ID=1,length=100,assembly=ref>",
            sample,
            "##contig=<length=9>",
            "##contig=<ID=2>",
            "##contig=<ID=\"chr 3\",length=12>",
            header("T1"),
            site("s1", "A", "C", "0/1"),
        )
        val bytes = lines.joinToString("\n").plus("\n").encodeToByteArray()
        val source = Buffer().also { it.write(bytes) }
        val parsed = Vcf.readHeader(source)
        assertEquals(
            listOf(VcfContig("1", 100L), VcfContig("2", null), VcfContig("chr 3", 12L)),
            parsed.contigs,
        )
        // The header reader stops on the #CHROM line, so the next line is the first record.
        assertEquals(site("s1", "A", "C", "0/1"), source.readLine())

        val table = Vcf.read(Buffer().also { it.write(bytes) })
        assertEquals(listOf("T1"), table.taxa.names)
        assertEquals(Genotype.of(Allele.A, Allele.C), table.genotype(0, 0))
    }

    @Test
    fun parsesARecordSplitAtEveryByteOffset() = runTest {
        // `\n` and `\r\n` both end a line. A split on the `\r` of `\r\n` has no line break yet.
        for (lineBreak in listOf("\n", "\r\n")) {
            val text = header("T1", "T2") + lineBreak + site("s1", "A", "C", "0/1", "1/1") + lineBreak
            val bytes = text.encodeToByteArray()
            val headerEnd = text.indexOf(lineBreak) + lineBreak.length
            val record = bytes.copyOfRange(headerEnd, bytes.size)
            val parsedHeader = Vcf.readHeader(Buffer().also { it.write(bytes, endIndex = headerEnd) })
            val table = Vcf.read(Buffer().also { it.write(bytes) })
            val expected = parseRecords(parsedHeader, record, 0, record.size).single()
            assertEquals(table.positions[0], expected.position, "whole read, break ${lineBreak.length}")
            assertContentEquals(table.siteGenotypes(0), expected.calls, "whole read, break ${lineBreak.length}")
            val newline = record.indexOf('\n'.code.toByte())
            for (split in 0..record.size) {
                val parser = VcfRecordParser(parsedHeader)
                val calls = ByteArray(parsedHeader.taxaNames.size)
                val first = parser.tryParse(record, 0, split, eof = false, calls)
                val actual = if (split <= newline) {
                    assertTrue(first is VcfParse.Incomplete, "split $split should need another refill")
                    parser.tryParse(record, 0, record.size, eof = true, calls)
                } else {
                    first
                }
                val site = actual as VcfParse.Site
                assertEquals(expected.position, site.position, "position at split $split")
                assertEquals(expected.position.knownVariants, site.position.knownVariants)
                assertContentEquals(expected.calls, calls, "genotypes at split $split")
            }
        }
    }

    @Test
    fun parserStartedAtAMiddleRecordMatchesTheRestOfTheFile() = runTest {
        val lines = arrayOf(
            header("T1", "T2"),
            site("s1", "A", "C", "0/0", "1/1", position = "10"),
            site("s2", "A", "G", "0/1", "1/1", position = "20"),
            site("s3", "C", "T", "1/1", "0/0", position = "30"),
        )
        val text = lines.joinToString("\n") + "\n"
        val bytes = text.encodeToByteArray()
        val headerEnd = text.indexOf('\n') + 1
        val secondSite = text.indexOf('\n', headerEnd) + 1

        val whole = read(*lines)
        val suffix = Vcf.read(
            Buffer().also {
                it.write(bytes, endIndex = headerEnd)
                it.write(bytes, startIndex = secondSite)
            },
        )
        val parsedHeader = Vcf.readHeader(Buffer().also { it.write(bytes, endIndex = headerEnd) })
        val parsed = parseRecords(parsedHeader, bytes, secondSite, bytes.size)

        assertEquals(listOf("s2", "s3"), parsed.map { it.position.name })
        assertEquals(suffix.positions.map { it.name }, parsed.map { it.position.name })
        assertEquals(whole.positions.drop(1).map { it.name }, parsed.map { it.position.name })
        parsed.forEachIndexed { index, site ->
            assertEquals(suffix.positions[index].knownVariants, site.position.knownVariants)
            assertEquals(whole.positions[index + 1].knownVariants, site.position.knownVariants)
            repeat(suffix.numberOfTaxa) { taxon ->
                assertEquals(suffix.genotype(taxon, index), site.calls[taxon], "suffix taxon $taxon")
                assertEquals(whole.genotype(taxon, index + 1), site.calls[taxon], "whole taxon $taxon")
            }
        }
    }

    @Test
    fun readsARecordLongerThanTheFillBuffer() = runTest {
        val info = "N".repeat(70_000)
        val table = read(header("T1"), "1\t10\ts1\tA\tC\t.\tPASS\t$info\tGT\t0/1")
        assertEquals("s1", table.positions[0].name)
        assertEquals(Genotype.of(Allele.A, Allele.C), table.genotype(0, 0))
    }

    @Test
    fun readsBackWhatHapMapWrote() = runTest {
        val written = Buffer().also { HapMap.write(tutorial(), it) }.readString()
        val reread = HapMap.read(Buffer().also { it.writeString(written) })
        assertEquals(tutorial().taxa.names, reread.taxa.names)
        assertEquals(tutorial().positions.map { it.name }, reread.positions.map { it.name })
    }

    private suspend fun read(vararg lines: String): GenotypeTable =
        Vcf.read(Buffer().also { buffer -> lines.forEach { buffer.writeString(it + "\n") } })

    private fun header(vararg taxa: String): String = (Vcf.FIXED_COLUMNS + taxa.toList()).joinToString("\t")

    private fun site(
        snpId: String,
        reference: String,
        alternate: String,
        vararg calls: String,
        chromosome: String = "1",
        position: String = "10",
        format: String = "GT",
    ) = (
        listOf(chromosome, position, snpId, reference, alternate, ".", "PASS", "NS=1", format) + calls.toList()
        ).joinToString("\t")

    private fun parseRecords(header: VcfHeader, window: ByteArray, from: Int, until: Int): List<ParsedRecord> {
        val parser = VcfRecordParser(header)
        val calls = ByteArray(header.taxaNames.size)
        val records = ArrayList<ParsedRecord>()
        var at = from
        while (true) {
            when (val parsed = parser.tryParse(window, at, until, eof = true, calls)) {
                is VcfParse.End -> return records
                is VcfParse.Skip -> at += parsed.consumed
                is VcfParse.Site -> {
                    records.add(ParsedRecord(parsed.position, calls.copyOf()))
                    at += parsed.consumed
                }
                VcfParse.Incomplete -> error("a full window at EOF should hold complete lines")
            }
        }
    }

    private class ParsedRecord(val position: Position, val calls: ByteArray)

    private fun reportLineDifference(line: Int, expected: String, actual: String): Nothing {
        val column = expected.commonPrefixWith(actual).length
        throw AssertionError(
            "Line ${line + 1} differs at column ${column + 1}:\n" +
                "  expected: ...${expected.drop(column).take(40)}\n" +
                "  actual:   ...${actual.drop(column).take(40)}",
        )
    }

    private companion object {
        const val TUTORIAL = "tutorial/mdp_genotype_w_depth.vcf"
        const val GOLDEN = "golden/vcf/mdp_genotype_w_depth.hmp.txt"
        const val TAXON_CI3A = 57
        val FIXED_HEADER = Vcf.FIXED_COLUMNS.joinToString("\t")

        val tutorial = SuspendLazy {
            Vcf.read(Buffer().also { it.write(TestData.readBytes(TUTORIAL)) })
        }
    }
}
