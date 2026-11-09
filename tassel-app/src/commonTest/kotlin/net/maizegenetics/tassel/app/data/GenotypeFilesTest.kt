package net.maizegenetics.tassel.app.data

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.app.testing.ByteGenotypeFile
import net.maizegenetics.tassel.app.testing.FakeGenotypeFiles
import net.maizegenetics.tassel.app.testing.SMALL_HAPMAP
import net.maizegenetics.tassel.app.testing.TextGenotypeFile
import net.maizegenetics.tassel.app.testing.hapMapText
import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.io.Vcf
import net.maizegenetics.tassel.ops.Progress
import net.maizegenetics.tassel.ops.ProgressReporter
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * The [GenotypeFiles] contract, exercised through a byte-backed [GenotypeFile] like the browser's. The two real
 * implementations need a window or a page, so what is tested here is the part that does not: that a file which hands
 * over bytes rather than a path is enough for the reader, which is what makes the web target possible at all.
 */
class GenotypeFilesTest {
    @Test
    fun aFileThatOnlyHasBytesIsEnoughToReadATable() = runTest {
        val file = TextGenotypeFile("mdp.hmp.txt", SMALL_HAPMAP)

        val table = file.readGenotypes()

        assertEquals(listOf("B73", "Mo17", "CML247"), table.taxa.names)
        assertEquals(2, table.numberOfSites)
        assertEquals(Nucleotide.genotype('A'), table.genotype(0, 0))
        assertEquals(Nucleotide.genotype('M'), table.genotype(2, 0))
        assertEquals(Nucleotide.genotype('T'), table.genotype(2, 1))
        assertEquals(1, file.opened)
    }

    @Test
    fun theBundledTutorialFileIsNamedTheWayBothTargetsBundleIt() = runTest {
        val files = FakeGenotypeFiles(tutorial = TextGenotypeFile(TUTORIAL_FILE_NAME, SMALL_HAPMAP))

        val tutorial = files.tutorial()

        assertEquals("mdp_genotype.hmp.txt", TUTORIAL_FILE_NAME)
        assertEquals(TUTORIAL_FILE_NAME, tutorial.name)
        assertEquals(1, files.tutorialRequests)
    }

    @Test
    fun aCancelledChooserReturnsNoFile() = runTest {
        val files = FakeGenotypeFiles(chosen = null)

        assertNull(files.choose())
        assertEquals(1, files.chooseRequests)
    }

    @Test
    fun progressIsReportedToTheCaller() = runTest {
        // HapMap.read reports every 4,096 sites, and once at the end with the real total.
        val file = TextGenotypeFile("wide.hmp.txt", hapMapText(listOf("a", "b"), List(5000) { "AC" }))
        val reports = mutableListOf<Progress>()

        val table = file.readGenotypes(ProgressReporter { reports.add(it) })

        assertEquals(5000, table.numberOfSites)
        assertTrue(reports.isNotEmpty())
        assertEquals(Progress(5000, 5000, "Reading HapMap"), reports.last())
    }

    @Test
    fun aFileThatIsNotHapMapFails() = runTest {
        val file = TextGenotypeFile("notes.txt", "this is not a genotype file\n")

        val failure = assertFailsWith<IllegalArgumentException> { file.readGenotypes() }

        assertTrue(failure.message!!.isNotEmpty())
    }

    @Test
    fun aVcfNameReadsVcfAndAGzippedCopyReadsTheSameTable() = runTest {
        val plain = vcfText(listOf("B73", "Mo17"), listOf("0/0", "0/1")).encodeToByteArray()

        val vcf = ByteGenotypeFile("sample.vcf", plain).readGenotypes()
        val gzipped = ByteGenotypeFile("sample.vcf.gz", storedGzip(plain)).readGenotypes()

        assertEquals(listOf("B73", "Mo17"), gzipped.taxa.names)
        assertEquals(1, gzipped.numberOfSites)
        assertEquals(Nucleotide.genotype("AA"), gzipped.genotype(0, 0))
        assertEquals(Genotype.of(Allele.A, Allele.C), gzipped.genotype(1, 0))
        assertEquals(vcf.genotype(0, 0), gzipped.genotype(0, 0))
        assertEquals(vcf.genotype(1, 0), gzipped.genotype(1, 0))
    }

    @Test
    fun aGzippedHapMapIsNotReadAsVcf() = runTest {
        val text = hapMapText(listOf("B73", "Mo17"), listOf("AC"))
        val file = ByteGenotypeFile("mdp.hmp.txt.gz", storedGzip(text.encodeToByteArray()))

        val table = file.readGenotypes()

        assertEquals(listOf("B73", "Mo17"), table.taxa.names)
        assertEquals(1, table.numberOfSites)
        assertEquals(Nucleotide.genotype('A'), table.genotype(0, 0))
        assertEquals(Nucleotide.genotype('C'), table.genotype(1, 0))
    }

    @Test
    fun aTxtNameStillReadsHapMap() = runTest {
        val text = hapMapText(listOf("a"), listOf("A"))

        val table = ByteGenotypeFile("notes.txt", text.encodeToByteArray()).readGenotypes()

        assertEquals(listOf("a"), table.taxa.names)
        assertEquals(Nucleotide.genotype('A'), table.genotype(0, 0))
    }

    @Test
    fun aVcfThatIsNotVcfFailsWithTheReaderMessage() = runTest {
        val file = ByteGenotypeFile("notes.vcf", "not a genotype file\n".encodeToByteArray())

        val failure = assertFailsWith<IllegalArgumentException> { file.readGenotypes() }

        assertTrue(failure.message!!.contains("VCF"))
    }

    @Test
    fun aVcfReportsProgressWhileItIsRead() = runTest {
        // Vcf.read reports every 256 sites, and once at the end with the real total.
        val file = ByteGenotypeFile("wide.vcf", vcfText(listOf("a", "b"), sites = 300).encodeToByteArray())
        val reports = mutableListOf<Progress>()

        val table = file.readGenotypes(ProgressReporter { reports.add(it) })

        assertEquals(300, table.numberOfSites)
        assertTrue(reports.any { it.message == "Reading VCF" && it.completed in 1 until table.numberOfSites })
        assertEquals(Progress(300, 300, "Reading VCF"), reports.last())
    }

    @Test
    fun indelsAndMissingCallsSurviveTheRoundTripThroughBytes() = runTest {
        val file = TextGenotypeFile("odd.hmp.txt", hapMapText(listOf("a", "b", "c"), listOf("N+-")))

        val table = file.readGenotypes()

        assertEquals(Genotype.UNKNOWN, table.genotype(0, 0))
        assertEquals(Nucleotide.genotype('+'), table.genotype(1, 0))
        assertEquals(Nucleotide.genotype('-'), table.genotype(2, 0))
    }
}

/** One site when [calls] is set, or [sites] homozygous `0/0` sites. Stored gzip of this is a plain member, not BGZF. */
private fun vcfText(taxa: List<String>, calls: List<String>): String = vcfText(taxa, calls, sites = null)

private fun vcfText(taxa: List<String>, sites: Int): String = vcfText(taxa, emptyList(), sites)

private fun vcfText(taxa: List<String>, calls: List<String>, sites: Int?): String {
    val header = (Vcf.FIXED_COLUMNS + taxa).joinToString("\t")
    val records = if (sites == null) {
        listOf(vcfSite(taxa, calls, "s1", "10"))
    } else {
        (1..sites).map { site -> vcfSite(taxa, List(taxa.size) { "0/0" }, "S$site", "$site") }
    }
    return (listOf(header) + records).joinToString("\n") + "\n"
}

private fun vcfSite(taxa: List<String>, calls: List<String>, snpId: String, position: String): String {
    require(calls.size == taxa.size) { "${calls.size} calls for ${taxa.size} taxa" }
    return (listOf("1", position, snpId, "A", "C", ".", "PASS", ".", "GT") + calls).joinToString("\t")
}

/**
 * A one-member gzip file whose deflate body is a single stored block. No compressor: the dispatch tests only need
 * bytes that [net.maizegenetics.tassel.io.decompressing] accepts.
 */
internal fun storedGzip(plain: ByteArray): ByteArray {
    require(plain.size <= 0xFFFF) { "a stored block holds at most 65535 bytes" }
    val stored = ByteArray(5 + plain.size)
    stored[0] = 1
    stored[1] = (plain.size and 0xFF).toByte()
    stored[2] = ((plain.size ushr 8) and 0xFF).toByte()
    val complement = plain.size xor 0xFFFF
    stored[3] = (complement and 0xFF).toByte()
    stored[4] = ((complement ushr 8) and 0xFF).toByte()
    plain.copyInto(stored, 5)
    val header = byteArrayOf(0x1f, 0x8b.toByte(), 8, 0, 0, 0, 0, 0, 0, 0)
    val out = ByteArray(header.size + stored.size + 8)
    header.copyInto(out)
    stored.copyInto(out, header.size)
    val trailer = header.size + stored.size
    putU32(out, trailer, gzipCrc32(plain))
    putU32(out, trailer + 4, plain.size)
    return out
}

/** CRC-32 from RFC 1952: reflected polynomial `0xEDB88320`, initial 0, final xor all-ones. */
private fun gzipCrc32(data: ByteArray): Int {
    var crc = 0.inv()
    for (byte in data) {
        crc = crc xor (byte.toInt() and 0xFF)
        repeat(8) {
            crc = if (crc and 1 != 0) (crc ushr 1) xor 0xEDB88320.toInt() else crc ushr 1
        }
    }
    return crc.inv()
}

private fun putU32(dest: ByteArray, at: Int, value: Int) {
    dest[at] = (value and 0xFF).toByte()
    dest[at + 1] = ((value ushr 8) and 0xFF).toByte()
    dest[at + 2] = ((value ushr 16) and 0xFF).toByte()
    dest[at + 3] = ((value ushr 24) and 0xFF).toByte()
}
