package net.maizegenetics.tassel.app.data

import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.app.testing.FakeGenotypeFiles
import net.maizegenetics.tassel.app.testing.SMALL_HAPMAP
import net.maizegenetics.tassel.app.testing.TextGenotypeFile
import net.maizegenetics.tassel.app.testing.hapMapText
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.Nucleotide
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
    fun indelsAndMissingCallsSurviveTheRoundTripThroughBytes() = runTest {
        val file = TextGenotypeFile("odd.hmp.txt", hapMapText(listOf("a", "b", "c"), listOf("N+-")))

        val table = file.readGenotypes()

        assertEquals(Genotype.UNKNOWN, table.genotype(0, 0))
        assertEquals(Nucleotide.genotype('+'), table.genotype(1, 0))
        assertEquals(Nucleotide.genotype('-'), table.genotype(2, 0))
    }
}
