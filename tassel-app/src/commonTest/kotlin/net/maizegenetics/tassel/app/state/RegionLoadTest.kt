package net.maizegenetics.tassel.app.state

import kotlinx.coroutines.CompletableDeferred
import kotlinx.coroutines.CoroutineScope
import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.app.data.RegionSource
import net.maizegenetics.tassel.app.data.storedGzip
import net.maizegenetics.tassel.app.testing.ByteGenotypeFile
import net.maizegenetics.tassel.app.testing.FakeGenotypeFiles
import net.maizegenetics.tassel.app.testing.NamedGenotypeFile
import net.maizegenetics.tassel.app.testing.RecordingResultFiles
import net.maizegenetics.tassel.app.testing.testEngine
import net.maizegenetics.tassel.app.testing.testTable
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.IndexedVcf
import net.maizegenetics.tassel.io.Vcf
import net.maizegenetics.tassel.io.VcfRegion
import net.maizegenetics.tassel.ops.ProgressReporter
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertIs
import kotlin.test.assertNull
import kotlin.test.assertTrue

/** The region step: pairing failures, validation before a read, and cancel. */
class RegionLoadTest {
    @Test
    fun aGzippedVcfWithoutAnIndexLoadsWhole() = runTest {
        val header = (Vcf.FIXED_COLUMNS + listOf("B73")).joinToString("\t")
        val site = listOf("1", "1", "s1", "A", "C", ".", "PASS", ".", "GT", "0/0").joinToString("\t")
        val gzipped = storedGzip("$header\n$site\n".encodeToByteArray())
        val files = FakeGenotypeFiles(chosen = ByteGenotypeFile("sample.vcf.gz", gzipped))
        val viewModel = AppViewModel(files, testEngine(), RecordingResultFiles(), this)

        viewModel.handle(AppEvent.OpenFile)
        testScheduler.advanceUntilIdle()

        assertEquals(Activity.Idle, viewModel.state.value.activity)
        assertEquals("sample.vcf.gz", viewModel.state.value.selected?.name)
        assertEquals(1, viewModel.state.value.selected?.table?.numberOfSites)
        assertNull(viewModel.state.value.region)
    }

    @Test
    fun anUnknownContigStaysOnThePromptAndDoesNotRead() = runTest {
        val source = fakeSource()
        val viewModel = viewModel(this, source)

        viewModel.handle(AppEvent.OpenFile)
        testScheduler.advanceUntilIdle()
        viewModel.handle(AppEvent.SetRegion("chrUn:1-1000"))

        val prompt = viewModel.state.value.region
        assertEquals("chr3B", prompt?.contigs?.single()?.name)
        assertEquals(851_934_019L, prompt?.contigs?.single()?.length)
        assertTrue(prompt?.message.orEmpty().contains("chrUn"))
        assertEquals(0, source.reads)

        viewModel.handle(AppEvent.ConfirmRegion)
        testScheduler.advanceUntilIdle()

        assertEquals(0, source.reads)
        assertTrue(viewModel.state.value.datasets.isEmpty())
        assertEquals("chrUn:1-1000", viewModel.state.value.region?.text)
        assertEquals(Activity.Idle, viewModel.state.value.activity)
    }

    @Test
    fun loadingAWholeContigAddsThatRegionAndClosesTheVcf() = runTest {
        val source = fakeSource { region ->
            assertEquals(VcfRegion("chr3B", null, null), region)
            testTable(2, 3)
        }
        val viewModel = viewModel(this, source)

        viewModel.handle(AppEvent.OpenFile)
        testScheduler.advanceUntilIdle()
        viewModel.handle(AppEvent.UseContig("chr3B"))
        assertNull(viewModel.state.value.region?.message)
        viewModel.handle(AppEvent.ConfirmRegion)
        testScheduler.advanceUntilIdle()

        assertEquals(listOf("wheat.vcf.gz [chr3B]"), viewModel.state.value.datasets.map { it.name })
        assertEquals("2 taxa x 3 sites", viewModel.state.value.selected?.dimensions)
        assertEquals(Activity.Idle, viewModel.state.value.activity)
        assertNull(viewModel.state.value.region)
        assertTrue(source.closed)
        assertEquals(1, source.reads)
    }

    @Test
    fun cancelStopsARegionLoad() = runTest {
        val gate = CompletableDeferred<Unit>()
        val source = fakeSource {
            gate.await()
            testTable(1, 1)
        }
        val viewModel = viewModel(this, source)

        viewModel.handle(AppEvent.OpenFile)
        testScheduler.advanceUntilIdle()
        viewModel.handle(AppEvent.UseContig("chr3B"))
        viewModel.handle(AppEvent.ConfirmRegion)

        assertIs<Activity.Loading>(viewModel.state.value.activity)
        assertEquals(1, source.reads)
        viewModel.handle(AppEvent.Cancel)
        testScheduler.advanceUntilIdle()

        assertEquals(Activity.Idle, viewModel.state.value.activity)
        assertTrue(viewModel.state.value.datasets.isEmpty())
        assertNull(viewModel.state.value.region)
        assertTrue(source.closed)
    }

    @Test
    fun cancellingThePromptReadsNothing() = runTest {
        val source = fakeSource()
        val viewModel = viewModel(this, source)

        viewModel.handle(AppEvent.OpenFile)
        testScheduler.advanceUntilIdle()
        viewModel.handle(AppEvent.CancelRegion)

        assertNull(viewModel.state.value.region)
        assertEquals(Activity.Idle, viewModel.state.value.activity)
        assertTrue(viewModel.state.value.datasets.isEmpty())
        assertEquals(0, source.reads)
        assertTrue(source.closed)
    }

    private fun viewModel(scope: CoroutineScope, source: FakeRegionSource): AppViewModel {
        val index = NamedGenotypeFile("wheat.vcf.gz.csi")
        val files = FakeGenotypeFiles(chosen = NamedGenotypeFile("wheat.vcf.gz", index))
        return AppViewModel(
            files,
            testEngine(),
            RecordingResultFiles(),
            scope,
            openIndexed = { _, _ -> source },
        )
    }

    private fun fakeSource(onRead: suspend (VcfRegion) -> GenotypeTable = { testTable(1, 1) }): FakeRegionSource =
        FakeRegionSource(
            contigs = listOf(IndexedVcf.Contig("chr3B", 851_934_019L)),
            onRead = onRead,
        )
}

/** An indexed VCF that never touches BGZF. [reads] counts [read] calls. */
private class FakeRegionSource(
    override val fileName: String = "wheat.vcf.gz",
    override val contigs: List<IndexedVcf.Contig>,
    private val onRead: suspend (VcfRegion) -> GenotypeTable,
) : RegionSource {
    var reads: Int = 0
        private set

    var closed: Boolean = false
        private set

    override suspend fun read(region: VcfRegion, progress: ProgressReporter): GenotypeTable {
        reads++
        return onRead(region)
    }

    override fun close() {
        closed = true
    }
}
