package net.maizegenetics.tassel.app.state

import kotlinx.coroutines.test.TestScope
import kotlinx.coroutines.test.runTest
import kotlinx.io.Sink
import net.maizegenetics.tassel.app.export.ResultFiles
import net.maizegenetics.tassel.app.testing.FakeGenotypeFiles
import net.maizegenetics.tassel.app.testing.RecordingResultFiles
import net.maizegenetics.tassel.app.testing.SMALL_HAPMAP
import net.maizegenetics.tassel.app.testing.TextGenotypeFile
import net.maizegenetics.tassel.app.testing.testEngine
import net.maizegenetics.tassel.io.HapMap
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/** Exporting a result, which is the last step of every analysis and the one that has to match TASSEL 5. */
class ExportTest {
    @Test
    fun exportingTheShownPaneWritesTheTasselTextForIt() = runTest {
        val saved = RecordingResultFiles()
        val viewModel = loaded(this, saved)
        viewModel.handle(AppEvent.ShowForm("kinship"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        viewModel.handle(AppEvent.ExportResult)
        testScheduler.advanceUntilIdle()

        assertEquals(listOf("tutorial.kinship.txt"), saved.names)
        val lines = saved.contents.single().lineSequence().toList()
        // TASSEL 5's square-matrix text: the annotations sorted by key, the taxa count, then a row per taxon.
        assertTrue(lines[0].startsWith("##Centered_IBS.SumPk="), lines[0])
        assertEquals("##Matrix_Type=Centered_IBS", lines[1])
        assertEquals("3", lines[2])
        assertTrue(lines[3].startsWith("B73\t"), lines[3])
        val text = saved.contents.single()
        assertTrue(text.endsWith("\n"))
    }

    @Test
    fun eachPaneOfAMultiPartResultExportsItsOwnFile() = runTest {
        val saved = RecordingResultFiles()
        val viewModel = loaded(this, saved)
        viewModel.handle(AppEvent.ShowForm("summary"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()
        val result = viewModel.state.value.results.single()

        result.parts.indices.forEach { part ->
            viewModel.handle(AppEvent.SelectResult(result.id, part))
            viewModel.handle(AppEvent.ExportResult)
            testScheduler.advanceUntilIdle()
        }

        assertEquals(
            listOf(
                "tutorial.summary.overall_summary.txt",
                "tutorial.summary.allele_summary.txt",
                "tutorial.summary.site_summary.txt",
                "tutorial.summary.taxa_summary.txt",
            ),
            saved.names,
        )
        // Each file is the core writer's own text for that table, header included.
        result.parts.forEachIndexed { index, part ->
            assertEquals(part.table.columnNames.joinToString("\t"), saved.contents[index].lineSequence().first())
        }
    }

    @Test
    fun aFilteredGenotypeTableIsExportedAsHapMap() = runTest {
        val saved = RecordingResultFiles()
        val viewModel = loaded(this, saved)
        viewModel.handle(AppEvent.ShowForm("filter-sites"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()
        val filtered = viewModel.state.value.datasets.last()

        viewModel.handle(AppEvent.ExportDataset(filtered.id))
        testScheduler.advanceUntilIdle()

        assertEquals(listOf("tutorial.hmp.txt + filter-sites.hmp.txt"), saved.names)
        assertEquals(HapMap.FIXED_COLUMNS.first(), saved.contents.single().substringBefore('\t'))
    }

    @Test
    fun exportingWithNothingShownWritesNothing() = runTest {
        val saved = RecordingResultFiles()
        val viewModel = loaded(this, saved)

        viewModel.handle(AppEvent.ExportResult)
        viewModel.handle(AppEvent.ExportDataset(-1))
        testScheduler.advanceUntilIdle()

        assertEquals(emptyList(), saved.names)
    }

    @Test
    fun aSaveThatFailsIsReportedRatherThanLost() = runTest {
        val viewModel = AppViewModel(
            FakeGenotypeFiles(tutorial = TextGenotypeFile("tutorial.hmp.txt", SMALL_HAPMAP)),
            testEngine(),
            FailingResultFiles(),
            this,
        )
        viewModel.handle(AppEvent.OpenTutorial)
        testScheduler.advanceUntilIdle()

        viewModel.handle(AppEvent.ExportDataset(viewModel.state.value.datasets.single().id))
        testScheduler.advanceUntilIdle()

        val failed = viewModel.state.value.activity
        assertTrue(failed is Activity.Failed, "$failed")
        assertEquals("tutorial.hmp.txt", failed.name)
    }

    @Test
    fun theFileNameFollowsTheClisOutputNaming() {
        assertEquals("mdp_genotype.kinship", baseName("mdp_genotype.hmp.txt", "kinship"))
        assertEquals("result.pca", baseName("", "pca"))
    }

    private fun loaded(scope: TestScope, saved: RecordingResultFiles): AppViewModel {
        val files = FakeGenotypeFiles(tutorial = TextGenotypeFile("tutorial.hmp.txt", SMALL_HAPMAP))
        val viewModel = AppViewModel(files, testEngine(), saved, scope)
        viewModel.handle(AppEvent.OpenTutorial)
        scope.testScheduler.advanceUntilIdle()
        return viewModel
    }
}

/** A [ResultFiles] whose save always fails, as a full disk would. */
private class FailingResultFiles : ResultFiles {
    override suspend fun save(fileName: String, write: (Sink) -> Unit): Boolean =
        throw IllegalStateException("no room on device")
}
