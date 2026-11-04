package net.maizegenetics.tassel.app.state

import kotlinx.coroutines.CompletableDeferred
import kotlinx.coroutines.test.TestScope
import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.app.testing.FakeGenotypeFiles
import net.maizegenetics.tassel.app.testing.RecordingResultFiles
import net.maizegenetics.tassel.app.testing.SMALL_HAPMAP
import net.maizegenetics.tassel.app.testing.TextGenotypeFile
import net.maizegenetics.tassel.app.testing.testEngine
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.ops.OperationRegistry
import net.maizegenetics.tassel.ops.TasselOperations
import net.maizegenetics.tassel.ops.operation
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.table.ResultTable
import net.maizegenetics.tassel.table.resultTable
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertIs
import kotlin.test.assertNotNull
import kotlin.test.assertNull
import kotlin.test.assertTrue

/**
 * The state machine around running an operation: what the UI draws while it runs, when it is cancelled, when its
 * parameters are wrong, and what it leaves in the tree.
 *
 * Runs go through a real `LocalEngine`, so these also check the engine contract the app depends on: text
 * parameters, and a progress report that doubles as a cancellation check (ADR-0004).
 */
class AnalysisRunTest {
    @Test
    fun theOperationsTheFormsAreGeneratedFromComeFromTheEngine() = runTest {
        val viewModel = loaded(this)

        assertEquals(
            TasselOperations.registry.operations.map { it.id },
            viewModel.state.value.operations.map { it.id },
        )
    }

    @Test
    fun openingAnAnalysisFillsItsFormWithTheSpecsDefaults() = runTest {
        val viewModel = loaded(this)

        viewModel.handle(AppEvent.ShowForm("kinship"))

        val form = assertNotNull(viewModel.state.value.form)
        assertEquals("kinship", form.operationId)
        assertEquals(TasselOperations.registry["kinship"]!!.parameters.map { it.name }, form.fields.map { it.name })
        assertEquals("Centered_IBS", form.values["method"])
    }

    @Test
    fun anUnknownOperationOpensNoForm() = runTest {
        val viewModel = loaded(this)

        viewModel.handle(AppEvent.ShowForm("no-such-operation"))

        assertNull(viewModel.state.value.form)
    }

    @Test
    fun runningAnAnalysisAddsItsResultUnderTheDatasetAndShowsIt() = runTest {
        val viewModel = loaded(this)

        viewModel.handle(AppEvent.ShowForm("kinship"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        val state = viewModel.state.value
        assertEquals(Activity.Idle, state.activity)
        // The form closes on success, so the pane the result opened in is what the user sees.
        assertNull(state.form)
        assertEquals(1, state.results.size)
        assertEquals("kinship", state.results[0].operationId)
        assertEquals(state.datasets[0].id, state.results[0].datasetId)
        assertEquals(Showing.Result(state.results[0].id, 0), state.showing)
        assertEquals("Centered_IBS", state.selectedPart?.title)
    }

    @Test
    fun everyRegisteredOperationRunsFromTheUiWithItsDefaults() = runTest {
        // Nothing here names an operation: the five come from the engine, as the form generator takes them.
        val viewModel = loaded(this)
        val operations = viewModel.state.value.operations.map { it.id }
        assertEquals(5, operations.size)

        for (operation in operations) {
            viewModel.handle(AppEvent.Select(viewModel.state.value.datasets[0].id))
            viewModel.handle(AppEvent.ShowForm(operation))
            viewModel.handle(AppEvent.Run)
            testScheduler.advanceUntilIdle()

            assertEquals(Activity.Idle, viewModel.state.value.activity, "$operation left the app busy")
            assertNull(viewModel.state.value.form, "$operation did not finish")
        }
        // summary, kinship, and pca leave result nodes; the two filters leave genotype tables in the tree.
        assertEquals(listOf("summary", "kinship", "pca"), viewModel.state.value.results.map { it.operationId })
        assertEquals(
            listOf("tutorial.hmp.txt", "tutorial.hmp.txt + filter-sites", "tutorial.hmp.txt + filter-taxa"),
            viewModel.state.value.datasets.map { it.name },
        )
    }

    @Test
    fun progressFromTheOperationIsWhatTheBannerShows() = runTest {
        val gate = CompletableDeferred<Unit>()
        val viewModel = loaded(this, OperationRegistry(gatedOperation(gate)))

        viewModel.handle(AppEvent.ShowForm("gated"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        val running = assertIs<Activity.Running>(viewModel.state.value.activity)
        assertEquals("gated", running.operationId)
        assertEquals(1L, running.progress.completed)
        assertEquals(2L, running.progress.total)
        assertEquals("halfway", running.progress.message)
        assertTrue(viewModel.state.value.busy)

        gate.complete(Unit)
        testScheduler.advanceUntilIdle()
        assertEquals(Activity.Idle, viewModel.state.value.activity)
        assertEquals(1, viewModel.state.value.results.size)
    }

    @Test
    fun cancellingMidRunStopsTheOperationAndLeavesNoResult() = runTest {
        val gate = CompletableDeferred<Unit>()
        val operation = gatedOperation(gate)
        val viewModel = loaded(this, OperationRegistry(operation))

        viewModel.handle(AppEvent.ShowForm("gated"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()
        assertIs<Activity.Running>(viewModel.state.value.activity)

        // There is no cancel call on the engine; cancelling the coroutine driving the run is the whole of it.
        viewModel.handle(AppEvent.Cancel)
        gate.complete(Unit)
        testScheduler.advanceUntilIdle()

        assertEquals(Activity.Idle, viewModel.state.value.activity)
        assertEquals(emptyList(), viewModel.state.value.results)
        // The operation stopped at its next progress report rather than running to the end.
        assertTrue(viewModel.state.value.form != null, "the form stays open so the run can be started again")
    }

    @Test
    fun aSecondRunReplacesTheFirstRatherThanRacingIt() = runTest {
        val gate = CompletableDeferred<Unit>()
        val viewModel = loaded(this, OperationRegistry(gatedOperation(gate)))

        viewModel.handle(AppEvent.ShowForm("gated"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()
        viewModel.handle(AppEvent.Run)
        gate.complete(Unit)
        testScheduler.advanceUntilIdle()

        // The cancelled first run does not clear the second one's banner, and only one result is produced.
        assertEquals(Activity.Idle, viewModel.state.value.activity)
        assertEquals(1, viewModel.state.value.results.size)
    }

    @Test
    fun aValueTheEngineRejectsComesBackOnTheFieldItBelongsTo() = runTest {
        val viewModel = loaded(this)

        viewModel.handle(AppEvent.ShowForm("kinship"))
        viewModel.handle(AppEvent.SetParameter("maxAlleles", "99"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        val form = assertNotNull(viewModel.state.value.form)
        assertEquals(setOf("maxAlleles"), form.problems.keys)
        assertEquals(Activity.Idle, viewModel.state.value.activity)
        assertEquals(emptyList(), viewModel.state.value.results)
    }

    @Test
    fun anOperationThatThrowsIsReportedAndDismissible() = runTest {
        val viewModel = loaded(this)

        viewModel.handle(AppEvent.ShowForm("pca"))
        // limitBy is declared for the legacy translator; only number_of_components is ported (tassel-6#19).
        viewModel.handle(AppEvent.SetParameter("limitBy", "min_eigenvalue"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        val failed = assertIs<Activity.Failed>(viewModel.state.value.activity)
        assertEquals("pca", failed.name)
        assertTrue(failed.message.isNotEmpty())

        viewModel.handle(AppEvent.DismissFailure)
        assertEquals(Activity.Idle, viewModel.state.value.activity)
    }

    @Test
    fun runningWithNoDatasetDoesNothing() = runTest {
        val viewModel = AppViewModel(FakeGenotypeFiles(), testEngine(), RecordingResultFiles(), this)
        testScheduler.advanceUntilIdle()

        viewModel.handle(AppEvent.ShowForm("kinship"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        assertEquals(Activity.Idle, viewModel.state.value.activity)
        assertEquals(emptyList(), viewModel.state.value.results)
    }

    @Test
    fun closingADatasetTakesItsResultsWithIt() = runTest {
        val viewModel = loaded(this)
        viewModel.handle(AppEvent.ShowForm("kinship"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        viewModel.handle(AppEvent.Close(viewModel.state.value.datasets[0].id))

        assertEquals(emptyList(), viewModel.state.value.results)
        assertEquals(Showing.Nothing, viewModel.state.value.showing)
    }

    @Test
    fun closingAResultGoesBackToTheDatasetItCameFrom() = runTest {
        val viewModel = loaded(this)
        viewModel.handle(AppEvent.ShowForm("kinship"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()
        val dataset = viewModel.state.value.datasets[0].id

        viewModel.handle(AppEvent.CloseResult(viewModel.state.value.results[0].id))

        assertEquals(emptyList(), viewModel.state.value.results)
        assertEquals(Showing.Dataset(dataset), viewModel.state.value.showing)
    }

    /** A view model with one small dataset loaded and shown. */
    private fun loaded(scope: TestScope, registry: OperationRegistry = TasselOperations.registry): AppViewModel {
        val files = FakeGenotypeFiles(tutorial = TextGenotypeFile("tutorial.hmp.txt", SMALL_HAPMAP))
        val viewModel = AppViewModel(files, testEngine(registry), RecordingResultFiles(), scope)
        viewModel.handle(AppEvent.OpenTutorial)
        scope.testScheduler.advanceUntilIdle()
        return viewModel
    }

    /**
     * An operation that reports progress, waits for [gate], then reports again. The pause is where a test can
     * look at the running state, and where a cancelled run stops: `LocalEngine` throws out of the second report.
     */
    private fun gatedOperation(gate: CompletableDeferred<Unit>) =
        operation<GenotypeTable, ResultTable>(id = "gated", help = "Waits, for the tests") {
            run { input, progress ->
                progress.report(1, 2, "halfway")
                gate.await()
                progress.report(2, 2, "done")
                resultTable("Gated", listOf("Taxa")) { row(input.numberOfTaxa) }
            }
        }
}
