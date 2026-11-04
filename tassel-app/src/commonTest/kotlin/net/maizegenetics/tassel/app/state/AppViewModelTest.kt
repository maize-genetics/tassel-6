package net.maizegenetics.tassel.app.state

import kotlinx.coroutines.CompletableDeferred
import kotlinx.coroutines.test.TestScope
import kotlinx.coroutines.test.advanceUntilIdle
import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.app.testing.FakeGenotypeFiles
import net.maizegenetics.tassel.app.testing.SMALL_HAPMAP
import net.maizegenetics.tassel.app.testing.TextGenotypeFile
import net.maizegenetics.tassel.app.testing.hapMapText
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertIs
import kotlin.test.assertNull
import kotlin.test.assertTrue

/** The state machine behind the data tree: what the UI draws after each event. */
class AppViewModelTest {
    @Test
    fun theTreeStartsEmpty() = runTest {
        val state = AppViewModel(FakeGenotypeFiles(), this).state.value

        assertEquals(emptyList(), state.datasets)
        assertNull(state.selectedId)
        assertNull(state.selected)
        assertEquals(Activity.Idle, state.activity)
    }

    @Test
    fun loadingTheTutorialDatasetAddsItAndSelectsIt() = runTest {
        val files = FakeGenotypeFiles(tutorial = TextGenotypeFile("mdp_genotype.hmp.txt", SMALL_HAPMAP))
        val viewModel = AppViewModel(files, this)

        viewModel.handle(AppEvent.OpenTutorial)
        advanceUntilIdle()

        val state = viewModel.state.value
        assertEquals(1, files.tutorialRequests)
        assertEquals(Activity.Idle, state.activity)
        assertEquals("mdp_genotype.hmp.txt", state.selected?.name)
        assertEquals("3 taxa x 2 sites", state.selected?.dimensions)
        assertEquals(3, state.selected?.table?.numberOfTaxa)
    }

    @Test
    fun loadingShowsProgressUntilTheReaderFinishes() = runTest {
        val gate = CompletableDeferred<Unit>()
        val file = TextGenotypeFile("big.hmp.txt", SMALL_HAPMAP) { gate.await() }
        val viewModel = AppViewModel(FakeGenotypeFiles(tutorial = file), this)

        viewModel.handle(AppEvent.OpenTutorial)

        val loading = viewModel.state.value.activity
        assertIs<Activity.Loading>(loading)
        assertEquals("big.hmp.txt", loading.name)
        assertTrue(viewModel.state.value.datasets.isEmpty())

        gate.complete(Unit)
        advanceUntilIdle()
        assertEquals(Activity.Idle, viewModel.state.value.activity)
        assertEquals(1, viewModel.state.value.datasets.size)
    }

    @Test
    fun cancellingTheChooserChangesNothing() = runTest {
        val files = FakeGenotypeFiles(chosen = null)
        val viewModel = AppViewModel(files, this)

        viewModel.handle(AppEvent.OpenFile)
        advanceUntilIdle()

        assertEquals(1, files.chooseRequests)
        assertEquals(AppState(), viewModel.state.value)
    }

    @Test
    fun aFileThatCannotBeReadIsReportedAndDismissible() = runTest {
        val files = FakeGenotypeFiles(chosen = TextGenotypeFile("notes.txt", "not a genotype file\n"))
        val viewModel = AppViewModel(files, this)

        viewModel.handle(AppEvent.OpenFile)
        advanceUntilIdle()

        val failed = viewModel.state.value.activity
        assertIs<Activity.Failed>(failed)
        assertEquals("notes.txt", failed.name)
        assertTrue(failed.message.isNotEmpty())
        assertTrue(viewModel.state.value.datasets.isEmpty())

        viewModel.handle(AppEvent.DismissFailure)
        assertEquals(Activity.Idle, viewModel.state.value.activity)
    }

    @Test
    fun aSecondLoadCancelsTheFirst() = runTest {
        // The browser has one thread; two readers racing for it would make both of them slow.
        val gate = CompletableDeferred<Unit>()
        val slow = TextGenotypeFile("slow.hmp.txt", SMALL_HAPMAP) { gate.await() }
        val quick = TextGenotypeFile("quick.hmp.txt", hapMapText(listOf("a"), listOf("A")))
        val viewModel = AppViewModel(FakeGenotypeFiles(tutorial = slow, chosen = quick), this)

        viewModel.handle(AppEvent.OpenTutorial)
        viewModel.handle(AppEvent.OpenFile)
        advanceUntilIdle()

        assertEquals(listOf("quick.hmp.txt"), viewModel.state.value.datasets.map { it.name })
        assertEquals(0, slow.opened)
    }

    @Test
    fun selectingADatasetChangesWhichOneTheGridShows() = runTest {
        val viewModel = loadedTwice(this)

        viewModel.handle(AppEvent.Select(viewModel.state.value.datasets.first().id))

        assertEquals("first.hmp.txt", viewModel.state.value.selected?.name)
    }

    @Test
    fun closingTheSelectedDatasetFallsBackToTheLastOneLeft() = runTest {
        val viewModel = loadedTwice(this)
        val second = viewModel.state.value.datasets.last()

        viewModel.handle(AppEvent.Close(second.id))

        assertEquals(listOf("first.hmp.txt"), viewModel.state.value.datasets.map { it.name })
        assertEquals("first.hmp.txt", viewModel.state.value.selected?.name)
    }

    @Test
    fun closingTheLastDatasetLeavesNothingSelected() = runTest {
        val viewModel = loadedTwice(this)
        viewModel.state.value.datasets.forEach { viewModel.handle(AppEvent.Close(it.id)) }

        assertEquals(emptyList(), viewModel.state.value.datasets)
        assertNull(viewModel.state.value.selected)
    }

    @Test
    fun closingAnUnselectedDatasetKeepsTheSelection() = runTest {
        val viewModel = loadedTwice(this)
        val first = viewModel.state.value.datasets.first()

        viewModel.handle(AppEvent.Close(first.id))

        assertEquals("second.hmp.txt", viewModel.state.value.selected?.name)
    }

    /** A view model with two datasets loaded, the second one selected. */
    private fun loadedTwice(scope: TestScope): AppViewModel {
        val first = TextGenotypeFile("first.hmp.txt", SMALL_HAPMAP)
        val second = TextGenotypeFile("second.hmp.txt", SMALL_HAPMAP)
        val viewModel = AppViewModel(FakeGenotypeFiles(tutorial = first, chosen = second), scope)
        viewModel.handle(AppEvent.OpenTutorial)
        scope.advanceUntilIdle()
        viewModel.handle(AppEvent.OpenFile)
        scope.advanceUntilIdle()
        return viewModel
    }
}
