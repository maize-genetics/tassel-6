package net.maizegenetics.tassel.app.state

import kotlinx.coroutines.test.TestScope
import kotlinx.coroutines.test.runTest
import net.maizegenetics.tassel.app.analysis.AnalysisEngine
import net.maizegenetics.tassel.app.export.GenotypeExportFormat
import net.maizegenetics.tassel.app.export.phenotypeFileName
import net.maizegenetics.tassel.app.testing.FakeGenotypeFiles
import net.maizegenetics.tassel.app.testing.RecordingResultFiles
import net.maizegenetics.tassel.app.testing.ReleaseRecordingEngine
import net.maizegenetics.tassel.app.testing.SMALL_HAPMAP
import net.maizegenetics.tassel.app.testing.TextGenotypeFile
import net.maizegenetics.tassel.app.testing.phenotypeText
import net.maizegenetics.tassel.app.testing.testEngine
import net.maizegenetics.tassel.dataset.DataSet
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertIs
import kotlin.test.assertNotNull
import kotlin.test.assertNull
import kotlin.test.assertTrue

/** Phenotype datasets in the tree (#167): what opens as one, what it offers, what joins take, and what it exports. */
class PhenotypeDatasetTest {
    private val traits = TextGenotypeFile("traits.txt", phenotypeText(mapOf("EarHT" to listOf("64.75", "NA", "30"))))
    private val extra = TextGenotypeFile("extra.txt", phenotypeText(mapOf("dpoll" to listOf("60", "61", "62"))))
    private val tutorialGenotype = TextGenotypeFile("mdp_genotype.hmp.txt", SMALL_HAPMAP)

    @Test
    fun aDatasetKnowsItsTypeAndDimensions() = runTest {
        val viewModel = opened(this, traits)

        val dataset = viewModel.state.value.datasets.single()
        assertIs<DatasetData.Phenotypes>(dataset.data)
        assertEquals("PhenotypeTable", dataset.inputType)
        assertEquals("3 observations x 2 attributes", dataset.dimensions)
        assertNotNull(dataset.phenotype)
        assertNull(dataset.genotypes)
    }

    @Test
    fun aGenotypeDatasetKeepsItsTypeAndDimensions() = runTest {
        val viewModel = opened(this, tutorialGenotype)

        val dataset = viewModel.state.value.datasets.single()
        assertIs<DatasetData.Genotypes>(dataset.data)
        assertEquals("GenotypeTable", dataset.inputType)
        assertEquals("3 taxa x 2 sites", dataset.dimensions)
        assertNull(dataset.phenotype)
    }

    @Test
    fun theTutorialOpensTheGenotypeAndBothTraitFilesWithTheGenotypeShown() = runTest {
        val files = FakeGenotypeFiles(
            tutorial = tutorialGenotype,
            tutorialPhenotypes = listOf(traits, extra),
        )
        val viewModel = AppViewModel(files, testEngine(), RecordingResultFiles(), this)

        viewModel.handle(AppEvent.OpenTutorial)
        testScheduler.advanceUntilIdle()

        val state = viewModel.state.value
        assertEquals(listOf("mdp_genotype.hmp.txt", "traits.txt", "extra.txt"), state.datasets.map { it.name })
        assertEquals(listOf("GenotypeTable", "PhenotypeTable", "PhenotypeTable"), state.datasets.map { it.inputType })
        assertEquals("mdp_genotype.hmp.txt", state.selected?.name)
        assertEquals(Activity.Idle, state.activity)
    }

    @Test
    fun aPickedFileIsSniffedNotJudgedByItsName() = runTest {
        // A phenotype called `.hmp.txt` is still a phenotype, and a HapMap called `.txt` is still a genotype.
        val phenotype = TextGenotypeFile("looks_like.hmp.txt", phenotypeText(mapOf("EarHT" to listOf("1", "2", "3"))))
        val genotype = TextGenotypeFile("plain.txt", SMALL_HAPMAP)

        assertEquals("PhenotypeTable", opened(this, phenotype).state.value.datasets.single().inputType)
        assertEquals("GenotypeTable", opened(this, genotype).state.value.datasets.single().inputType)
    }

    @Test
    fun aBadPhenotypeIsReportedWithTheReadersMessage() = runTest {
        val broken = TextGenotypeFile("broken.txt", "<Phenotype>\ntaxa\tdata\nTaxon\tX\nB73\tnot-a-number\n")

        val viewModel = opened(this, broken)

        val failed = viewModel.state.value.activity
        assertIs<Activity.Failed>(failed)
        assertEquals("broken.txt", failed.name)
        assertTrue(viewModel.state.value.datasets.isEmpty())
    }

    @Test
    fun operationsAreFilteredByTheInputTypeOfTheSelection() = runTest {
        val files = FakeGenotypeFiles(tutorial = tutorialGenotype, tutorialPhenotypes = listOf(traits, extra))
        val viewModel = AppViewModel(files, testEngine(), RecordingResultFiles(), this)
        viewModel.handle(AppEvent.OpenTutorial)
        testScheduler.advanceUntilIdle()
        val (genotype, first, second) = viewModel.state.value.datasets

        viewModel.handle(AppEvent.Select(genotype.id))
        val forGenotype = viewModel.state.value.availableOperations.map { it.id }
        viewModel.handle(AppEvent.Select(first.id))
        val forPhenotype = viewModel.state.value.availableOperations.map { it.id }
        viewModel.handle(AppEvent.ToggleChecked(second.id))
        val forBoth = viewModel.state.value.availableOperations.map { it.id }

        val operations = viewModel.state.value.operations
        assertEquals(operations.filter { it.inputType == "GenotypeTable" }.map { it.id }, forGenotype)
        assertEquals(operations.filter { it.inputType == "PhenotypeTable" }.map { it.id }, forPhenotype)
        assertEquals(operations.filter { it.inputType == "DataSet" }.map { it.id }, forBoth)
        assertTrue(forGenotype.isNotEmpty() && forPhenotype.isNotEmpty() && forBoth.isNotEmpty())
        assertTrue(forGenotype.intersect(forPhenotype.toSet()).isEmpty())
        assertTrue("join-phenotypes" in forBoth && "join-phenotypes" !in forPhenotype)
    }

    @Test
    fun nothingIsOfferedWithNothingSelected() = runTest {
        val viewModel = AppViewModel(FakeGenotypeFiles(), testEngine(), RecordingResultFiles(), this)
        testScheduler.advanceUntilIdle()

        assertEquals(emptyList(), viewModel.state.value.availableOperations)
        assertEquals(emptyList(), operationsFor(viewModel.state.value.operations, emptyList()))
    }

    @Test
    fun tickingBuildsTheSelectionInTickOrderAndStartsFromTheShownDataset() = runTest {
        val viewModel = twoPhenotypes(this)
        val (first, second) = viewModel.state.value.datasets
        viewModel.handle(AppEvent.Select(second.id))

        // `second` is shown, so ticking `first` makes the selection [second, first]: the shown one leads.
        viewModel.handle(AppEvent.ToggleChecked(first.id))
        assertEquals(listOf(second.id, first.id), viewModel.state.value.selection.map { it.id })

        viewModel.handle(AppEvent.ToggleChecked(second.id))
        assertEquals(listOf(first.id), viewModel.state.value.selection.map { it.id })

        // A plain click goes back to a single selection.
        viewModel.handle(AppEvent.ToggleChecked(second.id))
        viewModel.handle(AppEvent.Select(first.id))
        assertEquals(emptyList<Long>(), viewModel.state.value.checked)
        assertEquals(listOf(first.id), viewModel.state.value.selection.map { it.id })
    }

    @Test
    fun aMultiSelectionBecomesADataSetInSelectionOrder() = runTest {
        val viewModel = twoPhenotypes(this)
        val (first, second) = viewModel.state.value.datasets
        viewModel.handle(AppEvent.Select(second.id))
        viewModel.handle(AppEvent.ToggleChecked(first.id))

        val dataSet = dataSetOf(viewModel.state.value.selection)

        assertEquals(DataSet(listOf(second.phenotype!!, first.phenotype!!)), dataSet)
        assertEquals(listOf(second.phenotype, first.phenotype), dataSet.phenotypes())
    }

    @Test
    fun joiningTwoPhenotypesAddsAPhenotypeDatasetAndReleasesTheDataSet() = runTest {
        val engine = ReleaseRecordingEngine()
        val viewModel = twoPhenotypes(this, engine)
        val (first, second) = viewModel.state.value.datasets
        viewModel.handle(AppEvent.Select(first.id))
        viewModel.handle(AppEvent.ToggleChecked(second.id))

        viewModel.handle(AppEvent.ShowForm("join-phenotypes"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        val state = viewModel.state.value
        assertEquals(Activity.Idle, state.activity)
        val joined = state.datasets.last()
        assertEquals("traits + extra + join-phenotypes", joined.name)
        assertEquals(listOf("Taxon", "EarHT", "dpoll"), joined.phenotype!!.attributes.map { it.name })
        assertEquals("3 observations x 3 attributes", joined.dimensions)
        assertEquals(Showing.Dataset(joined.id), state.showing)
        assertEquals(emptyList(), state.checked)
        assertEquals(1, engine.released.size)
    }

    @Test
    fun aPhenotypeOperationRunsOnAPhenotypeDatasetAndItsTableBecomesADataset() = runTest {
        val viewModel = opened(this, traits)
        viewModel.handle(AppEvent.ShowForm("filter-phenotype"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        val state = viewModel.state.value
        assertEquals(Activity.Idle, state.activity)
        assertEquals(2, state.datasets.size)
        assertEquals("traits.txt + filter-phenotype", state.datasets.last().name)
        assertEquals("PhenotypeTable", state.datasets.last().inputType)
    }

    @Test
    fun aSummaryOfAPhenotypeIsAResultNotADataset() = runTest {
        val viewModel = opened(this, traits)
        viewModel.handle(AppEvent.ShowForm("phenotype-summary"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        val state = viewModel.state.value
        assertEquals(1, state.datasets.size)
        assertEquals(listOf("phenotype-summary"), state.results.map { it.operationId })
        assertEquals(3, state.results.single().parts.size)
    }

    @Test
    fun aJoinOfAGenotypeAndAPhenotypeIsAResultOfTheAlignedPhenotype() = runTest {
        val files = FakeGenotypeFiles(tutorial = tutorialGenotype, tutorialPhenotypes = listOf(traits))
        val viewModel = AppViewModel(files, testEngine(), RecordingResultFiles(), this)
        viewModel.handle(AppEvent.OpenTutorial)
        testScheduler.advanceUntilIdle()
        val (genotype, phenotype) = viewModel.state.value.datasets
        viewModel.handle(AppEvent.Select(genotype.id))
        viewModel.handle(AppEvent.ToggleChecked(phenotype.id))

        viewModel.handle(AppEvent.ShowForm("join-genotype-phenotype"))
        viewModel.handle(AppEvent.Run)
        testScheduler.advanceUntilIdle()

        val state = viewModel.state.value
        assertEquals(Activity.Idle, state.activity)
        val result = state.results.single()
        assertEquals(genotype.id, result.datasetId)
        assertEquals("Phenotype", result.parts.single().title)
    }

    @Test
    fun theOpenFormClosesWhenTheSelectionStopsFittingIt() = runTest {
        val files = FakeGenotypeFiles(tutorial = tutorialGenotype, tutorialPhenotypes = listOf(traits))
        val viewModel = AppViewModel(files, testEngine(), RecordingResultFiles(), this)
        viewModel.handle(AppEvent.OpenTutorial)
        testScheduler.advanceUntilIdle()
        val (genotype, phenotype) = viewModel.state.value.datasets
        viewModel.handle(AppEvent.Select(phenotype.id))
        viewModel.handle(AppEvent.ShowForm("phenotype-summary"))

        viewModel.handle(AppEvent.Select(genotype.id))

        assertNull(viewModel.state.value.form)
    }

    @Test
    fun aPhenotypeExportsAsTheCoreWritersBytes() = runTest {
        val saved = RecordingResultFiles()
        val viewModel = opened(this, traits, saved)
        val dataset = viewModel.state.value.datasets.single()

        viewModel.handle(AppEvent.ExportDataset(dataset.id))
        testScheduler.advanceUntilIdle()

        assertEquals(listOf("traits.txt"), saved.names)
        assertEquals(dataset.phenotype!!.toDelimitedString(), saved.contents.single())
        assertTrue(saved.contents.single().startsWith("<Phenotype>\n"))
    }

    @Test
    fun aFormatChosenForAPhenotypeIsIgnored() = runTest {
        val saved = RecordingResultFiles()
        val viewModel = opened(this, traits, saved)
        val dataset = viewModel.state.value.datasets.single()

        viewModel.handle(AppEvent.ExportDataset(dataset.id, GenotypeExportFormat.Bcf))
        testScheduler.advanceUntilIdle()

        assertEquals(listOf("traits.txt"), saved.names)
        assertEquals(dataset.phenotype!!.toDelimitedString(), saved.contents.single())
    }

    @Test
    fun aJoinedPhenotypeExportsUnderADerivedName() = runTest {
        assertEquals("mdp_traits.txt", phenotypeFileName("mdp_traits.txt"))
        assertEquals("mdp_traits.txt", phenotypeFileName("mdp_traits.TXT"))
        assertEquals(
            "mdp_traits + mdp_structure + join-phenotypes.txt",
            phenotypeFileName("mdp_traits + mdp_structure + join-phenotypes"),
        )
    }

    @Test
    fun closingATickedDatasetUnticksIt() = runTest {
        val viewModel = twoPhenotypes(this)
        val (first, second) = viewModel.state.value.datasets
        viewModel.handle(AppEvent.Select(first.id))
        viewModel.handle(AppEvent.ToggleChecked(second.id))

        viewModel.handle(AppEvent.Close(second.id))

        assertEquals(listOf(first.id), viewModel.state.value.checked)
        assertEquals(listOf(first.id), viewModel.state.value.selection.map { it.id })
    }

    /** `traits.txt` and `extra.txt`, picked together, so the dialog path sniffs both. */
    private fun twoPhenotypes(scope: TestScope, engine: AnalysisEngine = testEngine()): AppViewModel {
        val files = FakeGenotypeFiles(chosenFiles = listOf(traits, extra))
        val viewModel = AppViewModel(files, engine, RecordingResultFiles(), scope)
        viewModel.handle(AppEvent.OpenFile)
        scope.testScheduler.advanceUntilIdle()
        return viewModel
    }

    private fun opened(
        scope: TestScope,
        file: TextGenotypeFile,
        saved: RecordingResultFiles = RecordingResultFiles(),
    ): AppViewModel {
        val viewModel = AppViewModel(FakeGenotypeFiles(chosen = file), testEngine(), saved, scope)
        viewModel.handle(AppEvent.OpenFile)
        scope.testScheduler.advanceUntilIdle()
        return viewModel
    }
}
