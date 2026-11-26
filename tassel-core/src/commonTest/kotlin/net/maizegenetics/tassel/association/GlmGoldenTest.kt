package net.maizegenetics.tassel.association

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.filter.SiteFilter
import net.maizegenetics.tassel.filter.filterSites
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.phenotype.PhenotypeFiles
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.testing.SuspendLazy
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.assertMatchesGolden
import kotlin.test.Test

/**
 * `glm` against every TASSEL 5 GLM golden (#188, ADR-0009), through the TASSEL 5-compatible writers. The rounded reports
 * are compared in format and the full-precision twins at the MANIFEST tolerance (`rel 1e-8`). Departures from TASSEL 5
 * are named in [GlmDepartures] and excluded cell by cell, never by loosening a tolerance (ADR-0029).
 */
class GlmGoldenTest {
    private class Inputs {
        val tutorial = SuspendLazy { hapmap("tutorial/mdp_genotype.hmp.txt") }
        val tutorialMaf = SuspendLazy { filterSites(tutorial(), SiteFilter(minMaf = 0.05)) }
        val chr9and10 = SuspendLazy { hapmap("tutorial/mdp_genotype_chr9_10.hmp.txt") }
        val synthetic = SuspendLazy { hapmap("association/glm_synthetic.hmp.txt") }
        val traits = SuspendLazy { phenotype("tutorial/mdp_traits.txt") }
        val structure = SuspendLazy { phenotype("tutorial/mdp_population_structure.txt") }
        val syntheticTraits = SuspendLazy { phenotype("association/glm_synthetic_traits.txt") }
        val syntheticFactor = SuspendLazy { phenotype("association/glm_synthetic_factor.txt") }
        val syntheticCollinear = SuspendLazy { phenotype("association/glm_synthetic_collinear.txt") }
    }

    private companion object {
        val inputs = Inputs()

        /** `-excludeLastTrait` on the population structure file drops Q3, since Q1 + Q2 + Q3 = 1. */
        val tutorialSpec = AssociationSpec(covariates = listOf("Q1", "Q2"))

        suspend fun hapmap(path: String): GenotypeTable =
            HapMap.read(Buffer().also { it.write(TestData.readBytes(path)) })

        suspend fun phenotype(path: String): PhenotypeTable =
            PhenotypeFiles.read(Buffer().also { it.write(TestData.readBytes(path)) }).table
    }

    private suspend fun check(
        directory: String,
        name: String,
        input: DataSet,
        spec: AssociationSpec = AssociationSpec(),
        options: GlmOptions = GlmOptions(),
        dosage: Boolean = false,
    ) {
        val data = AssociationData.of(input, spec)
        val source = if (dosage) GenotypeMarkerSource.withDosage(data.data) else GenotypeMarkerSource.of(data.data)
        val table = glmTable(data, source, options)
        val layout = if (dosage) GlmLayout.ReferenceProbability else GlmLayout.Genotype
        fun stats(form: NumberForm) = buildString {
            table.writeGlmStatistics(this, layout, form, appendAddDom = options.appendAddDom)
        }

        fun effects(form: NumberForm) = buildString { table.writeGlmEffects(this, layout, form) }

        val label = "$directory/$name"
        assertMatchesGolden("golden/$label.stats.txt", stats(NumberForm.Tassel))
        assertMatchesGolden("golden/$label.effects.txt", effects(NumberForm.Tassel))
        GlmDepartures.assertFullMatches(label, "stats", stats(NumberForm.Full))
        GlmDepartures.assertFullMatches(label, "effects", effects(NumberForm.Full))
    }

    @Test
    fun tutorialWithPopulationStructure() = runTest {
        check(
            "glm-tutorial",
            "mdp_traits",
            DataSet.of(inputs.tutorial(), inputs.traits(), inputs.structure()),
            tutorialSpec,
        )
    }

    @Test
    fun tutorialWithMafFilter() = runTest {
        check(
            "glm-tutorial-maf",
            "mdp_traits",
            DataSet.of(inputs.tutorialMaf(), inputs.traits(), inputs.structure()),
            tutorialSpec,
        )
    }

    @Test
    fun noCovariates() = runTest {
        check("glm-no-covariates", "mdp_traits", DataSet.of(inputs.chr9and10(), inputs.traits()))
    }

    @Test
    fun referenceProbability() = runTest {
        check(
            "glm-refprob",
            "mdp_traits",
            DataSet.of(inputs.chr9and10(), inputs.traits()),
            options = GlmOptions(marker = MarkerOptions(encoding = MarkerEncoding.dosage)),
            dosage = true,
        )
    }

    private suspend fun options(name: String, options: GlmOptions) = check(
        "glm-options",
        name,
        DataSet.of(inputs.chr9and10(), inputs.traits(), inputs.structure()),
        tutorialSpec,
        options,
    )

    @Test
    fun minClassSize() = runTest { options("min-class-size-5", GlmOptions(marker = MarkerOptions(minClassSize = 5))) }

    @Test
    fun biallelicOnly() = runTest {
        options("biallelic-only", GlmOptions(marker = MarkerOptions(biallelicOnly = true)))
    }

    @Test
    fun maxP() = runTest { options("max-p-0.01", GlmOptions(maxP = 0.01)) }

    @Test
    fun appendAddDom() = runTest { options("append-add-dom", GlmOptions(appendAddDom = true)) }

    private suspend fun synthetic(
        name: String,
        options: GlmOptions = GlmOptions(),
        phenotype: SuspendLazy<PhenotypeTable> = inputs.syntheticTraits,
    ) = check("glm-synthetic", name, DataSet.of(inputs.synthetic(), phenotype()), options = options)

    @Test
    fun syntheticModel() = runTest { synthetic("model") }

    @Test
    fun syntheticMinClassSize() = runTest {
        synthetic("min-class-size-3", GlmOptions(marker = MarkerOptions(minClassSize = 3)))
    }

    @Test
    fun syntheticBiallelicOnly() = runTest {
        synthetic("biallelic-only", GlmOptions(marker = MarkerOptions(biallelicOnly = true)))
    }

    @Test
    fun syntheticAppendAddDom() = runTest { synthetic("append-add-dom", GlmOptions(appendAddDom = true)) }

    @Test
    fun syntheticFactorOnly() = runTest { synthetic("factor-only", phenotype = inputs.syntheticFactor) }

    @Test
    fun syntheticCollinearCovariates() = runTest { synthetic("collinear", phenotype = inputs.syntheticCollinear) }
}
