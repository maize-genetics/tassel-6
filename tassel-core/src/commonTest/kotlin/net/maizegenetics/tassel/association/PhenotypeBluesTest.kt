package net.maizegenetics.tassel.association

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeFiles
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.stats.fTestPValue
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.assertMatchesGolden
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * Phenotype-only GLM and BLUEs (#193): the three TASSEL 5 goldens in `glm-phenotype-only/` through the legacy writers,
 * TASSEL 5's own R-computed BLUEs (`PhenotypeLMTest`) as a second oracle, hand-solved cases, and the BLUEs as the
 * phenotype of a later `glm`.
 */
class PhenotypeBluesTest {
    private companion object {
        suspend fun phenotype(path: String): PhenotypeTable =
            PhenotypeFiles.read(Buffer().also { it.write(TestData.readBytes(path)) }).table

        /** The BLUEs `PhenotypeLMTest` checks, computed in R, in the order of the taxa of `earht_with_rep.txt`. */
        val R_BLUES = floatArrayOf(
            94.0315456025016f, 64.7254997617167f, 45.6519525059068f, 29.7524935013531f, 62.5555279926411f,
            38.8168399291231f, 67.7808311612335f, 70.8414557486813f, 45.2450936651001f, 63.8476743870665f,
            66.9296408321963f, 137.361806914807f, 83.4337287433577f, 78.2158230657591f, 66.1571330200281f,
            58.7176743870664f, 55.9297855278656f, 43.3691463840685f, 68.9777451837748f, 66.3652248691787f,
            77.0025258082017f, 52.3935891603886f, 93.7872183320187f, 72.7176743870665f, 45.7176743870664f,
            51.9557860580037f, 85.6218663548891f, 62.9402081302404f, 67.2176743870665f, 62.8193914236268f,
        )
    }

    private suspend fun read(text: String): PhenotypeTable =
        PhenotypeFiles.read(Buffer().also { it.write(text.encodeToByteArray()) }).table

    private suspend fun checkGolden(name: String, input: PhenotypeTable, spec: AssociationSpec = AssociationSpec()) {
        val result = phenotypeLm(input, spec)
        val label = "glm-phenotype-only/$name"
        fun anova(form: NumberForm) = buildString { result.anova.writePhenotypeAnova(this, form) }
        assertMatchesGolden("golden/$label.blues.txt", result.blues.toDelimitedString())
        assertMatchesGolden("golden/$label.anova.txt", anova(NumberForm.Tassel))
        GlmDepartures.assertFullMatches(label, "anova", anova(NumberForm.Full))
    }

    @Test
    fun replicatedTaxaMatchTassel5() = runTest {
        checkGolden("earht_with_rep", phenotype("phenotype/earht_with_rep.txt"))
    }

    @Test
    fun unreplicatedTaxaMatchTassel5ExceptTheNoiseCells() = runTest {
        checkGolden("mdp_traits", phenotype("tutorial/mdp_traits.txt"))
    }

    @Test
    fun factorCovariateAndMissingValuesMatchTassel5() = runTest {
        checkGolden("synthetic", phenotype("association/glm_phenotype_only.txt"))
    }

    @Test
    fun bluesMatchTheRValuesInTassel5sTest() = runTest {
        val blues = phenotypeLm(phenotype("phenotype/earht_with_rep.txt")).blues
        val observed = blues.attribute("earht") as NumericAttribute
        assertEquals(R_BLUES.size, observed.size)
        for (i in R_BLUES.indices) {
            assertTrue(abs(R_BLUES[i] - observed.floatValue(i)) <= 1e-4f, "BLUE $i: ${observed.floatValue(i)}")
        }
    }

    @Test
    fun unreplicatedTaxaGiveErrorDfZeroAndNaN() = runTest {
        val row = phenotypeLm(read("<Phenotype>\ntaxa\tdata\nTaxa\ty\nA\t1\nB\t4\nC\t9\n")).anova.single()
        assertEquals(0.0, row.errorDf)
        assertEquals(2.0, row.taxaDf)
        assertTrue(row.f.isNaN() && row.p.isNaN() && row.errorMs.isNaN())
        assertEquals(row.taxaMs, row.modelMs)
    }

    @Test
    fun oneWayAnovaIsSolvedByHand() = runTest {
        // Taxon means 2 and 6; SS(taxa) = 2 * 2^2 + 2 * 2^2 = 16 on 1 df; SS(error) = 4 on 2 df; F = 8.
        val result = phenotypeLm(read("<Phenotype>\ntaxa\tdata\nTaxa\ty\nA\t1\nA\t3\nB\t5\nB\t7\n"))
        val row = result.anova.single()
        assertEquals(8.0, row.f, 1e-12)
        assertEquals(fTestPValue(8.0, 1.0, 2.0), row.p, 1e-12)
        assertEquals(listOf(1.0, 2.0, 1.0), listOf(row.taxaDf, row.errorDf, row.modelDf))
        assertEquals(16.0, row.taxaMs, 1e-12)
        assertEquals(2.0, row.errorMs, 1e-12)
        assertEquals(16.0, row.modelMs, 1e-12)
        val blues = result.blues.attribute("y") as NumericAttribute
        assertEquals(listOf(2f, 6f), listOf(blues.floatValue(0), blues.floatValue(1)))
    }

    @Test
    fun bluesAverageOverFactorLevelsAndAtCovariateMeansOnlyWhenAsked() = runTest {
        val table = read(
            "<Phenotype>\ntaxa\tfactor\tcovariate\tdata\nTaxa\trep\tx\ty\n" +
                "A\t1\t1\t1\nA\t2\t3\t11\nB\t1\t2\t5\nB\t2\t5\t15\nC\t1\t4\t9\nC\t2\t6\t19.5\nC\t1\t7\t9.5\n",
        )
        fun blues(spec: AssociationSpec = AssociationSpec(), options: PhenotypeLmOptions = PhenotypeLmOptions()) =
            (phenotypeLm(table, spec, options).blues.attribute("y") as NumericAttribute).toFloatArray()

        val noCovariate = AssociationSpec(covariates = emptyList())
        // With no covariate the option has nothing to add.
        assertTrue(blues(noCovariate).contentEquals(blues(noCovariate, PhenotypeLmOptions(atCovariateMeans = true))))
        // With one, the default evaluates at x = 0 and the option adds the same slope times mean(x) to every taxon.
        val shifts = blues(options = PhenotypeLmOptions(atCovariateMeans = true)).zip(blues()) { a, b -> a - b }
        assertTrue(shifts.all { abs(it - shifts[0]) < 1e-4f } && abs(shifts[0]) > 1e-3f, "$shifts")
    }

    @Test
    fun atCovariateMeansAddsSlopeTimesMeanToEveryBlue() = runTest {
        val table = phenotype("association/glm_phenotype_only.txt")
        val spec = AssociationSpec(traits = listOf("y1"))
        val atZero = phenotypeLm(table, spec).blues.attribute("y1") as NumericAttribute
        val atMeans = phenotypeLm(table, spec, PhenotypeLmOptions(atCovariateMeans = true)).blues
            .attribute("y1") as NumericAttribute
        val shifts = (0 until atZero.size).filter { !atZero.isMissing(it) }.map { atMeans.value(it) - atZero.value(it) }
        // The data were made with a slope of 0.8 on a covariate with mean near 5.
        assertTrue(shifts.all { abs(it - shifts[0]) < 1e-4 }, "$shifts")
        assertTrue(shifts[0] in 2.0..6.0, "${shifts[0]}")
    }

    @Test
    fun aTaxonWithNoValueOfATraitHasAMissingBlue() = runTest {
        val result = phenotypeLm(phenotype("association/glm_phenotype_only.txt"))
        val y2 = result.blues.attribute("y2") as NumericAttribute
        val taxa = result.blues.taxaAttribute
        val t05 = (0 until taxa.size).first { taxa.taxon(it) == "T05" }
        assertTrue(y2.isMissing(t05))
        assertEquals(12, taxa.size)
        assertEquals(1, y2.missingCount)
    }

    @Test
    fun bluesFeedAGlm() = runTest {
        val genotypes = HapMap.read(
            Buffer().also {
                it.write(TestData.readBytes("tutorial/mdp_genotype_chr9_10.hmp.txt"))
            },
        )
        val blues = phenotypeLm(phenotype("phenotype/earht_with_rep.txt")).blues
        // The BLUEs survive the file format, so they can be saved and loaded as any phenotype.
        val reloaded = read(blues.toDelimitedString())
        assertEquals(blues.attributes.map { it.name }, reloaded.attributes.map { it.name })

        suspend fun table(phenotype: PhenotypeTable): AssociationTable {
            val data = AssociationData.of(DataSet.of(genotypes, phenotype))
            return glmTable(data, GenotypeMarkerSource.of(data.data))
        }

        val direct = table(blues)
        val viaFile = table(reloaded)
        assertEquals(listOf("earht"), direct.traits)
        assertTrue(direct.size > 0)
        // Taxa without genotypes are dropped, and the rest are unreplicated: a plain one-way test per marker.
        assertEquals(direct.size, viaFile.size)
        for (i in 0 until direct.size) {
            assertEquals(direct[i].statistic.isNaN(), viaFile[i].statistic.isNaN())
            if (!direct[i].statistic.isNaN()) {
                assertEquals(
                    direct[i].statistic,
                    viaFile[i].statistic,
                    1e-9 * abs(direct[i].statistic),
                )
            }
        }
    }
}
