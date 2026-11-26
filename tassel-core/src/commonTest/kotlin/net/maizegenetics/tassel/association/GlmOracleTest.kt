package net.maizegenetics.tassel.association

import kotlinx.coroutines.test.runTest
import kotlinx.io.Buffer
import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.io.HapMap
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeFiles
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.stats.model.Anova
import net.maizegenetics.tassel.stats.model.Contrast
import net.maizegenetics.tassel.stats.model.DesignMatrix
import net.maizegenetics.tassel.stats.model.ModelSpec
import net.maizegenetics.tassel.stats.model.ModelTerm
import net.maizegenetics.tassel.taxa.TaxaList
import net.maizegenetics.tassel.testing.OracleManifest
import net.maizegenetics.tassel.testing.TestData
import net.maizegenetics.tassel.testing.assertDelimitedMatches
import net.maizegenetics.tassel.testing.goldenText
import kotlin.math.abs
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * The GLM kernel against R (`testdata/oracle/models/`, ADR-0029 kind 2), which TASSEL 5 cannot be for what it does not
 * do or does wrongly: replicated taxa (decision Q1), allele-copy coding at multiallelic sites, `addEffect` and
 * `domEffect` from the marker's own coefficients, and a near-collinear covariate (departures D1 to D3 in
 * the `[glm-departures]` MANIFEST section). Statistics are compared at the manifest's `rel 1e-10`, degrees of freedom exactly.
 */
class GlmOracleTest {
    private fun fixture(name: String) = TestData.readOracleText("models/$name")

    private fun tolerance(name: String) = OracleManifest.load().tolerance("oracle/models/$name")

    private fun check(name: String, actual: String) =
        assertDelimitedMatches(fixture(name), actual, tolerance(name), source = "oracle/models/$name")

    private suspend fun hapmap(path: String): GenotypeTable =
        HapMap.read(Buffer().also { it.write(TestData.readBytes(path)) })

    private suspend fun phenotype(path: String): PhenotypeTable =
        PhenotypeFiles.read(Buffer().also { it.write(TestData.readBytes(path)) }).table

    // --- replicated taxa ---------------------------------------------------------------------------------------

    private fun tsvRows(text: String) = text.trim().lines().map { it.split('\t') }

    private fun replicateInput(): DataSet {
        val genotypeRows = tsvRows(fixture("glm-replicates.genotypes.tsv"))
        val sites = genotypeRows[0].drop(1)
        val taxa = genotypeRows.drop(1).map { it[0] }
        val calls = genotypeRows.drop(1).map { row -> row.drop(1).map { Nucleotide.genotype(it) } }
        val genotypes = InMemoryGenotypeTable.of(
            TaxaList.of(taxa),
            PositionList(sites.mapIndexed { i, name -> Position(Chromosome.of(1), (i + 1) * 100, name) }),
        ) { taxon, site -> calls[taxon][site] }
        val phenotypeRows = tsvRows(fixture("glm-replicates.phenotype.tsv")).drop(1)
        val phenotype = PhenotypeTable(
            listOf(
                PhenotypeTable.taxa(phenotypeRows.map { it[0] }),
                NumericAttribute(
                    "y",
                    AttributeType.data,
                    FloatArray(phenotypeRows.size) {
                        phenotypeRows[it][1].toFloat()
                    },
                ),
                NumericAttribute(
                    "c1",
                    AttributeType.covariate,
                    FloatArray(phenotypeRows.size) { phenotypeRows[it][2].toFloat() },
                ),
            ),
        )
        return DataSet.of(genotypes, phenotype)
    }

    @Test
    fun replicatedTaxaTestTheMarkerAgainstTaxaWithinMarkerLikeR() = runTest {
        val data = AssociationData.of(replicateInput())
        assertTrue(data.hasReplicates(0))
        val table = glmTable(data, GenotypeMarkerSource.of(data.data))
        check("glm-replicates.tests.tsv", testRows(table))
        // The error term is named in every record, and the unique taxa are counted apart from the observations.
        val main = table.records.filter { it.test == AssociationTestKind.Genotype }.toList()
        assertTrue(main.all { it.message == TAXA_WITHIN_MARKER })
        assertTrue(main.all { it.taxa in 1 until it.observations })
        // r3 has three taxa with no call: all of their replicates are dropped.
        val r3 = main.single { it.marker == "r3" }
        val r1 = main.single { it.marker == "r1" }
        assertEquals(53, r3.observations)
        assertEquals(62, r1.observations)
        assertEquals(24, r1.taxa)
        assertEquals(21, r3.taxa)
    }

    // --- the glm-synthetic departures --------------------------------------------------------------------------

    private suspend fun syntheticTable(
        traits: String = "association/glm_synthetic_traits.txt",
        options: GlmOptions = GlmOptions(),
    ): AssociationTable {
        val data = AssociationData.of(
            DataSet.of(hapmap("association/glm_synthetic.hmp.txt"), phenotype(traits)),
        )
        return glmTable(data, GenotypeMarkerSource.of(data.data), options)
    }

    @Test
    fun allelesAreCountedByCopyAtAMultiallelicSiteLikeR() = runTest {
        check("glm-synthetic.tests.tsv", testRows(syntheticTable()))
    }

    @Test
    fun addAndDomEffectsUseTheMarkersOwnCoefficientsLikeR() = runTest {
        val table = syntheticTable(options = GlmOptions(appendAddDom = true))
        val estimates = table.effects.groupBy { it.trait to it.markerIndex }
        val rows = ArrayList<List<String>>()
        for (row in table.recordOrder()) {
            val record = table[row]
            if (record.test != AssociationTestKind.Additive || record.effect.isNaN()) continue
            val dominance = table[row + 1]
            val classes = estimates.getValue(record.trait to record.markerIndex)
            rows += listOf(
                record.trait,
                record.marker,
                classes.size.toString(),
                number(classes[0].estimate),
                if (classes.size == 3) number(classes[1].estimate) else "NA",
                number(record.effect),
                number(dominance.effect),
            )
        }
        val header = "trait\tmarker\tclasses\tb0\tb1\taddEffect\tdomEffect"
        check(
            "glm-synthetic.effects.tsv",
            (listOf(header) + rows.map { it.joinToString("\t") }).joinToString("\n") + "\n",
        )
    }

    @Test
    fun aNearCollinearCovariateIsAliasedLikeR() = runTest {
        val table = syntheticTable("association/glm_synthetic_collinear.txt")
        check("glm-collinear.tests.tsv", testRows(table))
    }

    /** D1: TASSEL 5's coding (2 for the allele's homozygote, 1 for any heterozygote) reproduces its own add_F. */
    @Test
    fun tassel5sHeterozygoteCodingReproducesItsAddFAtAMultiallelicSite() = runTest {
        val genotypes = hapmap("association/glm_synthetic.hmp.txt")
        val traits = phenotype("association/glm_synthetic_traits.txt")
        val golden = tsvRows(goldenText("golden/glm-synthetic/model.stats.full.txt"))
        val header = golden[0]
        val site = (0 until genotypes.numberOfSites).first { genotypes.positions[it].name == "multiallelic" }
        for (trait in listOf("y1", "y2")) {
            val expected = golden.first {
                it[0] == trait && it[1] == "multiallelic"
            }[header.indexOf("add_F")].toDouble()
            assertEquals(expected, tassel5AddF(genotypes, traits, site, trait), abs(expected) * 1e-8)
        }
    }

    private fun assertEquals(expected: Double, actual: Double, delta: Double) {
        assertTrue(abs(expected - actual) <= delta, "expected $expected within $delta but was $actual")
    }

    private fun tassel5AddF(genotypes: GenotypeTable, traits: PhenotypeTable, site: Int, trait: String): Double {
        val calls = genotypes.siteGenotypes(site)
        val alleles = net.maizegenetics.tassel.genotype.AlleleStats.of(calls).alleles
        val observed = (0 until traits.numberOfObservations).map { traits.taxon(it) }
        fun coded(allele: Byte) = FloatArray(traits.numberOfObservations) { row ->
            val label = Nucleotide.iupac(calls[genotypes.taxa.indexOf(observed[row])]).toString()
            when {
                label == Nucleotide.alleleString(allele) -> 2f
                label in "RSYWKM0".map { it.toString() } -> 1f
                else -> 0f
            }
        }
        val extra = alleles.dropLast(1).mapIndexed { i, allele ->
            NumericAttribute("a$i", AttributeType.covariate, coded(allele))
        }
        val table = PhenotypeTable(traits.attributes + extra)
        val spec = ModelSpec(
            trait,
            listOf(
                ModelTerm.Intercept,
                ModelTerm.Factor("grp", Contrast.treatment),
                ModelTerm.Covariate("c1"),
                ModelTerm.Covariate("c2"),
            ) + extra.map { ModelTerm.Covariate(it.name) },
        )
        val anova = Anova.sequential(DesignMatrix.build(spec, table))
        val added = extra.map { anova.row(it.name)!! }
        val ss = added.sumOf { it.sumOfSquares }
        val df = added.sumOf { it.df }
        return (ss / df) / (anova.residualSumOfSquares / anova.residualDf)
    }

    // --- rendering ---------------------------------------------------------------------------------------------

    private fun number(value: Double): String = if (value.isNaN()) "NA" else value.toString()

    private fun count(value: Double): String = if (value.isNaN()) "NA" else value.toInt().toString()

    /** One row per (trait, marker) tested, in the columns of the R fixtures. */
    private fun testRows(table: AssociationTable): String {
        val classes = table.effects.groupBy { it.trait to it.markerIndex }.mapValues { it.value.size }
        val lines = ArrayList<String>()
        lines += "trait\tmarker\tn\tclasses\tmarker_df\terror_df\tmarker_F\tmarker_p\tmarker_log10p\tmarker_Rsq\t" +
            "marker_MS\terror_MS\tmodel_df\tmodel_MS\tadd_df\tadd_F\tadd_p\tdom_df\tdom_F\tdom_p"
        val order = table.recordOrder()
        var i = 0
        while (i < order.size) {
            val main = table[order[i]]
            var additive: AssociationRecord? = null
            var dominance: AssociationRecord? = null
            var j = i + 1
            while (j < order.size &&
                table[order[j]].let { it.trait == main.trait && it.markerIndex == main.markerIndex }
            ) {
                val next = table[order[j]]
                if (next.test == AssociationTestKind.Additive) additive = next
                if (next.test == AssociationTestKind.Dominance) dominance = next
                j++
            }
            i = j
            if (main.status != AssociationStatus.Ok) continue
            lines += listOf(
                main.trait,
                main.marker,
                main.observations.toString(),
                classes.getValue(main.trait to main.markerIndex).toString(),
                count(main.numeratorDf),
                count(main.denominatorDf),
                number(main.statistic),
                number(main.p),
                number(main.log10P),
                number(main.rSquared),
                number(main.markerMeanSquare),
                number(main.errorMeanSquare),
                count(main.modelDf),
                number(main.modelMeanSquare),
                count(additive?.numeratorDf ?: Double.NaN),
                number(additive?.statistic ?: Double.NaN),
                number(additive?.p ?: Double.NaN),
                count(dominance?.numeratorDf ?: Double.NaN),
                number(dominance?.statistic ?: Double.NaN),
                number(dominance?.p ?: Double.NaN),
            ).joinToString("\t")
        }
        return lines.joinToString("\n") + "\n"
    }
}
