package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.ops.FileKind
import net.maizegenetics.tassel.ops.operation

/** Which genotype data `glm` tests (TASSEL 5's `genotypeComponent`). `AlleleProbability` is broken in TASSEL 5 and not ported. */
internal enum class GlmGenotypeComponent { Genotype, ReferenceProbability }

/** Which reports the CLI writes for an association result (tassel-6#200). */
@Suppress("ktlint:standard:enum-entry-name-case")
internal enum class AssociationReport { legacy, modern, both }

/** What to do with a missing genotype. Only `drop` is implemented: the marker's test omits the observation. */
@Suppress("ktlint:standard:enum-entry-name-case")
internal enum class MissingGenotypes { drop, imputeMean }

/**
 * `glm`: TASSEL 5's fixed-effect linear model for every trait and marker of a data set.
 *
 * Declared here with its kernel ([glm]) and registered by tassel-6#200. The parameters that need a later kernel are
 * declared with their TASSEL 5 names, so the command line and the legacy translator keep working, and refused with the
 * task that adds them: permutations (#194), BLUEs (#193), streamed output files (#209), mean imputation.
 *
 * `siteStatsOut` and `siteStatFile` are dead parameters in TASSEL 5 (never read) and are not declared.
 *
 * @see net.maizegenetics.analysis.association.FixedEffectLMPlugin
 */
internal val glmOperation = operation<DataSet, AssociationTable>(
    id = "glm",
    legacyPlugin = "FixedEffectLMPlugin",
    help = "General linear model association of every trait with every marker",
) {
    val traits by stringList(
        help = "The traits to analyze; by default every `data` attribute of the phenotype",
    )
    val covariates by stringList(
        help = "Numeric fixed effects; by default every `covariate` attribute of the phenotype",
    )
    val factors by stringList(
        help = "Categorical fixed effects; by default every `factor` attribute of the phenotype",
    )
    val encoding by enum(
        MarkerEncoding.genotype,
        help = "How a marker enters the model: genotype (each class a level, TASSEL 5's test), additive (copies of " +
            "the minor allele, one df), additiveDominance, or dosage (a numeric value, see genotypeComponent)",
    )
    val genotypeComponent by enum(
        GlmGenotypeComponent.Genotype,
        legacy = "genotypeComponent",
        help = "If the genotype table contains more than one type of genotype data, choose the type to use for the " +
            "analysis. ReferenceProbability tests the reference allele probability as the dosage encoding.",
    )
    val minClassSize by int(
        default = 0,
        range = 0..Int.MAX_VALUE,
        legacy = "minClassSize",
        help = "The minimum acceptable genotype class size. Genotypes in a class with a smaller size will be set " +
            "to missing.",
    )
    val biallelicOnly by boolean(
        default = false,
        legacy = "biallelicOnly",
        help = "Only test sites that are bi-allelic. The alternative is to test sites with two or more alleles.",
    )
    val appendAddDom by boolean(
        default = false,
        legacy = "appendAddDom",
        help = "If true, additive and dominance effect estimates will be added to the stats report for bi-allelic " +
            "sites only. The additive effect will always be non-negative.",
    )
    val maxP by double(
        default = 1.0,
        range = 0.0..1.0,
        legacy = "maxP",
        help = "Only results with p <= maxP will be reported. Default = 1.0.",
    )
    val minMaf by double(
        default = 0.0,
        range = 0.0..0.5,
        help = "Skip markers whose minor allele frequency is below this; 0 is off",
    )
    val minMac by int(
        default = 0,
        range = 0..Int.MAX_VALUE,
        help = "Skip markers with fewer minor allele copies than this; 0 is off",
    )
    val missingGenotypes by enum(
        MissingGenotypes.drop,
        help = "drop omits an observation from the tests of a marker it has no genotype at, as TASSEL 5 does",
    )
    val report by enum(
        AssociationReport.legacy,
        help = "Which reports a command line writes: TASSEL 5's, the modern tab-separated report, or both",
    )
    val permute by boolean(
        default = false,
        legacy = "permute",
        help = "Should a permutation analysis be run? Not yet ported (tassel-6#194).",
    )
    val numberOfPermutations by int(
        default = 0,
        range = 0..Int.MAX_VALUE,
        legacy = "nperm",
        help = "The number of permutations to be run for the permutation analysis. Not yet ported (tassel-6#194).",
    )
    val phenotypeOnly by boolean(
        default = false,
        legacy = "phenoOnly",
        help = "Analyze the phenotype with no markers and generate BLUEs. Not yet ported (tassel-6#193).",
    )
    val saveToFile by boolean(
        default = false,
        legacy = "saveToFile",
        help = "Write the results to files as each site block is analyzed. Not yet ported (tassel-6#209).",
    )
    val siteFile by file(
        kind = FileKind.Output,
        legacy = "siteFile",
        help = "The name of the file to which the statistics are saved. Not yet ported (tassel-6#209).",
    )
    val alleleFile by file(
        kind = FileKind.Output,
        legacy = "alleleFile",
        help = "The name of the file to which the genotype effects are saved. Not yet ported (tassel-6#209).",
    )
    val bluesFile by file(
        kind = FileKind.Output,
        legacy = "bluesFile",
        help = "Name of file to which BLUEs values will be saved. Not yet ported (tassel-6#193).",
    )
    val anovaFile by file(
        kind = FileKind.Output,
        legacy = "anovaFile",
        help = "Name of file to which the ANOVA report will be saved. Not yet ported (tassel-6#193).",
    )

    run { input, progress ->
        notYet(permute() || numberOfPermutations() > 0, "permutations", "#194")
        notYet(phenotypeOnly() || bluesFile() != null || anovaFile() != null, "phenoOnly and BLUEs", "#193")
        notYet(saveToFile() || siteFile() != null || alleleFile() != null, "streamed output files", "#209")
        notYet(missingGenotypes() != MissingGenotypes.drop, "missingGenotypes = ${missingGenotypes()}", "a later task")

        // genotypeComponent ReferenceProbability is TASSEL 5's name for the dosage encoding.
        val markerEncoding = when {
            genotypeComponent() == GlmGenotypeComponent.ReferenceProbability -> {
                require(encoding() == MarkerEncoding.genotype || encoding() == MarkerEncoding.dosage) {
                    "glm: genotypeComponent = ReferenceProbability tests a dosage, not encoding = ${encoding()}"
                }
                MarkerEncoding.dosage
            }
            else -> {
                require(encoding() != MarkerEncoding.dosage) {
                    "glm: encoding = dosage needs genotypeComponent = ReferenceProbability"
                }
                encoding()
            }
        }
        val data = AssociationData.of(input, AssociationSpec(traits(), covariates(), factors()), progress)
        val source = if (markerEncoding == MarkerEncoding.dosage) {
            GenotypeMarkerSource.withDosage(data.data, progress)
        } else {
            GenotypeMarkerSource.of(data.data)
        }
        val options = GlmOptions(
            marker = MarkerOptions(
                encoding = markerEncoding,
                classOrder = ClassOrder.Glm,
                minClassSize = minClassSize(),
                biallelicOnly = biallelicOnly(),
                minMaf = minMaf(),
                minMac = minMac(),
            ),
            appendAddDom = appendAddDom(),
            maxP = maxP(),
        )
        val parameters = values.toMap().mapNotNull { (name, value) -> value?.let { name to it.toString() } }.toMap()
        glmTable(data, source, options, parameters, progress)
    }
}

private fun notYet(requested: Boolean, what: String, task: String) {
    require(!requested) { "glm: $what is not ported yet ($task)" }
}
