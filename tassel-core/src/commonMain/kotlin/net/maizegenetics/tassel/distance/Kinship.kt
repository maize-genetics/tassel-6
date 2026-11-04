package net.maizegenetics.tassel.distance

import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.ops.operation

/**
 * Kinship methods TASSEL 5 offers. The constant names are TASSEL 5's, because they are what
 * `-method` accepts on the command line and what the `Matrix_Type` annotation holds in the written
 * matrix.
 *
 * @see net.maizegenetics.analysis.distance.KinshipPlugin.KINSHIP_METHOD
 */
@Suppress("ktlint:standard:enum-entry-name-case")
enum class KinshipMethod {
    /**
     * Endelman's centered IBS, scaled to give a reasonable estimate of additive genetic variance.
     * [Equation 13](http://www.g3journal.org/content/2/11/1405.full.pdf).
     */
    Centered_IBS,
}

/**
 * Kinship from genotypes. TASSEL 5's `Normalized_IBS` and the two dominance methods are not ported
 * yet, so [KinshipMethod] has one constant and the parameter exists to keep the command line and
 * the legacy translator honest.
 *
 * @see net.maizegenetics.analysis.distance.KinshipPlugin
 */
val kinshipOperation = operation<GenotypeTable, DistanceMatrix>(
    id = "kinship",
    legacyPlugin = "KinshipPlugin",
    help = "Kinship matrix for every pair of taxa",
) {
    val method by enum(
        default = KinshipMethod.Centered_IBS,
        legacy = "method",
        help = "The Centered_IBS (Endelman - previously Scaled_IBS) method produces a kinship matrix that is " +
            "scaled to give a reasonable estimate of additive genetic variance. Uses algorithm " +
            "http://www.g3journal.org/content/2/11/1405.full.pdf Equation-13.",
    )
    val maxAlleles by int(
        default = DEFAULT_MAX_ALLELES,
        range = 2..DEFAULT_MAX_ALLELES,
        legacy = "maxAlleles",
        help = "Allele states per site to evaluate, most frequent first",
    )
    run { input, progress ->
        when (method()) {
            KinshipMethod.Centered_IBS -> centeredIbs(input, maxAlleles(), progress)
        }
    }
}
