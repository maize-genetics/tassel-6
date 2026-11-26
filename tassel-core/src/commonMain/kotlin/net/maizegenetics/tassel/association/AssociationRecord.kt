package net.maizegenetics.tassel.association

/** Which model produced a result. */
internal enum class AssociationModel { Glm, Mlm }

/** Which test a record holds (ADR-0027). The omnibus genotype test comes first among the records of a marker. */
internal enum class AssociationTestKind {
    /** The omnibus test of the genotype classes (TASSEL 5's `marker_F`, `k - 1` df). */
    Genotype,

    /** The test of the allele-copy (dosage) columns (TASSEL 5's `add_F`). */
    Additive,

    /** The test of the heterozygote deviation (TASSEL 5's `dom_F`). */
    Dominance,

    /** A single numeric value per observation, such as a reference probability. */
    Dosage,
}

/** Whether a (trait, marker, test) was tested. */
internal enum class AssociationStatus { Ok, Skipped, Failed }

/** How the variance components behind an MLM record were obtained. [None] for GLM. */
internal enum class VarianceEstimation {
    None,

    /** One estimate from the null model, reused for every marker (P3D). */
    P3d,

    /** Re-estimated at each marker. */
    EachMarker,
}

/**
 * One result per (trait, marker, test) with the fixed schema of ADR-0027. A number that does not apply is `NaN`, a
 * count that does not apply is `-1`, and text that does not apply is empty.
 *
 * The records of one (trait, marker) are produced adjacent, omnibus test first; sinks that filter by p-value rely on it.
 *
 * `modelDf` and the three mean squares carry TASSEL 5's `model_df`, `marker_MS`, `error_MS`, and `model_MS` so the
 * compatible GLM report can be written from the table. They are `NaN` outside GLM.
 */
internal data class AssociationRecord(
    // Identity.
    val trait: String,
    val markerIndex: Int,
    val marker: String,
    val chromosome: String,
    val position: Int,
    // Test.
    val model: AssociationModel = AssociationModel.Glm,
    val test: AssociationTestKind = AssociationTestKind.Genotype,
    val encoding: MarkerEncoding = MarkerEncoding.genotype,
    val effectAllele: String = "",
    val comparisonAllele: String = "",
    // Counts.
    val observations: Int = -1,
    val taxa: Int = -1,
    val alleleCount: Int = -1,
    val alleleFrequency: Double = Double.NaN,
    val minorAlleleCount: Int = -1,
    /** TASSEL 5's `minorObs`: the second largest genotype class among the observations. */
    val minorObs: Int = -1,
    // Estimates.
    val effect: Double = Double.NaN,
    val standardError: Double = Double.NaN,
    val ciLower: Double = Double.NaN,
    val ciUpper: Double = Double.NaN,
    // Statistics.
    val statistic: Double = Double.NaN,
    val numeratorDf: Double = Double.NaN,
    val denominatorDf: Double = Double.NaN,
    val p: Double = Double.NaN,
    val log10P: Double = Double.NaN,
    val rSquared: Double = Double.NaN,
    val markerMeanSquare: Double = Double.NaN,
    val errorMeanSquare: Double = Double.NaN,
    val modelDf: Double = Double.NaN,
    val modelMeanSquare: Double = Double.NaN,
    /** TASSEL 5's `perm_p`. */
    val permutationP: Double = Double.NaN,
    // Variance components (MLM).
    val geneticVariance: Double = Double.NaN,
    val residualVariance: Double = Double.NaN,
    val heritability: Double = Double.NaN,
    val varianceEstimation: VarianceEstimation = VarianceEstimation.None,
    // Status.
    val status: AssociationStatus = AssociationStatus.Ok,
    val skipReason: SkipReason? = null,
    val message: String = "",
) {
    init {
        require((status == AssociationStatus.Skipped) == (skipReason != null)) {
            "skipReason must be set exactly when the status is Skipped: $status, $skipReason"
        }
    }
}

/**
 * The effect of one genotype class or allele (TASSEL 5's genotype-effects report), keyed by trait and marker.
 *
 * @property allele The genotype class or allele, or empty when the effect is not per class (reference probability).
 * @property observations Observations in the class, or `-1` when not known.
 */
internal data class AssociationEffect(
    val trait: String,
    val markerIndex: Int,
    val marker: String,
    val chromosome: String,
    val position: Int,
    val allele: String = "",
    val observations: Int = -1,
    val estimate: Double,
    val standardError: Double = Double.NaN,
)
