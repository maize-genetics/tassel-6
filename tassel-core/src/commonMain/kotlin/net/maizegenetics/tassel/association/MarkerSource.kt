package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.position.Position

/**
 * How a marker becomes model columns (ADR-0026).
 *
 * The omnibus [genotype] test and the 1-df [additive] test are different tests, and a result is labeled by the encoding
 * that produced it.
 */
@Suppress("ktlint:standard:enum-entry-name-case")
internal enum class MarkerEncoding {
    /** Each observed genotype class is a factor level (TASSEL 5's omnibus test, `k - 1` df). */
    genotype,

    /** Copies of one named effect allele (0, 1, 2): one column, or one per non-reference allele when expanded. */
    additive,

    /** [additive] plus a heterozygote indicator, for biallelic sites (TASSEL 5's `add_` and `dom_` tests). */
    additiveDominance,

    /** A numeric value per observation from a numeric source (reference probability today). */
    dosage,
}

/** The order TASSEL 5 puts genotype classes in; the last class is the zero reference in both. */
internal enum class ClassOrder {
    /** GLM: IUPAC strings, `A C G T` first, then the rest by character (`LM/ModelEffectUtils`). */
    Glm,

    /** MLM: homozygotes first by signed byte, then the others by signed byte (`ModelEffectUtils` l.159-179). */
    Mlm,
}

/** Which allele an [MarkerEncoding.additive] column counts. TASSEL 6 has no reference allele (ADR-0026). */
internal sealed interface EffectAllele {
    /** The second most frequent allele among the analysis observations: the default. */
    data object Minor : EffectAllele

    /** The most frequent allele among the analysis observations. */
    data object Major : EffectAllele

    /** The allele with this code (see `Allele`). A marker without it among the observations is skipped. */
    data class Named(val allele: Byte) : EffectAllele
}

/** What [MarkerEncoding.additive] does at a site with more than two alleles; it never collapses them silently. */
internal enum class MultiallelicAdditive {
    /** Skip the marker with [SkipReason.Multiallelic]. */
    Skip,

    /** One column per allele except the most frequent, which is the reference. */
    Expand,
}

/**
 * Encoding and filter settings. Filters see the trait's analysis observations with a call, so a marker can be tested
 * for one trait and skipped for another.
 *
 * @property minClassSize TASSEL 5's `minClassSize`: `0` is off; otherwise a marker with fewer than two classes of at
 * least this size is skipped, and smaller classes become missing for that marker.
 * @property biallelicOnly TASSEL 5's `biallelicOnly`: keep two classes, or three with exactly one heterozygote.
 * @property minMaf Skip markers whose minor allele frequency is below this; `0.0` is off. New in TASSEL 6.
 * @property minMac Skip markers with fewer minor allele copies than this; `0` is off. New in TASSEL 6.
 */
internal data class MarkerOptions(
    val encoding: MarkerEncoding = MarkerEncoding.genotype,
    val classOrder: ClassOrder = ClassOrder.Glm,
    val effectAllele: EffectAllele = EffectAllele.Minor,
    val multiallelic: MultiallelicAdditive = MultiallelicAdditive.Skip,
    val minClassSize: Int = 0,
    val biallelicOnly: Boolean = false,
    val minMaf: Double = 0.0,
    val minMac: Int = 0,
) {
    init {
        require(minClassSize >= 0) { "minClassSize must not be negative: $minClassSize" }
        require(minMac >= 0) { "minMac must not be negative: $minMac" }
        require(minMaf in 0.0..0.5) { "minMaf must be between 0 and 0.5: $minMaf" }
    }
}

/** Why a marker is not tested for an observation set. Skipped markers are reported, never dropped silently. */
internal enum class SkipReason {
    /** No analysis observation has a usable value at the marker. */
    NoObservations,

    /** Fewer than two genotype classes (or alleles) remain, so there is nothing to test. */
    Monomorphic,

    /** `biallelicOnly` and the marker has more classes than a biallelic site can have. */
    NotBiallelic,

    /** `minClassSize` and fewer than two classes are big enough. */
    ClassSize,

    /** The minor allele frequency is below `minMaf`. */
    MinorAlleleFrequency,

    /** There are fewer minor allele copies than `minMac`. */
    MinorAlleleCount,

    /** The encoding needs two alleles and the marker has more. */
    Multiallelic,

    /** The named effect allele is absent from the analysis observations. */
    EffectAlleleAbsent,
}

/**
 * The model columns of one marker for one observation set, written into buffers that [MarkerBlock.load] reuses: read
 * what you need before the next call. Values past [size] are stale.
 *
 * Observations are positions in the set given to `scan`. [observation] lists the kept ones (those with a usable value
 * after filters), in set order, and every per-observation array here is parallel to it. A marker skipped as
 * [SkipReason.Monomorphic] still lists the observations that have a call (and, for the `genotype` encoding, its one
 * class), because a GLM with no marker columns is a model TASSEL 5 reports; every other skip leaves [size] at 0.
 */
internal class MarkerColumns(
    /** The largest observation set this can hold. */
    val capacity: Int,
) {
    /** Why the marker was skipped, or null when it is encoded. */
    var skip: SkipReason? = null
        private set

    /** Number of kept observations. */
    var size: Int = 0
        private set

    /** For each kept observation, its position in the observation set. */
    val observation = IntArray(capacity)

    /**
     * For [MarkerEncoding.genotype]: each kept observation's class, `0 until numberOfClasses`, in TASSEL 5's order.
     * The last class is the zero reference.
     */
    val classIndex = IntArray(capacity)

    /** Genotype classes present, as TASSEL 5's IUPAC strings (`A`, `M`, `0`), in model order. */
    var classLabels: List<String> = emptyList()
        private set

    /** The same classes as explicit allele strings, such as `A:C`. */
    var classAlleleLabels: List<String> = emptyList()
        private set

    /** Kept observations in each class, in model order. */
    var classCounts: IntArray = IntArray(0)
        private set

    /** Number of genotype classes present. */
    val numberOfClasses: Int get() = classLabels.size

    /**
     * TASSEL 5's `minorObs`: the second largest class count among the observations with a call, before
     * `minClassSize` sets small classes to missing; 0 with fewer than two classes.
     */
    var minorObs: Int = 0
        private set

    /** Number of model columns in [column]. */
    var columnCount: Int = 0
        private set

    /** One label per column: the effect allele string, or `dom` for the heterozygote indicator. */
    var columnLabels: List<String> = emptyList()
        private set

    /** The allele each column counts, or `Allele.UNKNOWN` for a column that is not an allele count. */
    var columnAlleles: ByteArray = ByteArray(0)
        private set

    /** The allele an additive column is measured against, or `Allele.UNKNOWN` when it is not a single allele. */
    var comparisonAllele: Byte = UNKNOWN_ALLELE
        private set

    private val buffers = ArrayList<DoubleArray>()

    /**
     * Column [index], one value per kept observation. For [MarkerEncoding.genotype] these are TASSEL 5's `add_F`
     * inputs: one allele-copy count per allele except the last, alleles in frequency order over the genotyped taxa.
     */
    fun column(index: Int): DoubleArray {
        require(index in 0 until columnCount) { "Column $index of $columnCount" }
        return buffers[index]
    }

    internal fun reset() {
        skip = null
        size = 0
        classLabels = emptyList()
        classAlleleLabels = emptyList()
        classCounts = EMPTY_COUNTS
        minorObs = 0
        columnCount = 0
        columnLabels = emptyList()
        columnAlleles = EMPTY_ALLELES
        comparisonAllele = UNKNOWN_ALLELE
    }

    internal fun skipped(reason: SkipReason, minorObs: Int = this.minorObs): MarkerColumns {
        reset()
        skip = reason
        this.minorObs = minorObs
        return this
    }

    /** Marks a marker skipped while keeping the observations and classes already set (a monomorphic marker). */
    internal fun flagSkipped(reason: SkipReason) {
        skip = reason
    }

    internal fun kept(size: Int, minorObs: Int) {
        this.size = size
        this.minorObs = minorObs
    }

    internal fun classes(labels: List<String>, alleleLabels: List<String>, counts: IntArray) {
        classLabels = labels
        classAlleleLabels = alleleLabels
        classCounts = counts
    }

    internal fun columns(labels: List<String>, alleles: ByteArray, comparison: Byte) {
        columnCount = labels.size
        columnLabels = labels
        columnAlleles = alleles
        comparisonAllele = comparison
        while (buffers.size < columnCount) buffers.add(DoubleArray(capacity))
    }

    private companion object {
        const val UNKNOWN_ALLELE: Byte = 0xF
        val EMPTY_COUNTS = IntArray(0)
        val EMPTY_ALLELES = ByteArray(0)
    }
}

/** A run of consecutive markers from [MarkerSource.scan], read from the underlying data once. */
internal interface MarkerBlock {
    /** Marker indices in this block. */
    val markers: IntRange

    /**
     * Encodes [marker] (one of [markers]) for observation set [set] (an index into the list given to `scan`) into
     * [into], which is reset first. Returns true when the marker is encoded, false when [MarkerColumns.skip] says why
     * not. Calls for different sets at the same marker reuse the one read of its genotypes.
     */
    fun load(marker: Int, set: Int, into: MarkerColumns): Boolean
}

/**
 * Turns the markers of some data set into model columns for chosen observation sets: the seam between association
 * models and what a marker is (ADR-0026). [GenotypeMarkerSource] is the only implementation in this epic; dosages,
 * probabilities, haplotypes, and structural variants plug in behind the same interface.
 *
 * Observations are positions in the analysis data (the rows of a `GenotypePhenotype`); an observation set is the
 * subset a trait uses, normally the rows where the trait and every covariate are present.
 */
internal interface MarkerSource {
    /** Number of markers. */
    val numberOfMarkers: Int

    /** The site of [marker], for names, chromosome, and position. */
    fun position(marker: Int): Position

    /**
     * Reads the markers in blocks, once, and lets each block be encoded for every set in [observationSets] (one per
     * trait that shares the data). Each set lists observations without repeats, in the order columns should follow.
     */
    fun scan(
        observationSets: List<IntArray>,
        options: MarkerOptions = MarkerOptions(),
        blockSize: Int = DEFAULT_MARKER_BLOCK_SIZE,
    ): Sequence<MarkerBlock>

    companion object {
        /** Default markers per block, as `GenotypeTable.DEFAULT_SITE_BLOCK_SIZE`. */
        const val DEFAULT_MARKER_BLOCK_SIZE: Int = 4096
    }
}
