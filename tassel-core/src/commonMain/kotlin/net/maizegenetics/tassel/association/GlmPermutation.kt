package net.maizegenetics.tassel.association

import kotlin.math.floor

/**
 * Settings of TASSEL 5's GLM permutation test: [permutations] shuffles of each trait's base-model residuals, drawn from a
 * generator seeded with [seed]. The seed is always known, so a run can be repeated: when the user gives none, the caller
 * draws one and it is recorded in the result's provenance (ADR-0027).
 */
internal data class PermutationOptions(val permutations: Int, val seed: Long) {
    init {
        require(permutations > 0) { "permutations must be positive: $permutations" }
    }
}

/**
 * The generator behind the permutations: SplitMix64, whose 64-bit state and `Long` arithmetic give the same sequence on
 * every target and in every Kotlin version, which `kotlin.random.Random` does not promise. It replaces TASSEL 5's unseeded
 * `java.util.Random` (departure D6, ADR-0029).
 */
internal class PermutationRandom(seed: Long) {
    private var state = seed

    fun nextLong(): Long {
        state += GOLDEN_GAMMA
        var z = state
        z = (z xor (z ushr 30)) * MIX_1
        z = (z xor (z ushr 27)) * MIX_2
        return z xor (z ushr 31)
    }

    /** A uniform integer in `0 until bound`, by rejection so that no value is likelier than another. */
    fun nextInt(bound: Int): Int {
        require(bound > 0) { "bound must be positive: $bound" }
        val limit = Long.MAX_VALUE - Long.MAX_VALUE % bound
        while (true) {
            val draw = nextLong() ushr 1
            if (draw < limit) return (draw % bound).toInt()
        }
    }

    /**
     * Shuffles [values] in place as `LinearModelUtils.shuffle` does (Fisher-Yates from the last element down, each one
     * swapped with a uniformly chosen element at or before it).
     */
    fun shuffle(values: IntArray) {
        for (i in values.size - 1 downTo 1) {
            val j = nextInt(i + 1)
            val kept = values[j]
            values[j] = values[i]
            values[i] = kept
        }
    }

    companion object {
        private val GOLDEN_GAMMA = 0x9E3779B97F4A7C15UL.toLong()
        private val MIX_1 = 0xBF58476D1CE4E5B9UL.toLong()
        private val MIX_2 = 0x94D049BB133111EBUL.toLong()
        private val STREAM = 0xD1B54A32D192ED03UL.toLong()

        /**
         * The generator of trait number [trait] under [seed]: its own stream, so a trait's permutations do not depend on
         * how many traits come before it or on how many workers fit the markers.
         */
        fun forTrait(seed: Long, trait: Int): PermutationRandom =
            PermutationRandom(PermutationRandom(seed xor ((trait + 1).toLong() * STREAM)).nextLong())
    }
}

/** The levels at which experiment-wise thresholds are reported. */
internal val PERMUTATION_ALPHAS: List<Double> = listOf(0.05, 0.01)

/**
 * What a permutation run learned, as it goes into the provenance of a result: the number of permutations and, per trait,
 * the experiment-wise threshold at each of [PERMUTATION_ALPHAS] (`NaN` where there are too few permutations to resolve
 * the level).
 */
internal data class PermutationSummary(val permutations: Int, val thresholds: Map<String, Map<Double, Double>>)

/**
 * The null distribution of a GLM scan: for each trait, the smallest marker p-value of each permutation, sorted. A
 * marker's `perm_p` is read from it, and so is the trait's threshold.
 *
 * @property traits Trait names, in analysis order.
 * @property minimumP One ascending array of `options.permutations` values per trait.
 */
internal class PermutationNull(
    val options: PermutationOptions,
    val traits: List<String>,
    val minimumP: Array<DoubleArray>,
) {
    init {
        require(minimumP.size == traits.size) { "One distribution per trait" }
        require(minimumP.all { it.size == options.permutations }) { "Each distribution has one value per permutation" }
    }

    /**
     * TASSEL 5's `perm_p` for a marker of trait [trait] with p-value [p]: `(index + 1) / nperm`, where `index` is where
     * `p` would be inserted in the sorted minimums, capped at 1. Where TASSEL 5's binary search lands on any one of
     * several equal minimums, this counts only the minimums strictly below `p`, so the two differ only for a `p` exactly
     * equal to a permutation's minimum. `NaN` for a `NaN` p.
     */
    fun permutationP(trait: Int, p: Double): Double {
        if (p.isNaN()) return Double.NaN
        val sorted = minimumP[trait]
        var low = 0
        var high = sorted.size
        while (low < high) {
            val middle = (low + high) ushr 1
            if (sorted[middle] < p) low = middle + 1 else high = middle
        }
        return minOf(1.0, (low + 1).toDouble() / sorted.size)
    }

    /**
     * The experiment-wise threshold of trait [trait] at [alpha]: the largest p-value whose `perm_p` is at most [alpha],
     * which is the `floor(alpha * nperm)`-th smallest minimum. `NaN` when that is fewer than one permutation.
     */
    fun threshold(trait: Int, alpha: Double): Double {
        require(alpha > 0.0 && alpha < 1.0) { "alpha must be between 0 and 1: $alpha" }
        val sorted = minimumP[trait]
        val k = floor(alpha * sorted.size + 1e-9).toInt()
        return if (k < 1) Double.NaN else sorted[k - 1]
    }

    fun summary(): PermutationSummary = PermutationSummary(
        options.permutations,
        traits.indices.associate { trait ->
            traits[trait] to PERMUTATION_ALPHAS.associateWith { threshold(trait, it) }
        },
    )
}

/**
 * Puts `perm_p` on the omnibus record (the first of each trait and marker, ADR-0027) of every chunk it passes on, from
 * [permutation]. Keeps no state, so chunks may arrive in any order.
 */
internal class PermutationPSink(private val permutation: PermutationNull, private val downstream: AssociationSink) :
    AssociationSink {
    private val traitIndex = permutation.traits.withIndex().associate { it.value to it.index }

    override suspend fun accept(chunk: AssociationChunk) {
        var groupTrait: String? = null
        var groupMarker = -1
        val records = chunk.records.map { record ->
            if (record.trait == groupTrait && record.markerIndex == groupMarker) {
                record
            } else {
                groupTrait = record.trait
                groupMarker = record.markerIndex
                val trait = traitIndex.getValue(record.trait)
                record.copy(permutationP = permutation.permutationP(trait, record.p))
            }
        }
        downstream.accept(AssociationChunk(chunk.sequence, records, chunk.effects))
    }
}
