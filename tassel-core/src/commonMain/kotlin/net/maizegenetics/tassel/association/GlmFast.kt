package net.maizegenetics.tassel.association

import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.coroutineScope
import net.maizegenetics.tassel.math.Cholesky
import net.maizegenetics.tassel.math.CholeskyResult
import net.maizegenetics.tassel.math.DenseMatrix
import net.maizegenetics.tassel.math.HouseholderQr
import net.maizegenetics.tassel.math.Triangular
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.stats.model.AliasAwareQr
import net.maizegenetics.tassel.stats.model.DesignMatrix
import net.maizegenetics.tassel.stats.model.ModelSpec
import kotlin.math.abs
import kotlin.math.sqrt

/** How the fast kernel fits a marker that drops observations (missing genotypes) from its trait's set. */
internal enum class DroppedRows {
    /**
     * The default: indicator columns for the dropped rows on the set's cached base, in closed form ([GramBasis]), at
     * `O(r rank^2)` per marker; a pattern base where the dropped rows carry so much leverage that the closed form would
     * be ill-conditioned or the kept rows' base could lose rank.
     */
    Indicators,

    /** A base QR of the kept rows for every marker that drops rows, cached by missingness pattern. */
    PatternBase,
}

/**
 * Settings of [glm] that change how it works, not what it computes: every choice gives the same results to rounding,
 * and none depends on the number of workers.
 *
 * @property patternCacheDoubles Doubles all workers' missingness-pattern caches may hold together.
 */
internal data class GlmFastTuning(
    val droppedRows: DroppedRows = DroppedRows.Indicators,
    val patternCacheDoubles: Long = 4_000_000L,
) {
    init {
        require(patternCacheDoubles >= 0) { "patternCacheDoubles must not be negative: $patternCacheDoubles" }
    }
}

/** What [glm] did, for tests and benchmarks: one count per (observation set, marker). */
internal class GlmFastStats {
    /** Fitted on the set's own cached base: no observation dropped. */
    var fullSet = 0L

    /** Fitted with indicator columns for the dropped rows, in closed form. */
    var indicators = 0L

    /** Fitted on a pattern base found in the cache. */
    var patternHits = 0L

    /** Fitted on a pattern base built for the marker (and cached). */
    var patternMisses = 0L

    /** Replicated traits, fitted by the reference path. */
    var reference = 0L

    /** Skipped markers. */
    var skipped = 0L

    internal fun add(other: GlmFastStats) {
        fullSet += other.fullSet
        indicators += other.indicators
        patternHits += other.patternHits
        patternMisses += other.patternMisses
        reference += other.reference
        skipped += other.skipped
    }
}

/**
 * The GLM fast path (tassel-6#192): the same tests, records, and order as [glmReference], computed by residualizing
 * against cached null models (Frisch-Waugh-Lovell, ADR-0025) and spread over [workers].
 *
 * - Traits whose base observation sets are equal form a group. A group's base design `X0` (mean, factors, covariates)
 *   is factored once, and each trait's response is rotated by `Q0'` once.
 * - A marker's block of columns `M` is rotated by `Q0'` and its rows past the base rank, `M_r = Q0_perp' M`, get a small
 *   QR of their own. The block's sum of squares, the error, and the effects come from that QR applied to each trait's
 *   rotated response, so the block work is done once per group and only `O(n k)` is per trait. With no dropped rows
 *   this is the arithmetic the reference does on `[X0 | M]`, reorganized, so the sums of squares agree to rounding.
 * - A marker that drops `r` rows (missing calls) is fitted on its kept rows, by [DroppedRows]: either `r` indicator
 *   columns on the set's base (an exact identity, ADR-0025), applied in closed form from the dropped rows of the base's
 *   orthonormal columns ([GramBasis]), or a base QR of the kept rows from a per-worker LRU cache keyed by the
 *   missingness pattern and verified on hit. The choice depends on the dropped rows' leverage only, never on what a
 *   cache holds, so it cannot depend on the worker count.
 * - Replicated traits (taxa nested in the marker, decision Q1) are fitted by the reference path, marker by marker.
 *
 * Each site block is split into up to [workers] contiguous runs of markers, each encoded through its own loader
 * ([MarkerBlock.fork]) and fitted on a coroutine that inherits the caller's dispatcher. A marker is fitted by one worker
 * from inputs no other marker changes, and results go back in table order, so any [workers] gives the same bytes
 * (ADR-0020, ADR-0025). Progress and cancellation are per site block, as in the reference.
 *
 * @see net.maizegenetics.analysis.association.FastMultithreadedAssociationPlugin
 * @see net.maizegenetics.analysis.modelfitter.AdditiveResidualForwardRegression
 */
internal suspend fun glm(
    data: AssociationData,
    source: MarkerSource,
    options: GlmOptions,
    sink: AssociationSink,
    progress: ProgressReporter = ProgressReporter.None,
    workers: Int = associationWorkerCount(),
    tuning: GlmFastTuning = GlmFastTuning(),
    stats: GlmFastStats? = null,
) {
    require(workers > 0) { "workers must be positive: $workers" }
    val plan = GlmPlan(data, options, tuning)
    val states = Array(workers) { GlmWorker(plan, tuning.patternCacheDoubles / workers) }
    val total = source.numberOfMarkers.toLong()
    var done = 0L
    var sequence = 0
    progress.report(0L, total, "glm")
    for (block in source.scan(plan.groups.map { it.observations }, options.marker, options.blockSize)) {
        val first = block.markers.first
        val count = block.markers.last - first + 1
        val positions = Array(count) { source.position(first + it) }
        val outputs = arrayOfNulls<MarkerOutput>(count)
        runBands(block, count, states) { state, loader, band -> state.fit(loader, band, first, positions, outputs) }
        val records = ArrayList<AssociationRecord>()
        val effects = ArrayList<AssociationEffect>()
        for (output in outputs) output!!.drainInto(records, effects)
        sink.accept(AssociationChunk(sequence++, records, effects))
        done += count
        progress.report(done, total, "glm")
    }
    if (stats != null) states.forEach { stats.add(it.stats) }
}

/**
 * The permutation null of TASSEL 5's GLM (`AbstractFixedEffectLM.createPermutedData` and `updateMinP`, tassel-6#194): for
 * each trait, [PermutationOptions.permutations] arrangements of its base-model residuals, each refitted at every marker
 * that is tested, with the smallest marker p-value of each arrangement kept.
 *
 * The residuals are shuffled cumulatively (arrangement `k` is arrangement `k - 1` shuffled again) from a generator of
 * the trait's own stream ([PermutationRandom.forTrait]), so the arrangements depend on the seed and the trait's position
 * alone. A marker's fit is a permuted response in the same residualized basis the observed fit uses, over the marker's
 * kept observations, so missing calls, aliased columns, and every encoding behave as in [glm]. The permuted response is the
 * shuffled residual alone, without TASSEL 5's fitted values added back: they lie in the base model's column space,
 * which the marker test projects out, so the p-values are the same. The statistic is the marker's omnibus test (the
 * first record of [glm]); a marker with no test degrees of freedom has no p and leaves the minimum alone.
 *
 * Minimums are merged by the smallest value, which does not depend on the order or the number of workers, so any [workers]
 * gives the same distribution. Replicated taxa are refused: shuffling residuals is only exchangeable when each
 * observation has its own taxon.
 *
 * @see net.maizegenetics.analysis.association.AbstractFixedEffectLM
 * @see net.maizegenetics.stats.linearmodels.LinearModelUtils.shuffle
 */
internal suspend fun glmPermutationNull(
    data: AssociationData,
    source: MarkerSource,
    options: GlmOptions,
    permutation: PermutationOptions,
    progress: ProgressReporter = ProgressReporter.None,
    workers: Int = associationWorkerCount(),
    tuning: GlmFastTuning = GlmFastTuning(),
): PermutationNull {
    require(workers > 0) { "workers must be positive: $workers" }
    val plan = GlmPlan(data, options, tuning)
    val replicated = plan.groups.firstOrNull { it.replicated }
    require(replicated == null) {
        "glm: permutations are not defined when taxa are replicated " +
            "(trait '${data.traits[replicated!!.traits.first()]}' has more than one observation of a taxon)"
    }
    val draws = plan.draw(permutation)
    val states = Array(workers) { GlmWorker(plan, tuning.patternCacheDoubles / workers, draws) }
    val total = source.numberOfMarkers.toLong()
    var done = 0L
    progress.report(0L, total, "glm permutations")
    for (block in source.scan(plan.groups.map { it.observations }, options.marker, options.blockSize)) {
        val first = block.markers.first
        val count = block.markers.last - first + 1
        runBands(block, count, states) { state, loader, band -> state.permute(loader, band, first) }
        done += count
        progress.report(done, total, "glm permutations")
    }
    val minimum = Array(data.traits.size) { DoubleArray(permutation.permutations) { 1.0 } }
    for (state in states) {
        val mine = state.minimumP ?: continue
        for (trait in minimum.indices) {
            for (k in 0 until permutation.permutations) {
                if (mine[trait][k] <
                    minimum[trait][k]
                ) {
                    minimum[trait][k] = mine[trait][k]
                }
            }
        }
    }
    minimum.forEach { it.sort() }
    return PermutationNull(permutation, data.traits, minimum)
}

/** A trait's base-model residuals, and each permutation's arrangement of them: `arrangement[row]` is the residual at `row`. */
private class PermutedTrait(val residuals: DoubleArray, val arrangements: Array<IntArray>)

/** The permuted residuals of every trait that can be tested, null for a trait with no observations. */
private class PermutationDraws(val permutations: Int, val traits: Array<PermutedTrait?>)

/** Draws [permutation]'s arrangements for each trait of the plan. */
private fun GlmPlan.draw(permutation: PermutationOptions): PermutationDraws {
    val traits = arrayOfNulls<PermutedTrait>(data.traits.size)
    for (group in groups) {
        val full = group.full ?: continue
        val n = full.basis.rows
        val rank = full.basis.baseRank
        val q = full.q1.data
        group.traits.forEachIndexed { slot, trait ->
            // The residual of the base fit: y minus its projection on the orthonormal base columns Q1.
            val values = group.values[slot]
            val y = DoubleArray(n) { values.value(group.observations[it]) }
            val residuals = y.copyOf()
            for (j in 0 until rank) {
                var dot = 0.0
                for (i in 0 until n) dot += q[j * n + i] * y[i]
                for (i in 0 until n) residuals[i] -= dot * q[j * n + i]
            }
            val random = PermutationRandom.forTrait(permutation.seed, trait)
            val arrangement = IntArray(n) { it }
            traits[trait] = PermutedTrait(
                residuals,
                Array(permutation.permutations) { arrangement.also(random::shuffle).copyOf() },
            )
        }
    }
    return PermutationDraws(permutation.permutations, traits)
}

/**
 * Runs [work] on up to `states.size` contiguous runs of the [count] markers of [block], each with its own worker state and
 * loader. Inherits the caller's dispatcher: the JVM and native engines run the bands on `Dispatchers.Default`, a test
 * dispatcher and the browser's one thread run them one after another.
 */
private suspend fun runBands(
    block: MarkerBlock,
    count: Int,
    states: Array<GlmWorker>,
    work: (GlmWorker, MarkerBlock, IntRange) -> Unit,
) {
    val bands = bands(count, states.size)
    if (bands.size == 1) {
        work(states[0], block, bands[0])
        return
    }
    val loaders = List(bands.size) { if (it == 0) block else block.fork() }
    coroutineScope {
        bands.mapIndexed { w, band -> async { work(states[w], loaders[w], band) } }.awaitAll()
    }
}

/** `0 until count` in up to [workers] contiguous, nearly equal runs, in order. */
private fun bands(count: Int, workers: Int): List<IntRange> {
    val n = minOf(count, workers).coerceAtLeast(1)
    return (0 until n).map { (count * it / n) until (count * (it + 1) / n) }
}

/** The records and effects of one marker, by trait, put back in trait order by [drainInto]. */
private class MarkerOutput(traits: Int) {
    val records = Array(traits) { ArrayList<AssociationRecord>(3) }
    val effects = Array(traits) { ArrayList<AssociationEffect>(0) }

    fun drainInto(records: MutableList<AssociationRecord>, effects: MutableList<AssociationEffect>) {
        this.records.forEach { records.addAll(it) }
        this.effects.forEach { effects.addAll(it) }
    }
}

/** Traits that share a base observation set, and therefore the base design and every marker's columns. */
private class GlmGroup(
    val index: Int,
    /** Trait indices, ascending. */
    val traits: IntArray,
    /** The shared base observation set. */
    val observations: IntArray,
    val replicated: Boolean,
    /** The base model; its response is the first trait, and the design is the same for every trait of the group. */
    val spec: ModelSpec,
    val values: List<NumericAttribute>,
) {
    /** The base over the whole set, or null for a replicated or empty set. Built once, then only read. */
    var full: FullSet? = null
}

/**
 * A group's base over its whole set, and what [GramBasis] needs from it: the orthonormal columns `Q1` (`n x rank`), each
 * row's leverage (its squared norm), and whether every base column is far enough from the earlier ones that no subset
 * of rows the closed form accepts can make it aliased.
 */
private class FullSet(val basis: QrBasis) {
    val q1: DenseMatrix = DenseMatrix(basis.rows, basis.baseRank).also { q ->
        for (j in 0 until basis.baseRank) q.data[j * basis.rows + j] = 1.0
        basis.qr.applyQ(q)
    }

    val leverage = DoubleArray(basis.rows) { i ->
        var sum = 0.0
        for (j in 0 until basis.baseRank) sum += q1.data[j * basis.rows + i] * q1.data[j * basis.rows + i]
        sum
    }

    /**
     * Every diagonal of the base's QR (its columns at unit length) is at least [SEPARATION] times the aliasing
     * tolerance: on kept rows whose `Q1` rows have smallest singular value `s`, a column's distance from the earlier
     * ones shrinks by at most `s`, and the closed form only takes rows with `s^2 >= GRAM_MIN_EIGENVALUE`.
     */
    val separated: Boolean = basis.qr.rDiagonal.all {
        abs(it) >= SEPARATION * GLM_RANK_TOLERANCE * basis.r11 / sqrt(GRAM_MIN_EIGENVALUE)
    }

    /** Whether rows [dropped] leave `G = Q1_K' Q1_K` with every eigenvalue at least [GRAM_MIN_EIGENVALUE]. */
    fun acceptsGram(dropped: IntArray): Boolean {
        if (!separated) return false
        // The eigenvalues of G are 1 minus those of Q1_D' Q1_D, whose largest is at most its trace: the leverage sum.
        var sum = 0.0
        for (row in dropped) sum += leverage[row]
        return sum <= 1.0 - GRAM_MIN_EIGENVALUE
    }

    private companion object {
        const val SEPARATION = 4.0
    }
}

/** The smallest eigenvalue [GramBasis] accepts for `G`: it bounds how much the closed form amplifies rounding. */
private const val GRAM_MIN_EIGENVALUE = 1e-3

/** The groups of a run and their bases over the whole set, shared read-only by the workers. */
private class GlmPlan(val data: AssociationData, val options: GlmOptions, val tuning: GlmFastTuning) {
    val groups: List<GlmGroup>

    init {
        val bySet = LinkedHashMap<SetKey, MutableList<Int>>()
        for (trait in data.traits.indices) bySet.getOrPut(SetKey(data.observations(trait))) { ArrayList() } += trait
        groups = bySet.values.mapIndexed { index, traits ->
            val first = traits.first()
            GlmGroup(
                index,
                traits.toIntArray(),
                data.observations(first),
                data.hasReplicates(first),
                data.baseSpec(first),
                traits.map { data.data.phenotype.attribute(data.traits[it]) as NumericAttribute },
            )
        }
        for (group in groups) {
            if (!group.replicated && group.observations.isNotEmpty()) {
                group.full = FullSet(qrBasis(group, group.observations))
            }
        }
    }

    /**
     * The base of [group] over [observations] (phenotype rows), with every trait's response rotated by it: the design
     * and QR the reference builds for those rows.
     */
    fun qrBasis(group: GlmGroup, observations: IntArray): QrBasis {
        val design = DesignMatrix.build(group.spec, data.data.phenotype, observations, GLM_RANK_TOLERANCE)
        val decomposition = design.decomposition
        // The intercept is a nonzero column, so with any observation the base has rank at least 1.
        val qr = checkNotNull(decomposition.qr) { "The base design of '${group.spec.response}' has no column" }
        val basis = QrBasis(qr, design.rows, decomposition.rank)
        basis.responses = Array(group.traits.size) { slot ->
            val values = group.values[slot]
            basis.response(DoubleArray(observations.size) { values.value(observations[it]) })
        }
        return basis
    }

    private class SetKey(val set: IntArray) {
        private val hash = set.contentHashCode()

        override fun hashCode(): Int = hash

        override fun equals(other: Any?): Boolean = other is SetKey && other.set.contentEquals(set)
    }
}

/** A trait's response in a basis: what the base explains past the mean, and the residual coordinates. */
private class Response(
    /** The response with the base projected out, [Basis.length] values. */
    val residual: DoubleArray,
    /** The base model's sum of squares past the mean, summed in the order the reference sums `Q'y` where it can be. */
    val baseSs: Double,
)

/**
 * The base model of one marker's kept observations: what a marker block is residualized against. [rows] is the number
 * of kept observations and [baseRank] the rank of the base design over them; a residual has [length] values and lies in
 * a space of [dimension] `= rows - baseRank`, which is what is left for the marker and the error.
 */
private sealed class Basis {
    abstract val rows: Int
    abstract val baseRank: Int
    abstract val length: Int
    val dimension: Int get() = rows - baseRank

    /** The first diagonal of the base QR: the scale of the reference's aliasing tolerance. */
    abstract val r11: Double

    /** One per trait of the group, in its order. */
    lateinit var responses: Array<Response>

    /** The residuals (`length x k`) of [block], `rows x k` with one row per kept observation, which it may overwrite. */
    abstract fun residualize(block: DenseMatrix): DenseMatrix

    /** A response over the kept observations as [responses] holds them: [y] has [rows] values, which this may overwrite. */
    abstract fun response(y: DoubleArray): Response
}

/**
 * A base factored over exactly the kept observations: the set's own, or one built for a missingness pattern. Its QR is
 * the one the reference builds for those rows, and residuals are its rotated coordinates past the base rank, so a fit
 * repeats the reference's arithmetic on `[X0 | M]`.
 */
private class QrBasis(val qr: HouseholderQr, override val rows: Int, override val baseRank: Int) : Basis() {
    override val length: Int get() = rows - baseRank
    override val r11: Double = abs(qr.rDiagonal[0])

    /** Doubles held: the factor and the responses, for the pattern cache's budget. */
    val size: Long get() = rows.toLong() * baseRank + responses.size.toLong() * length

    override fun residualize(block: DenseMatrix): DenseMatrix {
        qr.applyQt(block)
        return block.rowsFrom(baseRank)
    }

    /** Rotates the raw response [y], as the reference rotates `y` (not centered, so the sums repeat its order). */
    override fun response(y: DoubleArray): Response {
        qr.applyQt(y)
        return Response(y.copyOfRange(baseRank, rows), sumOfSquares(y, 1, baseRank))
    }
}

/**
 * The set's base with an indicator column per dropped row (ADR-0025): `[X0 | E]` over the whole set fits the kept rows
 * exactly, since each dropped row gets a parameter of its own. On the kept rows its residual projection is
 * `I - Q1_K G^-1 Q1_K'`, where `Q1_K` is the kept rows of the set base's orthonormal columns and
 * `G = Q1_K' Q1_K = I - Q1_D' Q1_D` is built from the dropped rows alone. So a marker costs `O(r rank^2 + rank^3)` and
 * then `O(n rank)` per column, against `O(n rank^2)` for a new QR of the kept rows. Residuals stay in kept-row
 * coordinates (length [rows], dimension `rows - rank`). Only built when [FullSet.acceptsGram] holds, so `G` is well
 * conditioned and the kept rows' base keeps the full rank.
 */
private class GramBasis(
    /** `Q1_K`, the kept rows of the set base's orthonormal columns, gathered once for the marker, column by column. */
    private val q1k: DoubleArray,
    override val rows: Int,
    override val baseRank: Int,
    /** The Cholesky factor of `G`. */
    private val lower: DenseMatrix,
    override val r11: Double,
) : Basis() {
    override val length: Int get() = rows

    override fun residualize(block: DenseMatrix): DenseMatrix {
        val w = DoubleArray(baseRank)
        for (j in 0 until block.columns) project(block.data, j * rows, w)
        return block
    }

    /** Centers [y] (one value per kept row) over the kept rows, and projects the base out of it. */
    override fun response(y: DoubleArray): Response {
        var sum = 0.0
        for (v in y) sum += v
        var mean = sum / y.size
        var correction = 0.0
        for (v in y) correction += v - mean
        mean += correction / y.size
        for (i in y.indices) y[i] -= mean
        val baseSs = project(y, 0, DoubleArray(baseRank))
        // A base of rank 1 is the mean alone and explains nothing past it: zero, as the reference's empty sum is.
        return Response(y, if (baseRank <= 1) 0.0 else baseSs)
    }

    /**
     * Replaces the [rows] values of [v] at [offset] with their residual from the base, and returns the squared length
     * of the projection removed, `c' G^-1 c = |L^-1 c|^2` with `c = Q1_K' v`, as a sum of squares rather than a
     * difference.
     */
    private fun project(v: DoubleArray, offset: Int, w: DoubleArray): Double {
        val rank = baseRank
        val n = rows
        val q = q1k
        // c = Q1_K' v, four columns per pass over v so that the sums are independent chains.
        var j = 0
        while (j < rank) {
            val o0 = j * n
            when (rank - j) {
                1 -> {
                    var s0 = 0.0
                    for (i in 0 until n) s0 += q[o0 + i] * v[offset + i]
                    w[j] = s0
                    j += 1
                }
                2 -> {
                    var s0 = 0.0
                    var s1 = 0.0
                    for (i in 0 until n) {
                        val x = v[offset + i]
                        s0 += q[o0 + i] * x
                        s1 += q[o0 + n + i] * x
                    }
                    w[j] = s0
                    w[j + 1] = s1
                    j += 2
                }
                else -> {
                    var s0 = 0.0
                    var s1 = 0.0
                    var s2 = 0.0
                    for (i in 0 until n) {
                        val x = v[offset + i]
                        s0 += q[o0 + i] * x
                        s1 += q[o0 + n + i] * x
                        s2 += q[o0 + 2 * n + i] * x
                    }
                    w[j] = s0
                    w[j + 1] = s1
                    w[j + 2] = s2
                    j += 3
                }
            }
        }
        Triangular.solveLower(lower, rank, w)
        val projected = sumOfSquares(w, 0, rank)
        Triangular.solveLowerTransposed(lower, rank, w)
        // v -= Q1_K w, three columns per pass.
        j = 0
        while (j < rank) {
            val o0 = j * n
            when (rank - j) {
                1 -> {
                    val w0 = w[j]
                    for (i in 0 until n) v[offset + i] -= w0 * q[o0 + i]
                    j += 1
                }
                2 -> {
                    val w0 = w[j]
                    val w1 = w[j + 1]
                    for (i in 0 until n) v[offset + i] -= w0 * q[o0 + i] + w1 * q[o0 + n + i]
                    j += 2
                }
                else -> {
                    val w0 = w[j]
                    val w1 = w[j + 1]
                    val w2 = w[j + 2]
                    for (i in 0 until n) {
                        v[offset + i] -= w0 * q[o0 + i] + w1 * q[o0 + n + i] + w2 * q[o0 + 2 * n + i]
                    }
                    j += 3
                }
            }
        }
        return projected
    }

    companion object {
        /**
         * The basis for [kept] and [dropped] rows (positions in the set), or null if `G` is not numerically positive
         * definite.
         */
        fun of(full: FullSet, kept: IntArray, dropped: IntArray): GramBasis? {
            val rank = full.basis.baseRank
            val n = full.basis.rows
            val q = full.q1.data
            val g = DenseMatrix(rank, rank)
            for (j in 0 until rank) g.data[j * rank + j] = 1.0
            // The lower triangle is what the Cholesky factorization reads: column a, rows b >= a.
            for (row in dropped) {
                for (a in 0 until rank) {
                    val qa = q[a * n + row]
                    for (b in a until rank) g.data[a * rank + b] -= qa * q[b * n + row]
                }
            }
            val factor = (Cholesky.of(g) as? CholeskyResult.Success)?.factor ?: return null
            val q1k = DoubleArray(kept.size * rank)
            for (j in 0 until rank) {
                for (i in kept.indices) q1k[j * kept.size + i] = q[j * n + kept[i]]
            }
            return GramBasis(q1k, kept.size, rank, factor.lower(), full.basis.r11)
        }
    }
}

/** A QR of the columns of a residualized matrix that are not aliased, in column order, and which those are. */
private class Factored(val kept: IntArray, val qr: HouseholderQr?) {
    val rank: Int get() = kept.size
}

/**
 * Factors the [candidates] columns of [residualized] in order, leaving out each column whose diagonal is not above
 * [threshold] (and any past [dimension], the room left in the residual space) and refactoring, as `AliasAwareQr` does
 * on the whole design. With [threshold] `GLM_RANK_TOLERANCE * r11` this is the reference's aliasing rule, so the same
 * marker columns are aliased.
 */
private fun factorInOrder(
    residualized: DenseMatrix,
    candidates: IntArray,
    threshold: Double,
    dimension: Int,
): Factored {
    val kept = candidates.toMutableList()
    val length = residualized.rows
    while (kept.isNotEmpty() && dimension > 0) {
        val matrix = if (kept.size == residualized.columns) {
            residualized
        } else {
            DenseMatrix(length, kept.size).also { matrix ->
                kept.forEachIndexed { target, source ->
                    residualized.data.copyInto(matrix.data, target * length, source * length, source * length + length)
                }
            }
        }
        // A zero relative tolerance only stops at an exactly zero diagonal; the threshold is applied here.
        val qr = HouseholderQr.of(matrix, pivot = false, relativeTolerance = 0.0)
        val diagonal = qr.rDiagonal
        var failing = (0 until minOf(diagonal.size, dimension)).firstOrNull { !(abs(diagonal[it]) > threshold) } ?: -1
        if (failing < 0 && kept.size > dimension) failing = dimension
        if (failing < 0) return Factored(kept.toIntArray(), qr)
        kept.removeAt(failing)
    }
    return Factored(IntArray(0), null)
}

/** A marker block of one group at one marker, residualized and factored once for all the group's traits. */
private class FactoredBlock(
    /** Columns in the block, before aliasing. */
    val columns: Int,
    /** `1 / norm` of each column over the kept rows; zero for an all-zero column. */
    val scale: DoubleArray,
    /** Kept columns, as indices into the block. */
    val factored: Factored,
) {
    /**
     * The fit of this block for [response]. The sums are taken in the reference's order: the base part, then the
     * block's entries, then the error's, so with a [QrBasis] they repeat its `Q'y` sums.
     */
    fun fit(basis: Basis, response: Response): GlmBlockFit {
        val z = response.residual.copyOf()
        val qr = factored.qr
        qr?.applyQt(z)
        val k = factored.rank
        var ss = 0.0
        var model = response.baseSs
        for (i in 0 until k) {
            ss += z[i] * z[i]
            model += z[i] * z[i]
        }
        var error = 0.0
        var total = model
        for (i in k until z.size) {
            error += z[i] * z[i]
            total += z[i] * z[i]
        }
        val errorDf = basis.dimension - k
        val estimate = {
            val coefficients = DoubleArray(columns) { Double.NaN }
            val standardErrors = DoubleArray(columns) { Double.NaN }
            if (qr != null) {
                val solved = qr.solveRotated(z)
                val rowNorms = qr.rInverseRowNorms()
                val sigma = if (errorDf > 0) sqrt(error / errorDf) else Double.NaN
                factored.kept.forEachIndexed { index, column ->
                    coefficients[column] = solved[index] * scale[column]
                    standardErrors[column] = sigma * rowNorms[index] * scale[column]
                }
            }
            GlmEstimates(coefficients, standardErrors, errorDf)
        }
        return GlmBlockFit(ss, k, error, errorDf, model, basis.baseRank - 1 + k, total, estimate)
    }

    companion object {
        /** Scales [values] (each with at least `basis.rows` values) as the reference does, residualizes, and factors. */
        fun of(basis: Basis, values: List<DoubleArray>): FactoredBlock {
            val n = basis.rows
            val scale = DoubleArray(values.size) {
                val norm = AliasAwareQr.columnNorm(values[it], 0, n)
                if (norm > 0.0 && norm.isFinite()) 1.0 / norm else 0.0
            }
            val candidates = values.indices.filter { scale[it] != 0.0 }.toIntArray()
            if (candidates.isEmpty()) return FactoredBlock(values.size, scale, Factored(IntArray(0), null))
            val scaled = DenseMatrix(n, candidates.size)
            candidates.forEachIndexed { target, column ->
                val source = values[column]
                val s = scale[column]
                for (i in 0 until n) scaled.data[target * n + i] = source[i] * s
            }
            val residualized = basis.residualize(scaled)
            val inOrder = factorInOrder(
                residualized,
                IntArray(candidates.size) { it },
                GLM_RANK_TOLERANCE * basis.r11,
                basis.dimension,
            )
            val kept = IntArray(inOrder.rank) { candidates[inOrder.kept[it]] }
            return FactoredBlock(values.size, scale, Factored(kept, inOrder.qr))
        }
    }
}

/**
 * One worker's state: its encoded columns, its missingness-pattern cache, a reference fitter for replicated traits, and
 * its counts. Encodes and fits the markers of a band and writes each marker's output to its own slot.
 */
private class GlmWorker(private val plan: GlmPlan, cacheDoubles: Long, private val draws: PermutationDraws? = null) {
    private val data = plan.data
    private val encoding = plan.options.marker.encoding
    private val builder = GlmRecords(plan.options)
    private val cache = PatternCache<QrBasis>(cacheDoubles)
    private val reference by lazy(LazyThreadSafetyMode.NONE) { ReferenceGlm(data, plan.options) }
    private val columns = Array(plan.groups.size) { MarkerColumns(plan.groups[it].observations.size) }
    private val loaded = BooleanArray(plan.groups.size)
    val stats = GlmFastStats()

    /** Per trait, the smallest marker p-value of each permutation among this worker's markers; null without [draws]. */
    val minimumP: Array<DoubleArray>? = draws?.let { d ->
        Array(plan.data.traits.size) { DoubleArray(d.permutations) { 1.0 } }
    }

    /** Encodes and fits markers `first + i` for `i` in [band], through [loader], into `outputs[i]`. */
    fun fit(
        loader: MarkerBlock,
        band: IntRange,
        first: Int,
        positions: Array<Position>,
        outputs: Array<MarkerOutput?>,
    ) {
        for (i in band) {
            for (group in plan.groups) loaded[group.index] = loader.load(first + i, group.index, columns[group.index])
            outputs[i] = analyze(first + i, positions[i])
        }
    }

    /** Encodes markers `first + i` for `i` in [band] through [loader] and folds their permuted fits into [minimumP]. */
    fun permute(loader: MarkerBlock, band: IntRange, first: Int) {
        for (i in band) {
            for (group in plan.groups) loaded[group.index] = loader.load(first + i, group.index, columns[group.index])
            for (group in plan.groups) {
                val markerColumns = columns[group.index]
                // A marker that is not encoded has no marker degrees of freedom, so no p-value.
                if (group.full == null || !loaded[group.index]) continue
                permuteGroup(group, markerColumns)
            }
        }
    }

    private fun permuteGroup(group: GlmGroup, markerColumns: MarkerColumns) {
        val draws = checkNotNull(draws)
        val basis = basis(group, markerColumns)
        val block = FactoredBlock.of(basis, blockColumns(encoding, GlmBlock.Marker, markerColumns))
        val n = markerColumns.size
        val kept = markerColumns.observation
        val y = DoubleArray(n)
        for (trait in group.traits) {
            val drawn = checkNotNull(draws.traits[trait])
            val minimum = checkNotNull(minimumP)[trait]
            for (k in 0 until draws.permutations) {
                val arrangement = drawn.arrangements[k]
                for (i in 0 until n) y[i] = drawn.residuals[arrangement[kept[i]]]
                val fit = block.fit(basis, basis.response(y))
                val p = Fstat.pValue(fit.ss, fit.df, fit.errorSs, fit.errorDf)
                if (p < minimum[k]) minimum[k] = p
            }
        }
    }

    private fun analyze(marker: Int, position: Position): MarkerOutput {
        val out = MarkerOutput(data.traits.size)
        for (group in plan.groups) {
            val markerColumns = columns[group.index]
            val isLoaded = loaded[group.index]
            if (group.replicated) {
                stats.reference++
                for (trait in group.traits) {
                    reference.analyze(
                        trait,
                        marker,
                        position,
                        markerColumns,
                        isLoaded,
                        out.records[trait],
                        out.effects[trait],
                    )
                }
                continue
            }
            if (!builder.isTested(markerColumns, isLoaded)) {
                stats.skipped++
                for (trait in group.traits) {
                    out.records[trait] +=
                        builder.skipped(GlmMarkerId(data.traits[trait], marker, position), markerColumns)
                }
                continue
            }
            val basis = basis(group, markerColumns)
            val blocks = arrayOfNulls<FactoredBlock>(GlmBlock.entries.size)
            group.traits.forEachIndexed { slot, trait ->
                val response = basis.responses[slot]
                val fits = GlmBlockFits { block ->
                    val factored = blocks[block.ordinal]
                        ?: FactoredBlock.of(basis, blockColumns(encoding, block, markerColumns))
                            .also { blocks[block.ordinal] = it }
                    factored.fit(basis, response)
                }
                val id = GlmMarkerId(data.traits[trait], marker, position)
                builder.tested(
                    id,
                    markerColumns,
                    isLoaded,
                    markerColumns.size,
                    "",
                    fits,
                    out.records[trait],
                    out.effects[trait],
                )
            }
        }
        return out
    }

    /** The base of the marker's kept observations: the set's own, indicator columns on it, or a pattern base. */
    private fun basis(group: GlmGroup, columns: MarkerColumns): Basis {
        val full = checkNotNull(group.full)
        val n = columns.size
        val setSize = group.observations.size
        if (n == setSize) {
            stats.fullSet++
            return full.basis
        }
        val kept = columns.observation.copyOf(n)
        val dropped = IntArray(setSize - n)
        var next = 0
        var d = 0
        for (position in 0 until setSize) {
            if (next < n && kept[next] == position) next++ else dropped[d++] = position
        }
        if (plan.tuning.droppedRows == DroppedRows.Indicators && full.acceptsGram(dropped)) {
            val basis = GramBasis.of(full, kept, dropped)
            if (basis != null) {
                stats.indicators++
                basis.responses = Array(group.traits.size) { slot ->
                    val values = group.values[slot]
                    basis.response(DoubleArray(n) { values.value(group.observations[kept[it]]) })
                }
                return basis
            }
        }
        val before = cache.misses
        val basis = cache.getOrPut(group.index, dropped, { it.size }) {
            plan.qrBasis(group, IntArray(n) { group.observations[kept[it]] })
        }
        if (cache.misses > before) stats.patternMisses++ else stats.patternHits++
        return basis
    }
}

/**
 * An LRU map from (observation set, missingness pattern) to a value, holding at most [budget] doubles by [getOrPut]'s
 * size function and at most [maxEntries] values. The key is the pattern's content: a lookup hashes it and an equal
 * hash is verified element by element before it counts as a hit. Not thread safe; each worker has its own.
 */
internal class PatternCache<V : Any>(private val budget: Long, private val maxEntries: Int = 1024) {
    private class Key(val set: Int, val pattern: IntArray) {
        private val hash = 31 * set + pattern.contentHashCode()

        override fun hashCode(): Int = hash

        override fun equals(other: Any?): Boolean =
            other is Key && other.set == set && other.pattern.contentEquals(pattern)
    }

    private class Entry<V>(val value: V, val size: Long)

    private val entries = LinkedHashMap<Key, Entry<V>>()
    private var held = 0L

    var hits = 0L
        private set
    var misses = 0L
        private set

    /** Values held now. */
    val size: Int get() = entries.size

    /** The value for [pattern] of [set], built by [build] on a miss and kept if it fits. [pattern] is not copied. */
    fun getOrPut(set: Int, pattern: IntArray, sizeOf: (V) -> Long, build: () -> V): V {
        val key = Key(set, pattern)
        val found = entries.remove(key)
        if (found != null) {
            hits++
            entries[key] = found
            return found.value
        }
        misses++
        val value = build()
        val size = sizeOf(value)
        if (size > budget) return value
        while (entries.isNotEmpty() && (held + size > budget || entries.size >= maxEntries)) {
            val eldest = entries.keys.iterator().next()
            held -= entries.remove(eldest)!!.size
        }
        entries[key] = Entry(value, size)
        held += size
        return value
    }
}

/** Rows [from] until the end of each column, as a new matrix. */
private fun DenseMatrix.rowsFrom(from: Int): DenseMatrix {
    val length = rows - from
    val result = DenseMatrix(length, columns)
    for (j in 0 until columns) data.copyInto(result.data, j * length, j * rows + from, j * rows + rows)
    return result
}

/** The sum of squares of [values] from [from] until [to], in order. */
private fun sumOfSquares(values: DoubleArray, from: Int, to: Int): Double {
    var sum = 0.0
    for (i in from until to) sum += values[i] * values[i]
    return sum
}
