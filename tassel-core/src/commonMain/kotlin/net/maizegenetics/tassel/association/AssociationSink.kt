package net.maizegenetics.tassel.association

import kotlinx.coroutines.sync.Mutex
import kotlinx.coroutines.sync.withLock

/**
 * The results of one unit of work, normally one site block for one trait. [sequence] numbers the chunks of a scan from
 * 0 in the order the report must list them; producers may finish in any order.
 *
 * Within a chunk, the records of one (trait, marker) are adjacent, omnibus test first.
 */
internal class AssociationChunk(
    val sequence: Int,
    val records: List<AssociationRecord>,
    val effects: List<AssociationEffect> = emptyList(),
)

/** Where kernels send results (ADR-0027). Sinks that wrap another sink pass chunks on in the order they get them. */
internal interface AssociationSink {
    suspend fun accept(chunk: AssociationChunk)
}

/**
 * Makes concurrent producers deterministic: chunks may arrive in any order from any number of coroutines, and
 * [downstream] sees them in sequence order, one at a time, so it need not be thread safe. Only chunks that arrive early
 * are buffered.
 */
internal class OrderedAssociationSink(private val downstream: AssociationSink) : AssociationSink {
    private val lock = Mutex()
    private val waiting = HashMap<Int, AssociationChunk>()
    private var next = 0

    override suspend fun accept(chunk: AssociationChunk) {
        lock.withLock {
            require(chunk.sequence >= next && chunk.sequence !in waiting) {
                "Chunk ${chunk.sequence} was already accepted"
            }
            waiting[chunk.sequence] = chunk
            while (true) {
                val ready = waiting.remove(next) ?: break
                next++
                downstream.accept(ready)
            }
        }
    }

    /** Fails if a chunk is missing, which means a producer stopped before it delivered. */
    suspend fun finish() {
        lock.withLock {
            check(waiting.isEmpty()) {
                "Chunk $next never arrived; ${waiting.size} later chunks are waiting"
            }
        }
    }
}

/** Which TASSEL 5 analysis's `maxP` rule applies. */
internal enum class MaxPRule {
    /** GLM keeps a row iff `maxP == 1.0 || p <= maxP`, so rows with no p are dropped when filtering. */
    Glm,

    /** MLM drops a row iff `p > maxP`, so rows with no p are kept. */
    Mlm,
    ;

    fun keeps(p: Double, maxP: Double): Boolean = when (this) {
        Glm -> maxP == 1.0 || p <= maxP
        Mlm -> !(p > maxP)
    }
}

/**
 * Keeps the (trait, marker) groups whose omnibus (first) record passes [rule]; the effects of a dropped group go with
 * it. A filter that keeps everything passes chunks through untouched.
 */
internal class MaxPSink(
    private val maxP: Double,
    private val rule: MaxPRule,
    private val downstream: AssociationSink,
) : AssociationSink {
    init {
        require(maxP in 0.0..1.0) { "maxP must be between 0 and 1: $maxP" }
    }

    private val passesAll = maxP == 1.0

    override suspend fun accept(chunk: AssociationChunk) {
        if (passesAll) return downstream.accept(chunk)
        val kept = ArrayList<AssociationRecord>(chunk.records.size)
        val dropped = HashSet<Pair<String, Int>>()
        var groupTrait: String? = null
        var groupMarker = -1
        var groupKept = true
        for (record in chunk.records) {
            if (record.trait != groupTrait || record.markerIndex != groupMarker) {
                groupTrait = record.trait
                groupMarker = record.markerIndex
                groupKept = rule.keeps(record.p, maxP)
                if (!groupKept) dropped.add(record.trait to record.markerIndex)
            }
            if (groupKept) kept.add(record)
        }
        val effects = if (dropped.isEmpty()) {
            chunk.effects
        } else {
            chunk.effects.filter { (it.trait to it.markerIndex) !in dropped }
        }
        downstream.accept(AssociationChunk(chunk.sequence, kept, effects))
    }
}

/** Counts every record, then passes the chunk on. Not thread safe; put it behind an [OrderedAssociationSink]. */
internal class CountingSink(private val downstream: AssociationSink) : AssociationSink {
    private var records = 0L
    private var tested = 0L
    private var failed = 0L
    private val skipped = LinkedHashMap<SkipReason, Long>()

    val summary: AssociationSummary
        get() = AssociationSummary(records, tested, skipped.toMap(), failed)

    override suspend fun accept(chunk: AssociationChunk) {
        for (record in chunk.records) {
            records++
            when (record.status) {
                AssociationStatus.Ok -> tested++
                AssociationStatus.Failed -> failed++
                AssociationStatus.Skipped -> record.skipReason!!.let { skipped[it] = (skipped[it] ?: 0L) + 1L }
            }
        }
        downstream.accept(chunk)
    }
}

/** An [AssociationTableBuilder] as a sink. */
internal class TableSink(val builder: AssociationTableBuilder = AssociationTableBuilder()) : AssociationSink {
    override suspend fun accept(chunk: AssociationChunk) {
        chunk.records.forEach(builder::append)
        chunk.effects.forEach(builder::append)
    }
}

/**
 * The default pipeline: chunks in any order, then in sequence order through a count of everything and the `maxP`
 * filter into an in-memory table. Call [finish] once every producer is done. [traits] fixes the order reports list traits in
 * (the order of the analysis), which a `maxP` filter that drops every row of the first trait would otherwise change.
 */
internal class AssociationCollector(
    maxP: Double = 1.0,
    rule: MaxPRule = MaxPRule.Glm,
    traits: List<String> = emptyList(),
) : AssociationSink {
    private val table = TableSink(AssociationTableBuilder().also { it.reserveTraits(traits) })
    private val counting = CountingSink(MaxPSink(maxP, rule, table))
    private val ordered = OrderedAssociationSink(counting)

    override suspend fun accept(chunk: AssociationChunk) = ordered.accept(chunk)

    suspend fun finish(provenance: AnalysisProvenance): AssociationTable {
        ordered.finish()
        return table.builder.build(provenance, counting.summary)
    }
}
