package net.maizegenetics.tassel.spike.gram

import kotlin.time.TimeSource

class RunStats(val minSeconds: Double, val medianSeconds: Double, val samplesSeconds: List<Double>)

/**
 * One discarded warm-up, then [runs] timed calls. The median rank matches `tassel dev kinship-bench`:
 * 1-based `int((n + 1) / 2)`.
 */
fun timeRuns(runs: Int = 3, block: () -> Unit): RunStats = timeRunsReturning(runs, block).first

/** Like [timeRuns], and keeps the value from the last timed call. */
fun <T : Any> timeRunsReturning(runs: Int = 3, block: () -> T): Pair<RunStats, T> {
    require(runs >= 1) { "runs must be positive: $runs" }
    block()
    val samples = DoubleArray(runs)
    lateinit var last: T
    repeat(runs) { index ->
        val started = TimeSource.Monotonic.markNow()
        last = block()
        samples[index] = started.elapsedNow().inWholeNanoseconds / 1e9
    }
    return statsOf(samples.toList()) to last
}

fun statsOf(samples: List<Double>): RunStats {
    require(samples.isNotEmpty()) { "no samples" }
    val sorted = samples.sorted()
    val medianAt = (sorted.size + 1) / 2 - 1
    return RunStats(sorted.first(), sorted[medianAt], samples)
}

fun secondsSince(mark: TimeSource.Monotonic.ValueTimeMark): Double = mark.elapsedNow().inWholeNanoseconds / 1e9
