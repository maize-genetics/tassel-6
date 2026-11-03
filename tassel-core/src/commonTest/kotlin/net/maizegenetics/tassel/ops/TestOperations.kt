package net.maizegenetics.tassel.ops

import kotlinx.coroutines.yield

enum class Rounding { NONE, FLOOR, CEILING }

/** Dummy operation using every parameter type; scales each input value. */
val scaleValues = operation<List<Double>, List<Double>>(
    id = "scale-values",
    legacyPlugin = "ScaleValuesPlugin",
    help = "Multiplies each value by a factor",
) {
    val factor by double(default = 1.0, range = 0.0..10.0, legacy = "scaleFactor", help = "Multiplier")
    val repeats by int(default = 1, range = 1..5, legacy = "numRepeats", help = "Times to apply the factor")
    val limit by int(range = 0..1000, legacy = "maxValues", help = "Keep at most this many values")
    val negate by boolean(legacy = "negateValues", help = "Negate the result")
    val rounding by enum(default = Rounding.NONE, legacy = "roundingMode", help = "Rounding")
    val report by file(kind = FileKind.Output, legacy = "reportFile", help = "Optional report")
    run { input, progress ->
        var values = input.take(limit() ?: input.size)
        repeat(repeats()) { i ->
            values = values.map { it * factor() }
            progress.report(i + 1, repeats(), report())
        }
        values.map { v ->
            val signed = if (negate()) -v else v
            when (rounding()) {
                Rounding.NONE -> signed
                Rounding.FLOOR -> kotlin.math.floor(signed)
                Rounding.CEILING -> kotlin.math.ceil(signed)
            }
        }
    }
}

/** Dummy operation with no input and a required file parameter. */
val readLines = operation<Unit, String>(id = "read-lines", legacyPlugin = "ReadLinesPlugin") {
    val path by file(legacy = "inputFile", help = "File to read").required()
    run { _, _ -> "read ${path()}" }
}

/** Counts to [Int.MAX_VALUE] with a suspending loop, so in practice it only stops when cancelled. */
val spin = operation<Unit, Int>(id = "spin") {
    run { _, progress ->
        var count = 0
        repeat(Int.MAX_VALUE) {
            progress.report(++count, 0)
            yield()
        }
        count
    }
}

/** Counts with a plain, non-suspending loop that only reports progress. */
val countBlocking = operation<Unit, Int>(id = "count-blocking") {
    val total by int(default = 10, range = 1..1000)
    run { _, progress -> countTo(total(), progress) }
}

private fun countTo(total: Int, progress: ProgressReporter): Int {
    for (i in 1..total) progress.report(i, total)
    return total
}
