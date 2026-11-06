package net.maizegenetics.tassel.distance

/**
 * How many workers kinship should use. One on wasmJs, where the dispatcher is a single thread;
 * the processor count on the JVM and on Native. Results do not depend on this (ADR-0020): bands
 * partition the triangle, they do not change which terms a cell adds.
 */
internal expect fun kinshipWorkerCount(): Int

/** Pairs in [rows] of an upper triangle over [numberOfTaxa] taxa, including the diagonal. */
internal fun pairCount(numberOfTaxa: Int, rows: IntRange): Int {
    var count = 0
    for (row in rows) count += numberOfTaxa - row
    return count
}

/**
 * Splits `0 until numberOfTaxa` into [workers] contiguous row bands whose pair counts are as even
 * as the row sizes allow. Row `first` owns `numberOfTaxa - first` pairs, so equal row counts are
 * not equal work. The bands are in row order and cover every row once.
 */
internal fun pairBands(numberOfTaxa: Int, workers: Int): List<IntRange> {
    if (numberOfTaxa <= 0) return emptyList()
    val count = workers.coerceIn(1, numberOfTaxa)
    if (count == 1) return listOf(0 until numberOfTaxa)
    val target = numberOfTaxa.toLong() * (numberOfTaxa + 1) / 2 / count
    val bands = ArrayList<IntRange>(count)
    var start = 0
    for (band in 0 until count - 1) {
        val rowsLeftForLaterBands = count - band - 1
        val lastAllowedEnd = numberOfTaxa - rowsLeftForLaterBands
        var end = start
        var pairs = 0L
        while (end < lastAllowedEnd) {
            val rowPairs = (numberOfTaxa - end).toLong()
            if (pairs > 0 && pairs + rowPairs > target) break
            pairs += rowPairs
            end++
            if (pairs >= target) break
        }
        if (end == start) end = start + 1
        bands += start until end
        start = end
    }
    bands += start until numberOfTaxa
    return bands
}
