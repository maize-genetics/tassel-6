package net.maizegenetics.tassel.io

/**
 * How many BGZF workers to use, for inflate and for deflate. One on wasmJs, where
 * [kotlinx.coroutines.Dispatchers.Default] is a single thread; the processor count on the JVM and
 * on Native.
 *
 * A sibling of `kinshipWorkerCount` (ADR-0020). `io` does not call into `distance`. Results do not
 * depend on this count (ADR-0021, ADR-0022). It is not an operation-spec parameter.
 */
expect fun bgzfWorkerCount(): Int

/**
 * Runs [block] on the caller thread until it finishes.
 *
 * The JVM and Native park here. wasmJs has no `runBlocking`; the one-worker writer never calls this.
 */
internal expect fun <T> bgzfBlocking(block: suspend () -> T): T

/** Blocks held at once: copied from the file, not yet returned to the caller. Twice [workers]. */
internal fun bgzfMaxBlocksInFlight(workers: Int): Int = workers.coerceAtLeast(1) * 2
