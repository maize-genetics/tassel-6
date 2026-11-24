package net.maizegenetics.tassel.io

/** The browser and Node have one wasm thread. [kotlinx.coroutines.Dispatchers.Default] does not add cores. */
actual fun bgzfWorkerCount(): Int = 1

@Suppress("UNUSED_PARAMETER")
internal actual fun <T> bgzfBlocking(block: suspend () -> T): T = error("BGZF compression does not block on wasmJs")
