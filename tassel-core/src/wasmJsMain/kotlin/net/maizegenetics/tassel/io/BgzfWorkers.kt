package net.maizegenetics.tassel.io

/** The browser and Node have one wasm thread. [kotlinx.coroutines.Dispatchers.Default] does not add cores. */
actual fun bgzfWorkerCount(): Int = 1
