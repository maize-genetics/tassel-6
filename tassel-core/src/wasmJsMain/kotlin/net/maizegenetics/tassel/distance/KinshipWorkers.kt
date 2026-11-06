package net.maizegenetics.tassel.distance

/** The browser and Node have one wasm thread. Dispatchers.Default does not add cores here. */
internal actual fun kinshipWorkerCount(): Int = 1
