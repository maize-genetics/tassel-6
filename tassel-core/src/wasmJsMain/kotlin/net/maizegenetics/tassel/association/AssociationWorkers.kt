package net.maizegenetics.tassel.association

/** The browser and Node have one wasm thread. Dispatchers.Default does not add cores here. */
internal actual fun associationWorkerCount(): Int = 1
