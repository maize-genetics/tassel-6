package net.maizegenetics.tassel.io

import kotlinx.coroutines.runBlocking

actual fun bgzfWorkerCount(): Int = Runtime.getRuntime().availableProcessors().coerceAtLeast(1)

internal actual fun <T> bgzfBlocking(block: suspend () -> T): T = runBlocking { block() }
