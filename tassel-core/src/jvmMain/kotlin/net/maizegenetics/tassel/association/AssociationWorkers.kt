package net.maizegenetics.tassel.association

internal actual fun associationWorkerCount(): Int = Runtime.getRuntime().availableProcessors().coerceAtLeast(1)
