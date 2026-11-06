package net.maizegenetics.tassel.distance

internal actual fun kinshipWorkerCount(): Int = Runtime.getRuntime().availableProcessors().coerceAtLeast(1)
