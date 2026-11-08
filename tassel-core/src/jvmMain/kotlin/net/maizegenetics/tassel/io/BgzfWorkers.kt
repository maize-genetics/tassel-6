package net.maizegenetics.tassel.io

actual fun bgzfWorkerCount(): Int = Runtime.getRuntime().availableProcessors().coerceAtLeast(1)
