package net.maizegenetics.tassel.io

/**
 * Opens a file through this target's [openRandomAccessSource]. Declared here because that function is not an
 * `expect`: `wasmJsMain` must not contain the Node implementation (ADR-0021).
 */
internal expect fun openRandomAccessFile(absolutePath: String): RandomAccessSource

/** A sparse file of [size] bytes with [writes] placed at absolute offsets. The hole is not filled. */
internal expect fun createSparseFixture(size: Long, writes: List<Pair<Long, ByteArray>>): String

internal expect fun deleteFixture(path: String)
