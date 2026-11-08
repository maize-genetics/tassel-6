@file:OptIn(ExperimentalWasmJsInterop::class)

package net.maizegenetics.tassel.io

internal actual fun openRandomAccessFile(absolutePath: String): RandomAccessSource =
    openRandomAccessSource(absolutePath)

internal actual fun createSparseFixture(size: Long, writes: List<Pair<Long, ByteArray>>): String {
    val path = nodeTempPath()
    val fd = nodeCreate(path)
    nodeTruncate(fd, size.toDouble())
    try {
        for ((offset, bytes) in writes) {
            if (bytes.isNotEmpty()) nodeWrite(fd, bytes.toLatin1(), offset.toDouble())
        }
    } finally {
        nodeCloseFixture(fd)
    }
    return path
}

internal actual fun deleteFixture(path: String) {
    nodeDelete(path)
}

private fun ByteArray.toLatin1(): String = buildString(size) {
    for (byte in this@toLatin1) append((byte.toInt() and 0xFF).toChar())
}

private fun nodeTempPath(): String = js(
    """(() => {
        const os = globalThis.process.getBuiltinModule('node:os');
        const path = globalThis.process.getBuiltinModule('node:path');
        const crypto = globalThis.process.getBuiltinModule('node:crypto');
        return path.join(os.tmpdir(), 'tassel-sparse-' + crypto.randomBytes(8).toString('hex'));
        })()""",
)

private fun nodeCreate(path: String): Int = js("globalThis.process.getBuiltinModule('node:fs').openSync(path, 'w')")

private fun nodeTruncate(fd: Int, size: Double) {
    js("globalThis.process.getBuiltinModule('node:fs').ftruncateSync(fd, size)")
}

private fun nodeWrite(fd: Int, text: String, position: Double) {
    js(
        """(() => {
        const fs = globalThis.process.getBuiltinModule('node:fs');
        const Buffer = globalThis.process.getBuiltinModule('node:buffer').Buffer;
        const buf = Buffer.from(text, 'latin1');
        fs.writeSync(fd, buf, 0, buf.length, position);
        })()""",
    )
}

private fun nodeCloseFixture(fd: Int) {
    js("globalThis.process.getBuiltinModule('node:fs').closeSync(fd)")
}

private fun nodeDelete(path: String) {
    js("globalThis.process.getBuiltinModule('node:fs').rmSync(path, { force: true })")
}
