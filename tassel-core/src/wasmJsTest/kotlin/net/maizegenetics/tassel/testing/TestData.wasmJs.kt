@file:OptIn(ExperimentalWasmJsInterop::class)

package net.maizegenetics.tassel.testing

// Tests run on Node as an ES module, where `require` is not defined; `process.getBuiltinModule`
// (Node 22.3+) is the way to reach `fs`.

private fun nodeEnvironmentVariable(name: String): String? = js("globalThis.process.env[name] ?? null")

// latin1 maps each byte to one UTF-16 unit, so file bytes cross the JS boundary as a string.
private fun nodeReadFileLatin1(path: String): String =
    js("globalThis.process.getBuiltinModule('node:fs').readFileSync(path, 'latin1')")

internal actual fun environmentVariable(name: String): String? = nodeEnvironmentVariable(name)

internal actual fun readFile(absolutePath: String): ByteArray {
    val text = nodeReadFileLatin1(absolutePath)
    return ByteArray(text.length) { text[it].code.toByte() }
}
