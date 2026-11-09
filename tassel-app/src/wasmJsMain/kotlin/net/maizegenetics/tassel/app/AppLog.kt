package net.maizegenetics.tassel.app

import kotlin.js.ExperimentalWasmJsInterop

@OptIn(ExperimentalWasmJsInterop::class)
internal actual fun logAppMessage(line: String) {
    logToConsole(line)
}

@OptIn(ExperimentalWasmJsInterop::class)
private fun logToConsole(line: String) {
    js("console.log(line)")
}
