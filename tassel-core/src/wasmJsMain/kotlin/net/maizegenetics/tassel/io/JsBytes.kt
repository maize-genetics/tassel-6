@file:OptIn(ExperimentalWasmJsInterop::class)

package net.maizegenetics.tassel.io

import kotlin.js.ExperimentalWasmJsInterop
import kotlin.js.JsAny
import kotlin.js.toJsString

/**
 * Copies a JavaScript string whose code units are bytes.
 *
 * latin1, as Node's `Buffer` encoding uses it, maps each byte to one UTF-16 unit, so values 0 through 255
 * survive. The string comes from one JavaScript call; this loop does not call back into JavaScript.
 * The Node file source and [arrayBufferToByteArray] both use it.
 */
fun latin1ToByteArray(latin1: String): ByteArray = ByteArray(latin1.length) { latin1[it].code.toByte() }

/**
 * Copies an `ArrayBuffer` into a [ByteArray].
 *
 * One JavaScript call turns the buffer into a latin1 string with chunked `String.fromCharCode`, which keeps
 * every byte. `TextDecoder('latin1')` is windows-1252 and would remap `0x80` through `0x9F`. No `node:fs`:
 * this source set is linked into the browser (ADR-0021).
 */
fun arrayBufferToByteArray(buffer: JsAny): ByteArray = latin1ToByteArray(arrayBufferLatin1(buffer))

/**
 * Copies [bytes] into an `ArrayBuffer`.
 *
 * Each byte becomes one UTF-16 code unit, the reverse of [arrayBufferToByteArray]. Chunks stay at
 * [LATIN1_CHUNK] so `String.fromCharCode` on the way back, and this copy on the way out, never pass a
 * longer argument list than the read path does. `TextEncoder` would turn bytes `0x80` through `0xFF`
 * into multi-byte UTF-8. The browser export builds its Blob from this buffer.
 */
fun byteArrayToArrayBuffer(bytes: ByteArray): JsAny {
    val buffer = allocateArrayBuffer(bytes.size)
    var offset = 0
    while (offset < bytes.size) {
        val end = minOf(offset + LATIN1_CHUNK, bytes.size)
        val chunk = CharArray(end - offset) { (bytes[offset + it].toInt() and 0xFF).toChar() }.concatToString()
        copyLatin1(buffer, offset, chunk.toJsString())
        offset = end
    }
    return buffer
}

private const val LATIN1_CHUNK: Int = 32_768

private fun allocateArrayBuffer(size: Int): JsAny = js("new ArrayBuffer(size)")

private fun copyLatin1(buffer: JsAny, offset: Int, latin1: JsAny) {
    js(
        """(() => {
            const bytes = new Uint8Array(buffer);
            for (let i = 0; i < latin1.length; i++) bytes[offset + i] = latin1.charCodeAt(i) & 255;
        })()""",
    )
}

private fun arrayBufferLatin1(buffer: JsAny): String = js(
    """(() => {
        const bytes = new Uint8Array(buffer);
        const parts = [];
        const chunk = 32768;
        for (let i = 0; i < bytes.length; i += chunk) {
            const end = Math.min(i + chunk, bytes.length);
            parts.push(String.fromCharCode.apply(null, bytes.subarray(i, end)));
        }
        return parts.join('');
    })()""",
)
