@file:OptIn(ExperimentalWasmJsInterop::class)

package net.maizegenetics.tassel.io

import kotlin.js.ExperimentalWasmJsInterop
import kotlin.js.JsAny

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
