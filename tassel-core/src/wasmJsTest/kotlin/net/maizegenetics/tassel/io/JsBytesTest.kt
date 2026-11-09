@file:OptIn(ExperimentalWasmJsInterop::class)

package net.maizegenetics.tassel.io

import kotlin.js.ExperimentalWasmJsInterop
import kotlin.js.JsAny
import kotlin.test.Test
import kotlin.test.assertEquals

/**
 * The byte copy the browser file picker uses. Node has no `FileReader`, so the test builds the `ArrayBuffer`
 * the reader would have produced. Bytes 0 through 255 are the ones a UTF-8 text read would replace.
 */
class JsBytesTest {
    @Test
    fun bytesZeroThrough255SurviveAnArrayBuffer() {
        val bytes = arrayBufferToByteArray(filledArrayBuffer(256))

        assertEquals(256, bytes.size)
        for (index in 0..255) assertEquals(index.toByte(), bytes[index])
    }

    @Test
    fun aBufferLongerThanOneChunkKeepsEveryByte() {
        val bytes = arrayBufferToByteArray(filledArrayBuffer(CHUNK_BOUNDARY))

        assertEquals(CHUNK_BOUNDARY, bytes.size)
        assertEquals(0x80.toByte(), bytes[0x80])
        assertEquals(0x9F.toByte(), bytes[0x9F])
        assertEquals((CHUNK - 1).toByte(), bytes[CHUNK - 1])
        assertEquals(0.toByte(), bytes[CHUNK])
        assertEquals((CHUNK_BOUNDARY - 1 and 0xFF).toByte(), bytes[CHUNK_BOUNDARY - 1])
    }

    @Test
    fun latin1CodeUnitsZeroThrough255Survive() {
        val text = CharArray(256) { it.toChar() }.concatToString()

        val bytes = latin1ToByteArray(text)

        for (index in 0..255) assertEquals(index.toByte(), bytes[index])
    }

    private companion object {
        const val CHUNK = 32_768
        const val CHUNK_BOUNDARY = CHUNK + 256
    }
}

/** An `ArrayBuffer` of [length] bytes whose value at index `i` is `i` modulo 256. */
private fun filledArrayBuffer(length: Int): JsAny = js(
    """(() => {
        const buffer = new ArrayBuffer(length);
        const bytes = new Uint8Array(buffer);
        for (let i = 0; i < length; i++) bytes[i] = i & 255;
        return buffer;
    })()""",
)
