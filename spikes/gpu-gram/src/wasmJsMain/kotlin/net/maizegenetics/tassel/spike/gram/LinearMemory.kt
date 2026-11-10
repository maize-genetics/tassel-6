@file:OptIn(ExperimentalWasmJsInterop::class, UnsafeWasmMemoryApi::class)

package net.maizegenetics.tassel.spike.gram

import kotlin.js.ExperimentalWasmJsInterop
import kotlin.js.JsAny
import kotlin.wasm.unsafe.UnsafeWasmMemoryApi
import kotlin.wasm.unsafe.wasmMemory
import kotlin.wasm.unsafe.withScopedMemoryAllocator

/**
 * Copies [panel] into Wasm linear memory and issues one `writeBuffer` from that view.
 *
 * Kotlin arrays live on the WasmGC heap, so this is an explicit staging copy. [allocate] is
 * 8-byte aligned, which satisfies `Float32Array`. The view is only used inside the allocator
 * scope; `writeBuffer` copies it before the scope returns.
 */
internal fun writePanelLinear(queue: JsAny, buffer: JsAny, panel: FloatArray, onStaged: () -> Unit) {
    withScopedMemoryAllocator { allocator ->
        val pointer = allocator.allocate(panel.size * 4)
        var cursor = pointer
        for (value in panel) {
            cursor.storeInt(value.toRawBits())
            cursor += 4
        }
        onStaged()
        writeLinearFloats(queue, buffer, wasmMemory, pointer.address.toInt(), panel.size)
    }
}

internal fun readFloats(source: JsAny, dest: DoubleArray, add: Boolean) {
    withScopedMemoryAllocator { allocator ->
        val pointer = allocator.allocate(dest.size * 4)
        copyFloatsToMemory(source, wasmMemory, pointer.address.toInt(), dest.size)
        var cursor = pointer
        for (index in dest.indices) {
            val value = Float.fromBits(cursor.loadInt()).toDouble()
            if (add) dest[index] += value else dest[index] = value
            cursor += 4
        }
    }
}

private fun writeLinearFloats(queue: JsAny, buffer: JsAny, memory: JsAny, byteOffset: Int, length: Int): Unit =
    js("queue.writeBuffer(buffer, 0, new Float32Array(memory.buffer, byteOffset, length))")

private fun copyFloatsToMemory(source: JsAny, memory: JsAny, byteOffset: Int, length: Int): Unit =
    js("new Float32Array(memory.buffer, byteOffset, length).set(source)")
