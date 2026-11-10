@file:OptIn(ExperimentalWasmJsInterop::class)

package net.maizegenetics.tassel.spike.gram

import kotlin.coroutines.resume
import kotlin.coroutines.resumeWithException
import kotlin.coroutines.suspendCoroutine
import kotlin.js.ExperimentalWasmJsInterop
import kotlin.js.JsAny

internal fun isBrowser(): Boolean = js("typeof window !== 'undefined'")

internal fun hasNavigatorGpu(): Boolean = js("typeof navigator !== 'undefined' && !!navigator.gpu")

internal fun userAgent(): String = js("String(navigator.userAgent)")

internal fun pageQuery(name: String): String? = js("new URLSearchParams(window.location.search).get(name)")

internal fun logLine(line: String): Unit = js(
    """
    (() => {
        const el = document.getElementById('log');
        if (el) el.textContent += line + '\n';
        console.log(line);
    })()
    """,
)

internal fun nodeEnv(name: String): String? = js("globalThis.process?.env?.[name] ?? null")

internal fun nodeOs(): String = js("globalThis.process.platform + ' ' + globalThis.process.arch")

internal fun nodeWrite(path: String, text: String): Unit = js(
    """
    (() => {
        const fs = globalThis.process.getBuiltinModule('node:fs');
        const pathMod = globalThis.process.getBuiltinModule('node:path');
        fs.mkdirSync(pathMod.dirname(path), { recursive: true });
        fs.writeFileSync(path, text);
    })()
    """,
)

private fun jsThen(promise: JsAny, onOk: (JsAny?) -> JsAny?, onErr: (JsAny?) -> JsAny?): Unit =
    js("promise.then(onOk, onErr)")

private fun stringify(value: JsAny?): String = js("value == null ? 'null' : String(value)")

internal suspend fun awaitJs(promise: JsAny): JsAny? = suspendCoroutine { continuation ->
    jsThen(
        promise,
        onOk = { value ->
            continuation.resume(value)
            null
        },
        onErr = { error ->
            continuation.resumeWithException(IllegalStateException(stringify(error)))
            null
        },
    )
}

internal fun requestAdapter(): JsAny = js("navigator.gpu.requestAdapter({ powerPreference: 'high-performance' })")

internal fun requestAdapterInfo(adapter: JsAny): JsAny = js(
    """
    (() => {
        if (adapter.info) return Promise.resolve(adapter.info);
        if (adapter.requestAdapterInfo) return adapter.requestAdapterInfo();
        return Promise.resolve(null);
    })()
    """,
)

internal fun requestDevice(adapter: JsAny): JsAny = js(
    """
    (() => {
        const limits = adapter.limits;
        const wanted = {};
        const names = [
            'maxBufferSize',
            'maxStorageBufferBindingSize',
            'maxComputeWorkgroupSizeX',
            'maxComputeInvocationsPerWorkgroup',
            'maxComputeWorkgroupsPerDimension',
            'maxComputeWorkgroupStorageSize'
        ];
        for (const name of names) wanted[name] = limits[name];
        return adapter.requestDevice({ requiredLimits: wanted });
    })()
    """,
)

internal fun limitOf(source: JsAny, name: String): Double = js("Number(source.limits[name])")

internal fun infoField(info: JsAny, field: String): String = js("String(info?.[field] ?? '')")

internal fun infoFlag(info: JsAny, field: String): Boolean = js("!!info?.[field]")

internal fun createShaderModule(device: JsAny, code: String): JsAny = js("device.createShaderModule({ code: code })")

internal fun compilationInfo(module: JsAny): JsAny = js("module.getCompilationInfo()")

internal fun shaderErrorCount(info: JsAny): Int =
    js("info.messages.filter(function (message) { return message.type === 'error'; }).length")

internal fun createPipeline(device: JsAny, module: JsAny, entryPoint: String): JsAny = js(
    """
    device.createComputePipelineAsync({
        layout: 'auto',
        compute: { module: module, entryPoint: entryPoint }
    })
    """,
)

internal fun createBuffer(device: JsAny, size: Int, usage: Int): JsAny =
    js("device.createBuffer({ size: size, usage: usage })")

internal fun createBindGroup(device: JsAny, pipeline: JsAny, panel: JsAny, triangle: JsAny, uniform: JsAny): JsAny = js(
    """
    device.createBindGroup({
        layout: pipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: panel } },
            { binding: 1, resource: { buffer: triangle } },
            { binding: 2, resource: { buffer: uniform } }
        ]
    })
    """,
)

internal fun writeParams(queue: JsAny, uniform: JsAny, taxa: Int, columns: Int, stride: Int, pairBase: Int): Unit =
    js("queue.writeBuffer(uniform, 0, new Uint32Array([taxa, columns, stride, pairBase]))")

internal fun queueOf(device: JsAny): JsAny = js("device.queue")

internal fun submitDispatch(
    device: JsAny,
    pipeline: JsAny,
    bindGroup: JsAny,
    groups: Int,
    triangle: JsAny,
    clear: Boolean,
    readback: JsAny,
    doRead: Boolean,
    readBytes: Int,
): JsAny = js(
    """
    (() => {
        const encoder = device.createCommandEncoder();
        if (clear) encoder.clearBuffer(triangle);
        const pass = encoder.beginComputePass();
        pass.setPipeline(pipeline);
        pass.setBindGroup(0, bindGroup);
        pass.dispatchWorkgroups(groups);
        pass.end();
        if (doRead) encoder.copyBufferToBuffer(triangle, 0, readback, 0, readBytes);
        device.queue.submit([encoder.finish()]);
        return device.queue.onSubmittedWorkDone();
    })()
    """,
)

internal fun mapRead(buffer: JsAny): JsAny = js("buffer.mapAsync(1)")

internal fun mappedFloats(buffer: JsAny): JsAny = js(
    """
    (() => {
        const range = new Float32Array(buffer.getMappedRange());
        const copy = new Float32Array(range.length);
        copy.set(range);
        buffer.unmap();
        return copy;
    })()
    """,
)

internal fun destroyJs(value: JsAny): Unit = js("value.destroy()")

internal fun deviceLost(device: JsAny): JsAny = js("device.lost")

internal fun lostMessage(info: JsAny): String = js("String(info?.message ?? '')")

internal fun fetchText(url: String): JsAny = js(
    """
    fetch(url).then(function (response) {
        if (!response.ok) throw new Error('HTTP ' + response.status + ' for ' + url);
        return response.text();
    })
    """,
)

internal fun postJson(json: String): JsAny = js(
    """
    fetch('/report', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: json
    }).then(function (response) {
        if (!response.ok) throw new Error('report HTTP ' + response.status);
        return response.text();
    })
    """,
)

internal fun jsText(value: JsAny): String = js("String(value)")

internal fun createFloat32(length: Int): JsAny = js("new Float32Array(length)")

internal fun setFloat32(array: JsAny, index: Int, value: Float): Unit = js("array[index] = value")

internal fun writeFloats(queue: JsAny, buffer: JsAny, source: JsAny): Unit = js("queue.writeBuffer(buffer, 0, source)")

internal const val GPU_MAP_READ: Int = 0x0001

internal const val GPU_COPY_SRC: Int = 0x0004

internal const val GPU_COPY_DST: Int = 0x0008

internal const val GPU_UNIFORM: Int = 0x0040

internal const val GPU_STORAGE: Int = 0x0080

internal val LIMIT_NAMES: List<String> = listOf(
    "maxBufferSize",
    "maxStorageBufferBindingSize",
    "maxComputeWorkgroupSizeX",
    "maxComputeInvocationsPerWorkgroup",
    "maxComputeWorkgroupsPerDimension",
    "maxComputeWorkgroupStorageSize",
)
