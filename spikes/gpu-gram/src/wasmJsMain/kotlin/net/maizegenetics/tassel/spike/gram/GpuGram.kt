@file:OptIn(ExperimentalWasmJsInterop::class)

package net.maizegenetics.tassel.spike.gram

import kotlin.js.ExperimentalWasmJsInterop
import kotlin.js.JsAny
import kotlin.time.TimeSource

internal class OpenGpu(
    val device: JsAny,
    val naive: JsAny,
    val tiled: JsAny,
    val adapter: AdapterReportJson,
    val setupSeconds: Double,
)

internal fun splitShader(source: String): Pair<String, String> {
    val naiveAt = source.indexOf("// @@NAIVE")
    val tiledAt = source.indexOf("// @@TILED")
    require(naiveAt > 0 && tiledAt > naiveAt) { "gram_f32.wgsl is missing entry markers" }
    val common = source.substring(0, naiveAt)
    return common + source.substring(naiveAt, tiledAt) to common + source.substring(tiledAt)
}

internal fun dispatchSpans(itemCount: Int, itemsPerGroup: Int): List<Pair<Int, Int>> {
    val spans = ArrayList<Pair<Int, Int>>()
    var base = 0
    while (base < itemCount) {
        val groups = minOf(65_535, (itemCount - base + itemsPerGroup - 1) / itemsPerGroup)
        spans += base to groups
        base += groups * itemsPerGroup
    }
    return spans
}

internal suspend fun openGpu(shaderSource: String): OpenGpu {
    val started = TimeSource.Monotonic.markNow()
    val adapter = awaitJs(requestAdapter()) ?: error("requestAdapter returned null")
    val info = awaitJs(requestAdapterInfo(adapter))
    val device = awaitJs(requestDevice(adapter)) ?: error("requestDevice returned null")
    val (naiveSource, tiledSource) = splitShader(shaderSource)
    val naiveModule = createShaderModule(device, naiveSource)
    val tiledModule = createShaderModule(device, tiledSource)
    val naive = awaitJs(createPipeline(device, naiveModule, "gram_naive")) ?: error("naive pipeline was null")
    val tiled = awaitJs(createPipeline(device, tiledModule, "gram_tiled")) ?: error("tiled pipeline was null")
    val requested = LIMIT_NAMES.associateWith { limitOf(adapter, it).toLong() }
    val granted = LIMIT_NAMES.associateWith { limitOf(device, it).toLong() }
    val report = if (info == null) {
        AdapterReportJson(
            vendor = "",
            architecture = "",
            device = "",
            description = "",
            isFallbackAdapter = false,
            requestedLimits = requested,
            grantedLimits = granted,
        )
    } else {
        AdapterReportJson(
            vendor = infoField(info, "vendor"),
            architecture = infoField(info, "architecture"),
            device = infoField(info, "device"),
            description = infoField(info, "description"),
            isFallbackAdapter = infoFlag(info, "isFallbackAdapter"),
            requestedLimits = requested,
            grantedLimits = granted,
        )
    }
    return OpenGpu(device, naive, tiled, report, secondsSince(started))
}

internal class PhaseClock {
    val staging = ArrayList<Double>()
    val writeBuffer = ArrayList<Double>()
    val dispatch = ArrayList<Double>()
    val mapAsync = ArrayList<Double>()
    val copyBack = ArrayList<Double>()
    val cpuAccumulation = ArrayList<Double>()
    val endToEnd = ArrayList<Double>()

    fun add(
        stagingSeconds: Double,
        writeSeconds: Double,
        dispatchSeconds: Double,
        mapSeconds: Double,
        copySeconds: Double,
        cpuSeconds: Double,
    ) {
        staging += stagingSeconds
        writeBuffer += writeSeconds
        dispatch += dispatchSeconds
        mapAsync += mapSeconds
        copyBack += copySeconds
        cpuAccumulation += cpuSeconds
        endToEnd += stagingSeconds + writeSeconds + dispatchSeconds + mapSeconds + copySeconds + cpuSeconds
    }

    fun toReport(): PhaseReport = PhaseReport(
        staging = statsOf(staging).toJson(),
        writeBuffer = statsOf(writeBuffer).toJson(),
        dispatch = statsOf(dispatch).toJson(),
        mapAsync = statsOf(mapAsync).toJson(),
        copyBack = statsOf(copyBack).toJson(),
        cpuAccumulation = statsOf(cpuAccumulation).toJson(),
        endToEnd = statsOf(endToEnd).toJson(),
    )
}

internal class GramRun(val phases: PhaseReport, val triangle: DoubleArray)

/**
 * One configuration of the Gram. [nextPanel] is restarted for every timed run.
 * [accumulateOnCpu] reads each panel back and adds it in f64. Otherwise the shader adds f32
 * partials and the triangle is read once at the end.
 */
internal suspend fun runConfiguredGram(
    gpu: OpenGpu,
    pipeline: JsAny,
    taxa: Int,
    shader: String,
    reset: () -> Unit,
    nextPanel: () -> FloatPanel?,
    linearMemory: Boolean,
    accumulateOnCpu: Boolean,
    runs: Int,
    warmup: Boolean = true,
): GramRun {
    val pairs = pairCount(taxa)
    val triangleBytes = pairs * 4
    val panelBytes = taxa * PANEL_COLUMNS * 4
    val spans = if (shader == "tiled") {
        dispatchSpans(upperTiles(taxa), itemsPerGroup = 1)
    } else {
        dispatchSpans(pairs, itemsPerGroup = 256)
    }
    val device = gpu.device
    val panelBuffer = createBuffer(device, panelBytes, GPU_STORAGE or GPU_COPY_DST)
    val triangleBuffer = createBuffer(device, triangleBytes, GPU_STORAGE or GPU_COPY_SRC or GPU_COPY_DST)
    val readBuffer = createBuffer(device, triangleBytes, GPU_MAP_READ or GPU_COPY_DST)
    val uniform = createBuffer(device, 16, GPU_UNIFORM or GPU_COPY_DST)
    val bindGroup = createBindGroup(device, pipeline, panelBuffer, triangleBuffer, uniform)
    val queue = queueOf(device)
    val clock = PhaseClock()
    lateinit var last: DoubleArray
    try {
        val iterations = if (warmup) runs + 1 else runs
        repeat(iterations) { iteration ->
            reset()
            val triangle = DoubleArray(pairs)
            var stagingSeconds = 0.0
            var writeSeconds = 0.0
            var dispatchSeconds = 0.0
            var mapSeconds = 0.0
            var copySeconds = 0.0
            var cpuSeconds = 0.0
            var wroteAny = false
            while (true) {
                val stageMark = TimeSource.Monotonic.markNow()
                val panel = nextPanel() ?: break
                val staged = if (linearMemory) {
                    var stagedAt = stageMark
                    writePanelLinear(queue, panelBuffer, panel.values) { stagedAt = TimeSource.Monotonic.markNow() }
                    stagingSeconds += secondsSince(stageMark) - secondsBetween(stagedAt, TimeSource.Monotonic.markNow())
                    writeSeconds += secondsBetween(stagedAt, TimeSource.Monotonic.markNow())
                    true
                } else {
                    val array = createFloat32(panel.values.size)
                    for (index in panel.values.indices) setFloat32(array, index, panel.values[index])
                    stagingSeconds += secondsSince(stageMark)
                    val writeMark = TimeSource.Monotonic.markNow()
                    writeFloats(queue, panelBuffer, array)
                    writeSeconds += secondsSince(writeMark)
                    true
                }
                if (!staged) continue
                for ((offset, groups) in spans) {
                    val clear = if (accumulateOnCpu) {
                        offset == spans.first().first
                    } else {
                        !wroteAny &&
                            offset == spans.first().first
                    }
                    val dispatchMark = TimeSource.Monotonic.markNow()
                    writeParams(queue, uniform, taxa, panel.columns, panel.columns, offset)
                    val done = submitDispatch(
                        device,
                        pipeline,
                        bindGroup,
                        groups,
                        triangleBuffer,
                        clear,
                        readBuffer,
                        doRead = false,
                        readBytes = triangleBytes,
                    )
                    awaitJs(done)
                    dispatchSeconds += secondsSince(dispatchMark)
                    wroteAny = true
                }
                if (accumulateOnCpu) {
                    val read = readback(device, triangleBuffer, readBuffer, triangleBytes)
                    mapSeconds += read.mapSeconds
                    val copyMark = TimeSource.Monotonic.markNow()
                    val mapped = mappedFloats(readBuffer)
                    copySeconds += secondsSince(copyMark)
                    val cpuMark = TimeSource.Monotonic.markNow()
                    readFloats(mapped, triangle, add = true)
                    cpuSeconds += secondsSince(cpuMark)
                }
            }
            if (!accumulateOnCpu && wroteAny) {
                val read = readback(device, triangleBuffer, readBuffer, triangleBytes)
                mapSeconds += read.mapSeconds
                val copyMark = TimeSource.Monotonic.markNow()
                val mapped = mappedFloats(readBuffer)
                readFloats(mapped, triangle, add = false)
                copySeconds += secondsSince(copyMark)
            }
            if (warmup && iteration == 0) return@repeat
            clock.add(stagingSeconds, writeSeconds, dispatchSeconds, mapSeconds, copySeconds, cpuSeconds)
            last = triangle
        }
    } finally {
        destroyJs(panelBuffer)
        destroyJs(triangleBuffer)
        destroyJs(readBuffer)
        destroyJs(uniform)
    }
    return GramRun(clock.toReport(), last)
}

private class ReadTiming(val mapSeconds: Double)

private fun copyToReadback(device: JsAny, triangle: JsAny, readBuffer: JsAny, readBytes: Int): JsAny = js(
    """
    (() => {
        const encoder = device.createCommandEncoder();
        encoder.copyBufferToBuffer(triangle, 0, readBuffer, 0, readBytes);
        device.queue.submit([encoder.finish()]);
        return device.queue.onSubmittedWorkDone();
    })()
    """,
)

private suspend fun readback(device: JsAny, triangle: JsAny, readBuffer: JsAny, readBytes: Int): ReadTiming {
    val mark = TimeSource.Monotonic.markNow()
    awaitJs(copyToReadback(device, triangle, readBuffer, readBytes))
    awaitJs(mapRead(readBuffer))
    return ReadTiming(secondsSince(mark))
}

private fun secondsBetween(start: TimeSource.Monotonic.ValueTimeMark, end: TimeSource.Monotonic.ValueTimeMark): Double =
    (end - start).inWholeNanoseconds / 1e9

private fun upperTiles(taxa: Int): Int {
    val tiles = (taxa + 15) / 16
    return tiles * (tiles + 1) / 2
}
