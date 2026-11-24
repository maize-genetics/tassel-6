@file:OptIn(ExperimentalWasmJsInterop::class)

package net.maizegenetics.tassel.spike.gram

import kotlin.js.ExperimentalWasmJsInterop
import kotlin.js.JsAny

internal suspend fun runBrowser() {
    val host = machineHost()
    val quick = pageQuery("quick") == "1"
    val runs = pageQuery("runs")?.toIntOrNull()?.coerceIn(1, 5) ?: if (quick) 1 else 3
    val shapes = if (quick) listOf(GRAM_SHAPES.first()) else GRAM_SHAPES
    logLine("GPU Gram spike  quick=$quick  runs=$runs  host=$host")
    val fallback = demonstrateFallbacks()
    logLine("forced failure: ${fallback.forcedFailure}")
    logLine("lost device: ${fallback.lostDevice}")
    val opened = openMeasuredGpu()
    var gpu = opened.gpu
    var status = opened.status
    if (gpu != null) {
        val selfTest = selfTest(gpu)
        if (selfTest != null) {
            status = "self-test failed: $selfTest"
            logLine(status)
            destroyJs(gpu.device)
            gpu = null
        } else {
            logLine("self-test passed")
        }
    }
    val sizes = ArrayList<BrowserSizeJson>()
    for (shape in shapes) {
        logLine("shape ${shape.label}")
        sizes += measureShape(shape, runs, gpu, warmup = !quick)
        postReport(browserReport(host, status, gpu, fallback, sizes, quick = quick, complete = false))
    }
    val finalReport = browserReport(host, status, gpu, fallback, sizes, quick = quick, complete = true)
    postReport(finalReport)
    logLine("break-even: ${finalReport.breakEven}")
    logLine("done")
    gpu?.let { destroyJs(it.device) }
}

private suspend fun machineHost(): String = try {
    awaitJs(fetchText("/host"))?.let { jsText(it) }?.trim().orEmpty()
} catch (error: Throwable) {
    logLine("host lookup failed: ${error.message}")
    ""
}

private suspend fun postReport(report: BrowserReportJson) {
    val json = browserReportJson(report)
    try {
        val response = awaitJs(postJson(json))
        logLine("posted ${json.length} bytes ${response?.let { jsText(it) } ?: ""}")
    } catch (error: Throwable) {
        logLine("report post failed: ${error.message}")
    }
}

private fun browserReport(
    host: String,
    status: String,
    gpu: OpenGpu?,
    fallback: FallbackDemoJson,
    sizes: List<BrowserSizeJson>,
    quick: Boolean,
    complete: Boolean,
): BrowserReportJson = BrowserReportJson(
    host = host,
    userAgent = userAgent(),
    quick = quick,
    complete = complete,
    setupSeconds = gpu?.setupSeconds,
    adapter = gpu?.adapter,
    gpuStatus = status,
    fallback = fallback,
    sizes = sizes,
    breakEven = breakEvenLabel(sizes),
)

private fun adapterLabel(adapter: AdapterReportJson): String =
    listOf(adapter.vendor, adapter.architecture, adapter.device, adapter.description)
        .map { it.trim() }
        .filter { it.isNotEmpty() }
        .distinct()
        .joinToString(" ")
        .ifEmpty { "unnamed adapter" }

private class OpenedGpu(val gpu: OpenGpu?, val status: String)

private suspend fun openMeasuredGpu(): OpenedGpu {
    if (!hasNavigatorGpu()) return OpenedGpu(null, "no navigator.gpu")
    return try {
        val shader = awaitJs(fetchText("/gram_f32.wgsl"))?.let { jsText(it) } ?: return OpenedGpu(null, "empty shader")
        val gpu = openGpu(shader)
        val label = adapterLabel(gpu.adapter)
        val status = if (gpu.adapter.isFallbackAdapter) {
            "fallback adapter $label; CPU only"
        } else {
            "ok $label"
        }
        logLine("adapter $status  setup ${gpu.setupSeconds}s")
        if (gpu.adapter.isFallbackAdapter) {
            destroyJs(gpu.device)
            OpenedGpu(null, status)
        } else {
            OpenedGpu(gpu, status)
        }
    } catch (error: Throwable) {
        OpenedGpu(null, "GPU setup failed: ${error.message}")
    }
}

private suspend fun selfTest(gpu: OpenGpu): String? {
    val codes = byteArrayOf(2, 0, 1, 2)
    val expected = dosageGram(codes, 2, 2, panelColumns = 2)
    val cursor = DosageCursor(codes, 2, 2)
    val shaders = listOf("naive" to gpu.naive, "tiled" to gpu.tiled)
    for ((name, pipeline) in shaders) {
        for (linear in listOf(true, false)) {
            val run = runConfiguredGram(
                gpu = gpu,
                pipeline = pipeline,
                taxa = 2,
                shader = name,
                reset = cursor::reset,
                nextPanel = cursor::next,
                linearMemory = linear,
                accumulateOnCpu = false,
                runs = 1,
                warmup = false,
            )
            val error = compareTriangles(run.triangle, expected)
            if (!error.bitExact) {
                return "$name linear=$linear maxAbs=${error.maxAbs}"
            }
        }
    }
    return null
}

private suspend fun measureShape(shape: GramShape, runs: Int, gpu: OpenGpu?, warmup: Boolean): BrowserSizeJson {
    val codes = generateDosages(shape.taxa, shape.sites)
    val checksum = dosageChecksum(codes).toString()
    logLine("${shape.label} checksum $checksum")
    val (cpu, dosageReference) = timeRunsReturning(runs) {
        dosageGram(codes, shape.taxa, shape.sites)
    }
    logLine("${shape.label} wasm cpu median ${cpu.medianSeconds}s")
    val centered = centeredGram(codes, shape.taxa, shape.sites)
    val results = ArrayList<GpuResultJson>()
    if (gpu != null) {
        val genotypes = encodeGenotypes(codes)
        val dosage = DosageCursor(codes, shape.taxa, shape.sites)
        val centeredCursor = CenteredCursor(genotypes, shape.taxa, shape.sites)
        val shaders = listOf("naive" to gpu.naive, "tiled" to gpu.tiled)
        for ((shaderName, pipeline) in shaders) {
            for ((linear, transfer) in listOf(true to "linear-memory", false to "per-element")) {
                for ((onCpu, accumulation) in listOf(false to "f32", true to "f64-cpu")) {
                    results += gpuResult(
                        gpu, pipeline, shaderName, transfer, linear, accumulation, onCpu,
                        shape, runs, warmup, dosage::reset, dosage::next, dosageReference, centered, panel = "dosage",
                    )
                    results += gpuResult(
                        gpu, pipeline, shaderName, transfer, linear, accumulation, onCpu,
                        shape, runs, warmup, centeredCursor::reset, centeredCursor::next, centered.values, centered,
                        panel = "centered",
                    )
                }
            }
        }
    }
    return BrowserSizeJson(shape.taxa, shape.sites, checksum, cpu.toJson(), results)
}

private suspend fun gpuResult(
    gpu: OpenGpu,
    pipeline: JsAny,
    shader: String,
    transfer: String,
    linear: Boolean,
    accumulation: String,
    onCpu: Boolean,
    shape: GramShape,
    runs: Int,
    warmup: Boolean,
    reset: () -> Unit,
    nextPanel: () -> FloatPanel?,
    reference: DoubleArray,
    centered: CenteredGram,
    panel: String,
): GpuResultJson {
    val label = "${shape.label} $shader $transfer $accumulation $panel"
    logLine(label)
    return try {
        val run = runConfiguredGram(
            gpu, pipeline, shape.taxa, shader, reset, nextPanel, linear, onCpu, runs, warmup,
        )
        val compared = compareTriangles(run.triangle, reference)
        if (panel == "dosage") {
            logLine("$label bitExact=${compared.bitExact} maxAbs=${compared.maxAbs}")
            GpuResultJson(
                shader = shader,
                transfer = transfer,
                accumulation = accumulation,
                panel = panel,
                phases = run.phases,
                dosageBitExact = compared.bitExact,
            )
        } else {
            val scaled = scaledMaxAbs(run.triangle, reference, centered.sumPk)
            logLine("$label maxAbs=${compared.maxAbs} scaled=$scaled")
            GpuResultJson(
                shader = shader,
                transfer = transfer,
                accumulation = accumulation,
                panel = panel,
                phases = run.phases,
                centeredMaxAbs = compared.maxAbs,
                centeredMaxRel = compared.maxRel,
                scaledMaxAbs = scaled,
                within8e6 = scaled <= SCALED_TOLERANCE,
            )
        }
    } catch (error: Throwable) {
        logLine("$label failed: ${error.message}")
        GpuResultJson(
            shader = shader,
            transfer = transfer,
            accumulation = accumulation,
            panel = panel,
            error = (error.message ?: error.toString()).take(500),
        )
    }
}

private suspend fun demonstrateFallbacks(): FallbackDemoJson {
    val cpu = dosageGram(byteArrayOf(2, 0, 1, 2), 2, 2, panelColumns = 2).toList()
    if (!hasNavigatorGpu()) return FallbackDemoJson("no navigator.gpu", "no navigator.gpu", cpu)
    val forced = try {
        forcedShaderFailure()
    } catch (error: Throwable) {
        "cpu fallback after shader failure: ${error.message}"
    }
    val lost = try {
        lostDeviceDemo()
    } catch (error: Throwable) {
        "cpu fallback after lost device: ${error.message}"
    }
    return FallbackDemoJson(forced, lost, cpu)
}

private suspend fun freshDevice(): JsAny? {
    // Chrome consumes an adapter on the first requestDevice, so each demo asks for its own.
    val adapter = awaitJs(requestAdapter()) ?: return null
    return awaitJs(requestDevice(adapter))
}

private suspend fun forcedShaderFailure(): String {
    val device = freshDevice() ?: return "cpu fallback: no device for forced failure"
    try {
        val module = createShaderModule(device, "@compute @workgroup_size(1) fn bad() { let x: f32 = \"no\"; }")
        val info = awaitJs(compilationInfo(module))
        val errors = if (info == null) 0 else shaderErrorCount(info)
        if (errors > 0) return "cpu fallback after shader failure ($errors errors)"
        return try {
            awaitJs(createPipeline(device, module, "bad"))
            "forced shader failure did not fail"
        } catch (error: Throwable) {
            "cpu fallback after shader failure: ${error.message}"
        }
    } finally {
        destroyJs(device)
    }
}

private suspend fun lostDeviceDemo(): String {
    val device = freshDevice() ?: return "cpu fallback: no device to lose"
    val lost = deviceLost(device)
    destroyJs(device)
    val info = awaitJs(lost)
    val message = if (info == null) "destroyed" else lostMessage(info).ifBlank { "destroyed" }
    return "cpu fallback after lost device: $message"
}
