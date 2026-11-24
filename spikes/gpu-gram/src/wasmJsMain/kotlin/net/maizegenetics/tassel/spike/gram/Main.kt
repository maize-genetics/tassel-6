package net.maizegenetics.tassel.spike.gram

import kotlinx.coroutines.DelicateCoroutinesApi
import kotlinx.coroutines.GlobalScope
import kotlinx.coroutines.launch

@OptIn(DelicateCoroutinesApi::class)
fun main() {
    if (isBrowser()) {
        logLine("wasm main")
        GlobalScope.launch {
            try {
                runBrowser()
            } catch (error: Throwable) {
                logLine("spike failed: ${error.message}")
            }
        }
    } else {
        val path = nodeEnv("TASSEL_GPU_GRAM_OUT") ?: "build/reports/gpu-gram/node-cpu.json"
        val report = CpuReportJson(
            target = "wasm-node",
            os = nodeOs(),
            sizes = timeDosageShapes(runs = 3),
        )
        nodeWrite(path, cpuReportJson(report))
        println("wrote $path")
    }
}
