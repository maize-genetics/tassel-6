package net.maizegenetics.tassel.spike.gram

import kotlinx.serialization.Serializable
import kotlinx.serialization.encodeToString
import kotlinx.serialization.json.Json

val ReportJson: Json = Json {
    prettyPrint = true
    explicitNulls = false
    encodeDefaults = true
}

@Serializable
data class GeneratorInfo(
    val name: String = "park-miller",
    val multiplier: Int = PARK_MILLER_MULTIPLIER,
    val modulus: Int = PARK_MILLER_MODULUS,
    val seed: Int = PARK_MILLER_SEED,
    val missingPercent: Int = MISSING_PERCENT,
    val mafMinPercent: Int = MAF_MIN_PERCENT,
    val mafMaxPercent: Int = MAF_MAX_PERCENT,
    val panelColumns: Int = PANEL_COLUMNS,
    val checksum: String = "u64 wrapping, start 0, checksum = checksum * 131 + code, site-major",
)

@Serializable
data class RunStatsJson(val minSeconds: Double, val medianSeconds: Double, val samplesSeconds: List<Double>)

fun RunStats.toJson(): RunStatsJson = RunStatsJson(minSeconds, medianSeconds, samplesSeconds)

@Serializable
data class PhaseReport(
    val staging: RunStatsJson,
    val writeBuffer: RunStatsJson,
    val dispatch: RunStatsJson,
    val mapAsync: RunStatsJson,
    val copyBack: RunStatsJson,
    val cpuAccumulation: RunStatsJson,
    val endToEnd: RunStatsJson,
)

@Serializable
data class GpuResultJson(
    val shader: String,
    val transfer: String,
    val accumulation: String,
    val panel: String,
    val phases: PhaseReport? = null,
    val dosageBitExact: Boolean? = null,
    val centeredMaxAbs: Double? = null,
    val centeredMaxRel: Double? = null,
    val scaledMaxAbs: Double? = null,
    val within8e6: Boolean? = null,
    val error: String? = null,
)

@Serializable
data class CpuSizeJson(
    val taxa: Int,
    val sites: Int,
    val checksum: String,
    val preview: List<Int>,
    val dosage: RunStatsJson,
    val centeredIbsMaxAbs: Double? = null,
    val centeredIbsSumPkAbs: Double? = null,
)

@Serializable
data class CpuReportJson(
    val target: String,
    val os: String,
    val generator: GeneratorInfo = GeneratorInfo(),
    val sizes: List<CpuSizeJson>,
)

@Serializable
data class AdapterReportJson(
    val vendor: String,
    val architecture: String,
    val device: String,
    val description: String,
    val isFallbackAdapter: Boolean,
    val requestedLimits: Map<String, Long>,
    val grantedLimits: Map<String, Long>,
)

@Serializable
data class FallbackDemoJson(val forcedFailure: String, val lostDevice: String, val cpuSelfTest: List<Double>)

@Serializable
data class BrowserSizeJson(
    val taxa: Int,
    val sites: Int,
    val checksum: String,
    val cpuDosage: RunStatsJson,
    val gpu: List<GpuResultJson>,
    val note: String? = null,
)

@Serializable
data class BrowserReportJson(
    val host: String = "",
    val userAgent: String,
    val quick: Boolean,
    val complete: Boolean,
    val setupSeconds: Double? = null,
    val adapter: AdapterReportJson? = null,
    val gpuStatus: String,
    val fallback: FallbackDemoJson,
    val generator: GeneratorInfo = GeneratorInfo(),
    val sizes: List<BrowserSizeJson>,
    val breakEven: String,
)

fun cpuReportJson(report: CpuReportJson): String = ReportJson.encodeToString(report)

fun browserReportJson(report: BrowserReportJson): String = ReportJson.encodeToString(report)

/**
 * Smallest measured shapes where the linear-memory f32 dosage Gram beats the wasm CPU median.
 * `"none"` is a valid result.
 */
fun breakEvenLabel(sizes: List<BrowserSizeJson>): String {
    val hits = sizes.mapNotNull { size ->
        val cpu = size.cpuDosage.medianSeconds
        val gpu = size.gpu
            .filter {
                it.panel == "dosage" &&
                    it.transfer == "linear-memory" &&
                    it.accumulation == "f32" &&
                    it.error == null &&
                    it.phases != null
            }
            .minOfOrNull { it.phases!!.endToEnd.medianSeconds }
        if (gpu != null && gpu < cpu) size.taxa.toString() + "x" + size.sites else null
    }
    return if (hits.isEmpty()) "none" else hits.joinToString(", ")
}
