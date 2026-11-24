package net.maizegenetics.tassel.spike.gram

import kotlinx.coroutines.runBlocking
import net.maizegenetics.tassel.distance.DistanceMatrix
import net.maizegenetics.tassel.distance.centeredIbs
import java.nio.file.Files
import java.nio.file.Path
import kotlin.math.abs

fun main(args: Array<String>) {
    val output = Path.of(args.firstOrNull() ?: "jvm-cpu.json")
    val report = cpuReport("jvm", osDescription())
    val text = cpuReportJson(report)
    Files.createDirectories(output.parent ?: Path.of("."))
    Files.writeString(output, text)
    println("wrote $output")
    for (size in report.sizes) {
        println(
            "taxa=${size.taxa} sites=${size.sites} checksum=${size.checksum} " +
                "median=${size.dosage.medianSeconds} centeredIbsMaxAbs=${size.centeredIbsMaxAbs}",
        )
    }
}

internal fun osDescription(): String = listOf("os.name", "os.arch", "os.version")
    .joinToString(" ") { System.getProperty(it) }

internal fun cpuReport(target: String, os: String): CpuReportJson {
    val sizes = timeDosageShapes(runs = 3, log = ::println).map { size ->
        if (size.taxa == 281 && size.sites == 3_093) {
            val delta = centeredIbsDelta(generateDosages(size.taxa, size.sites), size.taxa, size.sites)
            size.copy(centeredIbsMaxAbs = delta.first, centeredIbsSumPkAbs = delta.second)
        } else {
            size
        }
    }
    return CpuReportJson(target = target, os = os, sizes = sizes)
}

private fun centeredIbsDelta(codes: ByteArray, taxa: Int, sites: Int): Pair<Double, Double> = runBlocking {
    val ours = centeredGram(codes, taxa, sites)
    val matrix = centeredIbs(genotypeTable(codes, taxa, sites))
    var maxAbs = 0.0
    var index = 0
    for (first in 0 until taxa) {
        for (second in first until taxa) {
            val scaled = (ours.values[index] / ours.sumPk).toFloat()
            val diff = abs(scaled - matrix.distance(first, second)).toDouble()
            if (diff > maxAbs) maxAbs = diff
            index++
        }
    }
    val theirs = matrix.annotations.getValue(DistanceMatrix.CENTERED_IBS_SUMPK).toDouble()
    maxAbs to abs(ours.sumPk - theirs)
}
