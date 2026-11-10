package net.maizegenetics.tassel.spike.gram

import kotlin.math.abs
import kotlin.math.max

/**
 * Adds one panel's `Z Zᵀ` into a packed upper triangle.
 *
 * [panel] is row-major, one row per taxon, with [stride] stored columns and [columns] of them live.
 * Each pair uses one f64 accumulator across those columns, then one add into [triangle]. A full
 * matrix is the sum of successive panels of [PANEL_COLUMNS] columns, in column order. That panel
 * width is part of the centered-value reference, because f64 addition is not associative. Raw
 * `0/1/2` dosages are exact either way: every pair sum in the measured shapes is at most `80000`,
 * well below `2^24` and `2^53`.
 *
 * The triangle is packed like [triangleIndex]: `(0,0), (0,1), ... (0,n-1), (1,1), ...`.
 */
fun addPanel(panel: DoubleArray, taxa: Int, columns: Int, stride: Int, triangle: DoubleArray) {
    if (taxa == 0 || columns == 0) return
    var index = 0
    for (first in 0 until taxa) {
        val row = first * stride
        var second = first
        // Four pairs at a time. Each pair still adds columns in order into its own accumulator,
        // so the sums match a single-accumulator walk.
        val limit = taxa - 3
        while (second < limit) {
            val other0 = second * stride
            val other1 = other0 + stride
            val other2 = other1 + stride
            val other3 = other2 + stride
            var sum0 = 0.0
            var sum1 = 0.0
            var sum2 = 0.0
            var sum3 = 0.0
            var column = 0
            while (column < columns) {
                val value = panel[row + column]
                sum0 += value * panel[other0 + column]
                sum1 += value * panel[other1 + column]
                sum2 += value * panel[other2 + column]
                sum3 += value * panel[other3 + column]
                column++
            }
            triangle[index] += sum0
            triangle[index + 1] += sum1
            triangle[index + 2] += sum2
            triangle[index + 3] += sum3
            index += 4
            second += 4
        }
        while (second < taxa) {
            val other = second * stride
            var sum = 0.0
            var column = 0
            while (column < columns) {
                sum += panel[row + column] * panel[other + column]
                column++
            }
            triangle[index] += sum
            index++
            second++
        }
    }
}

/** Raw dosage Gram. Missing codes contribute `0`. */
fun dosageGram(codes: ByteArray, taxa: Int, sites: Int, panelColumns: Int = PANEL_COLUMNS): DoubleArray {
    val triangle = DoubleArray(pairCount(taxa))
    val panel = DoubleArray(taxa * panelColumns)
    var site = 0
    while (site < sites) {
        val columns = minOf(panelColumns, sites - site)
        fillDosagePanel(codes, taxa, site, columns, panelColumns, panel)
        addPanel(panel, taxa, columns, panelColumns, triangle)
        site += columns
    }
    return triangle
}

fun fillDosagePanel(codes: ByteArray, taxa: Int, siteStart: Int, columns: Int, stride: Int, panel: DoubleArray) {
    for (taxon in 0 until taxa) {
        val row = taxon * stride
        var column = 0
        while (column < columns) {
            val code = codes[(siteStart + column) * taxa + taxon].toInt()
            panel[row + column] = if (code == MISSING_CODE) 0.0 else code.toDouble()
            column++
        }
    }
}

/** Packs a tight f32 dosage panel (`stride == columns`) for upload. Missing contributes `0`. */
fun dosageFloatPanel(codes: ByteArray, taxa: Int, siteStart: Int, columns: Int): FloatArray {
    val panel = FloatArray(taxa * columns)
    for (taxon in 0 until taxa) {
        val row = taxon * columns
        var column = 0
        while (column < columns) {
            val code = codes[(siteStart + column) * taxa + taxon].toInt()
            panel[row + column] = if (code == MISSING_CODE) 0f else code.toFloat()
            column++
        }
    }
    return panel
}

class ValueError(val maxAbs: Double, val maxRel: Double, val bitExact: Boolean)

fun compareTriangles(actual: DoubleArray, reference: DoubleArray): ValueError {
    require(actual.size == reference.size) { "triangles differ: ${actual.size} vs ${reference.size}" }
    var maxAbs = 0.0
    var maxRel = 0.0
    var bitExact = true
    for (index in reference.indices) {
        val diff = abs(actual[index] - reference[index])
        if (diff > maxAbs) maxAbs = diff
        val relative = diff / max(abs(reference[index]), RELATIVE_FLOOR)
        if (relative > maxRel) maxRel = relative
        if (actual[index] != reference[index]) bitExact = false
    }
    return ValueError(maxAbs, maxRel, bitExact)
}

/**
 * Max absolute difference of kinship values after dividing by [sumPk] and rounding to [Float],
 * which is the ADR-0020 comparison (the bound is [SCALED_TOLERANCE], not the unscaled Gram).
 */
fun scaledMaxAbs(actual: DoubleArray, reference: DoubleArray, sumPk: Double): Double {
    require(actual.size == reference.size) { "triangles differ: ${actual.size} vs ${reference.size}" }
    var maxAbs = 0.0
    for (index in reference.indices) {
        val diff = abs((actual[index] / sumPk).toFloat() - (reference[index] / sumPk).toFloat()).toDouble()
        if (diff > maxAbs) maxAbs = diff
    }
    return maxAbs
}
