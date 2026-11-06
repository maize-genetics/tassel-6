package net.maizegenetics.tassel.distance

import net.maizegenetics.tassel.genotype.AlleleStats
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.table.doubleToString

/** TASSEL 5's `KinshipPlugin` default for the number of allele states per site to evaluate. */
const val DEFAULT_MAX_ALLELES: Int = 6

/** Columns of centered dosage kept for one product. A few megabytes, independent of the host. */
private const val PANEL_BYTES: Long = 4L * 1024L * 1024L

private const val MIN_PANEL_COLUMNS: Long = 8L

private const val MAX_PANEL_COLUMNS: Long = 4096L

/** Floor on sites between progress reports. The tutorial file reports about every 100 sites. */
internal const val KINSHIP_REPORT_FLOOR: Int = 100

/**
 * Columns of storage reserved for one dosage panel at [numberOfTaxa]. The kernel flushes earlier
 * when a progress boundary arrives, so this is the ceiling, not the width of every product.
 */
fun centeredIbsPanelColumns(numberOfTaxa: Int): Int {
    require(numberOfTaxa > 0) { "centeredIbsPanelColumns: numberOfTaxa must be positive: $numberOfTaxa" }
    val fit = PANEL_BYTES / (numberOfTaxa.toLong() * Double.SIZE_BYTES)
    return fit.coerceIn(MIN_PANEL_COLUMNS, MAX_PANEL_COLUMNS).toInt()
}

/**
 * Sites read between progress reports and panel flushes. A function of the table shape only, so
 * the `Double` sums do not depend on the machine (ADR-0020).
 */
fun centeredIbsPanelSites(numberOfTaxa: Int, numberOfSites: Int): Int {
    require(numberOfTaxa > 0) { "centeredIbsPanelSites: numberOfTaxa must be positive: $numberOfTaxa" }
    require(numberOfSites >= 0) { "centeredIbsPanelSites: numberOfSites must not be negative: $numberOfSites" }
    if (numberOfSites == 0) return 0
    val cap = maxOf(KINSHIP_REPORT_FLOOR, numberOfSites / 40)
    return minOf(centeredIbsPanelColumns(numberOfTaxa), cap).coerceAtLeast(1)
}

/**
 * Centered IBS (Endelman) kinship for every pair of taxa in [table].
 *
 * Each evaluated allele becomes a column of centered dosages, `copies - 2p` for a called genotype
 * and `0` when the genotype is not comparable, and the triangle accumulates `Z Z^T`
 * ([KinshipTriangleUpdate]). Missing calls therefore contribute nothing to either taxon of a pair.
 * The most frequent `min(alleles - 1, maxAlleles - 1)` alleles are evaluated, monomorphic sites
 * contribute no column, and `SumPk` gains `2 p (1 - p)` per column. The matrix is scaled by
 * `SumPk` and stored as [Float], which is the text format's precision (ADR-0020).
 *
 * [maxAlleles] is how many allele states per site to evaluate, 2 to 6 inclusive.
 *
 * @see net.maizegenetics.analysis.distance.EndelmanDistanceMatrix
 */
suspend fun centeredIbs(
    table: GenotypeTable,
    maxAlleles: Int = DEFAULT_MAX_ALLELES,
    progress: ProgressReporter = ProgressReporter.None,
): DistanceMatrix {
    require(maxAlleles in 2..DEFAULT_MAX_ALLELES) {
        "centeredIbs: maxAlleles must be between 2 and $DEFAULT_MAX_ALLELES inclusive: $maxAlleles"
    }
    val numberOfTaxa = table.numberOfTaxa
    val numberOfSites = table.numberOfSites
    val triangle = DoubleArray(DistanceMatrix.triangleSize(numberOfTaxa))
    if (numberOfTaxa == 0 || numberOfSites == 0) {
        return finishCenteredIbs(table, triangle, 0.0)
    }
    val stride = centeredIbsPanelColumns(numberOfTaxa)
    val interval = centeredIbsPanelSites(numberOfTaxa, numberOfSites)
    val panel = DoubleArray(numberOfTaxa * stride)
    val genotypes = ByteArray(numberOfTaxa)
    var filled = 0
    var sumPi = 0.0
    var sitesDone = 0
    var sitesInPanel = 0

    suspend fun flush(reportProgress: Boolean) {
        if (filled > 0) {
            kinshipTriangleUpdate.add(panel, numberOfTaxa, filled, stride, triangle)
            filled = 0
        }
        if (reportProgress) progress.report(sitesDone, numberOfSites)
    }

    for (block in table.siteBlocks()) {
        for (site in block.sites) {
            block.copySiteGenotypes(site, genotypes)
            val stats = AlleleStats.of(genotypes)
            val evaluated = minOf(stats.alleles.size - 1, maxAlleles - 1)
            if (evaluated > 0) {
                if (filled + evaluated > stride) flush(reportProgress = false)
                val gametes = stats.alleleGameteCount.toFloat()
                for (alleleIndex in 0 until evaluated) {
                    val frequency = stats.alleleCounts[alleleIndex].toFloat() / gametes
                    sumPi += frequency * (1.0 - frequency)
                    val twiceP = frequency.toDouble() * 2.0
                    val allele = stats.alleles[alleleIndex].toInt()
                    val column = filled
                    for (taxon in 0 until numberOfTaxa) {
                        panel[taxon * stride + column] =
                            centeredDosage(genotypes[taxon].toInt(), allele, twiceP)
                    }
                    filled++
                }
            }
            sitesDone++
            sitesInPanel++
            if (sitesInPanel == interval || sitesDone == numberOfSites) {
                flush(reportProgress = true)
                sitesInPanel = 0
            }
        }
    }
    return finishCenteredIbs(table, triangle, sumPi * 2.0)
}

/** Scales an upper triangle of unnormalized pair sums into TASSEL 5's `Float` matrix. */
internal fun finishCenteredIbs(table: GenotypeTable, unscaled: FloatArray, sumPk: Double): DistanceMatrix {
    val numberOfTaxa = table.numberOfTaxa
    val scaled = FloatArray(unscaled.size)
    var index = 0
    for (first in 0 until numberOfTaxa) {
        for (second in first until numberOfTaxa) {
            scaled[second * (second + 1) / 2 + first] = (unscaled[index] / sumPk).toFloat()
            index++
        }
    }
    return centeredIbsMatrix(table, scaled, sumPk)
}

/** Scales an upper triangle of unnormalized pair sums into TASSEL 5's `Float` matrix. */
internal fun finishCenteredIbs(table: GenotypeTable, unscaled: DoubleArray, sumPk: Double): DistanceMatrix {
    val numberOfTaxa = table.numberOfTaxa
    val scaled = FloatArray(unscaled.size)
    var index = 0
    for (first in 0 until numberOfTaxa) {
        for (second in first until numberOfTaxa) {
            scaled[second * (second + 1) / 2 + first] = (unscaled[index] / sumPk).toFloat()
            index++
        }
    }
    return centeredIbsMatrix(table, scaled, sumPk)
}

private fun centeredIbsMatrix(table: GenotypeTable, scaled: FloatArray, sumPk: Double): DistanceMatrix = DistanceMatrix(
    taxa = table.taxa,
    lowerTriangle = scaled,
    annotations = mapOf(
        DistanceMatrix.MATRIX_TYPE to KinshipMethod.Centered_IBS.name,
        DistanceMatrix.CENTERED_IBS_SUMPK to doubleToString(sumPk),
    ),
)
