package net.maizegenetics.tassel.app.grid

import net.maizegenetics.tassel.genotype.GenotypeTable
import kotlin.math.ceil
import kotlin.math.floor
import kotlin.math.min

/**
 * The block of cells a scrolled grid shows, as table indices: taxa [firstTaxon] through `firstTaxon + taxonCount - 1`
 * and sites [firstSite] through `firstSite + siteCount - 1`.
 *
 * A window is the unit of work for one frame. It covers only what is on screen, so its size depends on the viewport
 * and the cell size, never on the size of the table.
 */
data class GridWindow(val firstTaxon: Int, val taxonCount: Int, val firstSite: Int, val siteCount: Int) {
    init {
        require(firstTaxon >= 0 && firstSite >= 0) { "window starts before the table: $this" }
        require(taxonCount >= 0 && siteCount >= 0) { "window has a negative size: $this" }
    }

    /** Table taxon indices in the window. */
    val taxa: IntRange get() = firstTaxon until firstTaxon + taxonCount

    /** Table site indices in the window. */
    val sites: IntRange get() = firstSite until firstSite + siteCount

    /** Number of cells in the window. */
    val cellCount: Int get() = taxonCount * siteCount

    companion object {
        /** The window of a grid with nothing to show. */
        val Empty: GridWindow = GridWindow(0, 0, 0, 0)
    }
}

/**
 * The window a grid of [numberOfTaxa] by [numberOfSites] shows when it is scrolled [taxonOffset] pixels down and
 * [siteOffset] pixels right inside a viewport of [width] by [height] pixels, with cells [cellWidth] by [cellHeight].
 *
 * Cells that are only partly visible are included at both edges, so the window is one wider and one taller than the
 * viewport covers exactly. Offsets past the end of the table clamp to the last cell rather than throwing, because a
 * fling can overshoot.
 */
fun gridWindow(
    numberOfTaxa: Int,
    numberOfSites: Int,
    taxonOffset: Float,
    siteOffset: Float,
    width: Float,
    height: Float,
    cellWidth: Float,
    cellHeight: Float,
): GridWindow {
    require(cellWidth > 0f && cellHeight > 0f) { "cells must have a positive size: $cellWidth x $cellHeight" }
    if (numberOfTaxa <= 0 || numberOfSites <= 0 || width <= 0f || height <= 0f) return GridWindow.Empty
    val firstTaxon = firstVisible(taxonOffset, cellHeight, numberOfTaxa)
    val firstSite = firstVisible(siteOffset, cellWidth, numberOfSites)
    return GridWindow(
        firstTaxon = firstTaxon,
        taxonCount = visibleCount(taxonOffset, height, cellHeight, numberOfTaxa, firstTaxon),
        firstSite = firstSite,
        siteCount = visibleCount(siteOffset, width, cellWidth, numberOfSites, firstSite),
    )
}

/** Index of the first cell the viewport touches, clamped into the table. */
private fun firstVisible(offset: Float, cellSize: Float, count: Int): Int =
    floor(offset / cellSize).toInt().coerceIn(0, count - 1)

/** Number of cells from [first] that the viewport touches, clamped to the end of the table. */
private fun visibleCount(offset: Float, viewportSize: Float, cellSize: Float, count: Int, first: Int): Int {
    val end = ceil((offset + viewportSize) / cellSize).toInt()
    return (min(end, count) - first).coerceAtLeast(1)
}

/**
 * Genotype calls of one [window], site-major: the call for window row `row` at window column `column` is at index
 * `column * window.taxonCount + row`, the same layout `SiteBlock` uses.
 */
class WindowGenotypes(val window: GridWindow, private val genotypes: ByteArray) {
    init {
        require(genotypes.size == window.cellCount) {
            "expected ${window.cellCount} genotypes for $window, got ${genotypes.size}"
        }
    }

    /** The call at window row [row], window column [column]. */
    operator fun get(row: Int, column: Int): Byte {
        if (row < 0 || row >= window.taxonCount) throw IndexOutOfBoundsException("row $row of ${window.taxonCount}")
        if (column < 0 || column >= window.siteCount) {
            throw IndexOutOfBoundsException("column $column of ${window.siteCount}")
        }
        return genotypes[column * window.taxonCount + row]
    }

    /** The call at table taxon [taxon], table site [site], both of which must be inside [window]. */
    fun atTable(taxon: Int, site: Int): Byte = get(taxon - window.firstTaxon, site - window.firstSite)
}

/**
 * Reads exactly the cells of [window] out of this table, with random access.
 *
 * Random access is the right tool here, even though analyses use `siteBlocks`. `siteBlocks` always starts at site 0
 * and walks forward, and its default implementation allocates and fills a `ByteArray` per block, so reaching a scroll
 * offset with `siteBlocks(...).drop(n)` would do n allocations and n block fills on every frame, because
 * `Sequence.drop` still evaluates what it discards. `GenotypeTable`'s KDoc sanctions random access "for code that
 * needs it", and a virtualized viewport needs it: this reads `window.cellCount` genotypes, a few thousand for any
 * viewport, whatever the size of the table and wherever the grid is scrolled to.
 */
fun GenotypeTable.readWindow(window: GridWindow): WindowGenotypes {
    if (window.cellCount == 0) return WindowGenotypes(window, ByteArray(0))
    require(window.taxa.last < numberOfTaxa && window.sites.last < numberOfSites) {
        "$window is outside a $numberOfTaxa x $numberOfSites table"
    }
    val genotypes = ByteArray(window.cellCount)
    var at = 0
    for (site in window.sites) {
        for (taxon in window.taxa) genotypes[at++] = genotype(taxon, site)
    }
    return WindowGenotypes(window, genotypes)
}

/**
 * A run of consecutive window rows in one column that share a call, which is what the grid paints as a single
 * rectangle. Genotype data repeats down a column, so merging runs cuts the number of draw calls well below the number
 * of cells without changing a pixel.
 */
data class ColumnRun(val column: Int, val firstRow: Int, val rowCount: Int, val genotype: Byte)

/** The runs of [column], top to bottom. */
fun WindowGenotypes.columnRuns(column: Int): List<ColumnRun> {
    if (window.taxonCount == 0) return emptyList()
    val runs = ArrayList<ColumnRun>()
    var start = 0
    var current = get(0, column)
    for (row in 1 until window.taxonCount) {
        val genotype = get(row, column)
        if (genotype != current) {
            runs.add(ColumnRun(column, start, row - start, current))
            start = row
            current = genotype
        }
    }
    runs.add(ColumnRun(column, start, window.taxonCount - start, current))
    return runs
}
