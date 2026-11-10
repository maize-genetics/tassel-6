package net.maizegenetics.tassel.spike.gram

class FloatPanel(val values: FloatArray, val columns: Int)

/** Tight f32 dosage panels of [PANEL_COLUMNS], restartable for repeated timings. */
class DosageCursor(private val codes: ByteArray, private val taxa: Int, private val sites: Int) {
    private var site: Int = 0

    fun reset() {
        site = 0
    }

    fun next(): FloatPanel? {
        if (site >= sites) return null
        val columns = minOf(PANEL_COLUMNS, sites - site)
        val panel = dosageFloatPanel(codes, taxa, site, columns)
        site += columns
        return FloatPanel(panel, columns)
    }
}

/** Centered columns in the same order as [centeredGram], packed tight for upload. */
class CenteredCursor(private val genotypes: ByteArray, private val taxa: Int, private val sites: Int) {
    private var walker = CenteredWalker(genotypes, taxa, sites)

    fun reset() {
        walker = CenteredWalker(genotypes, taxa, sites)
    }

    fun next(): FloatPanel? {
        val panel = nextCenteredFloatPanel(walker, taxa, PANEL_COLUMNS) ?: return null
        return FloatPanel(panel, panel.size / taxa)
    }
}
