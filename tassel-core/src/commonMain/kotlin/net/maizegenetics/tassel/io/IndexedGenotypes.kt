package net.maizegenetics.tassel.io

import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.ops.ProgressReporter

/**
 * A bgzipped VCF or BCF opened with its index.
 *
 * [open] reads the first inflated bytes and chooses the reader. BCF 2.2 magic selects [IndexedBcf]. Anything else,
 * including a VCF, stays on [IndexedVcf.open], so a TBI is still a VCF index. This type does not search a filesystem.
 * [IndexedVcf.open] and [Bcf.read] stay the entry points they are (ADR-0022).
 *
 * [close] closes the BGZF reader only. The caller owns both sources.
 */
class IndexedGenotypes private constructor(private val reader: Reader) : AutoCloseable {
    /** Contigs the opened reader reports, in that reader's order. */
    val contigs: List<Contig> get() = reader.contigs

    /**
     * One contig named by the opened reader.
     *
     * For a BCF, [name] is a header-dictionary contig and [length] is present only when that `##contig` line has a
     * length. For a VCF, this is [IndexedVcf.Contig].
     */
    data class Contig(val name: String, val length: Long?)

    /**
     * The sites in [region].
     *
     * An unknown contig and an inverted region are errors. A known contig with nothing in range, or a region past
     * the index's coordinate limit, is a table with the file's taxa and no sites.
     */
    suspend fun read(region: VcfRegion, progress: ProgressReporter = ProgressReporter.None): IndexedGenotypeRead =
        reader.read(region, progress)

    override fun close() = reader.close()

    private interface Reader : AutoCloseable {
        val contigs: List<Contig>

        suspend fun read(region: VcfRegion, progress: ProgressReporter): IndexedGenotypeRead
    }

    private class BcfReader(private val indexed: IndexedBcf) : Reader {
        override val contigs: List<Contig> get() = indexed.contigs

        override suspend fun read(region: VcfRegion, progress: ProgressReporter): IndexedGenotypeRead =
            indexed.read(region, progress)

        override fun close() = indexed.close()
    }

    private class VcfReader(private val indexed: IndexedVcf) : Reader {
        override val contigs: List<Contig> = indexed.contigs.map { Contig(it.name, it.length) }

        override suspend fun read(region: VcfRegion, progress: ProgressReporter): IndexedGenotypeRead {
            val read = indexed.read(region, progress)
            return IndexedGenotypeRead(read.table, read.warnings)
        }

        override fun close() = indexed.close()
    }

    companion object {
        /**
         * Opens [data] with [index].
         *
         * Neither source is closed. The magic peek uses its own reader and leaves [data] positioned for the reader
         * that then takes over.
         */
        suspend fun open(data: RandomAccessSource, index: RandomAccessSource): IndexedGenotypes {
            val reader = if (isBcf22(data)) {
                BcfReader(IndexedBcf.open(data, index))
            } else {
                VcfReader(IndexedVcf.open(data, index))
            }
            return IndexedGenotypes(reader)
        }
    }
}

/** A region read and the warnings found while inflating it. */
class IndexedGenotypeRead(val table: GenotypeTable, val warnings: List<String>)

/**
 * True when the inflated bytes at offset 0 are the BCF 2.2 magic.
 *
 * The peek reader is closed here. [exactBlocks] keeps the check from inflating past the first block.
 */
private suspend fun isBcf22(data: RandomAccessSource): Boolean {
    val reader = BgzfReader(data)
    try {
        reader.exactBlocks(true)
        val magic = ByteArray(BCF_MAGIC_LENGTH)
        var at = 0
        while (at < magic.size) {
            val count = reader.read(magic, at, magic.size - at)
            if (count <= 0) break
            at += count
        }
        reader.seek(0)
        return isBcf22Magic(magic)
    } finally {
        reader.close()
    }
}
