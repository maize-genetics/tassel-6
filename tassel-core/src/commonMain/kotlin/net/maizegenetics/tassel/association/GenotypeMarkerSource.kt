package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.analysis.ReferenceProbability
import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.AlleleStats
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.genotype.SiteBlock
import net.maizegenetics.tassel.join.GenotypePhenotype
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.position.Position

/**
 * The [MarkerSource] over the sites of a [GenotypePhenotype]. Genotypes are read through
 * `siteBlocks(taxa = genotypeTaxonIndices())`, one block at a time and once per `scan` however many observation sets
 * share it; nothing is allocated per observation per marker beyond the reusable [MarkerColumns].
 *
 * Genotype classes and what counts as missing follow TASSEL 5's GLM, which works on IUPAC genotype strings:
 *
 * - A class is an IUPAC string, so phased heterozygotes (`A/C`, `C/A`) are one class, in filters and in the model
 *   (TASSEL 5 splits them in its filters), and every indel heterozygote is the class `0`.
 * - An observation is missing at a marker when it has no genotype row or either allele is unknown (`N/N`, `A/N`).
 *   TASSEL 5 reads `A/N` as the string `A`, a homozygote; that is a recorded departure.
 * - For `additive` and `additive-dominance`, an observation with a rare or undefined allele (`Z`, `X`) is also missing,
 *   since it has no copy count. For `genotype` it is an ordinary class.
 *
 * Instances are not thread safe: a scan reuses buffers. A block's loaders from [MarkerBlock.fork] have buffers of their
 * own and may load at the same time.
 *
 * @see net.maizegenetics.analysis.association.DiscreteSitesFELM
 * @see net.maizegenetics.analysis.association.ReferenceProbabilityFELM
 * @see net.maizegenetics.analysis.association.AbstractFixedEffectLM
 */
internal class GenotypeMarkerSource private constructor(
    private val data: GenotypePhenotype,
    private val reference: ReferenceProbability?,
) : MarkerSource {
    override val numberOfMarkers: Int get() = data.genotypes.numberOfSites

    override fun position(marker: Int): Position = data.genotypes.positions[marker]

    override fun scan(observationSets: List<IntArray>, options: MarkerOptions, blockSize: Int): Sequence<MarkerBlock> {
        require(blockSize > 0) { "blockSize must be positive: $blockSize" }
        if (options.encoding == MarkerEncoding.dosage) {
            requireNotNull(reference) { "The dosage encoding needs a source made with GenotypeMarkerSource.withDosage" }
        }
        val genotyped = data.genotypeTaxonIndices()
        val blockPosition = IntArray(data.genotypes.numberOfTaxa) { -1 }
        genotyped.forEachIndexed { position, row -> blockPosition[row] = position }
        val sets = observationSets.map { observations ->
            val taxon = IntArray(observations.size) {
                val observation = observations[it]
                require(observation in 0 until data.numberOfObservations) {
                    "Observation $observation is outside 0 until ${data.numberOfObservations}"
                }
                data.genotypeRow(observation).let { row -> if (row < 0) -1 else blockPosition[row] }
            }
            ObservationSet(taxon)
        }
        val encoder = Encoder(options, sets, genotyped, reference)
        return data.genotypes.siteBlocks(blockSize, genotyped).map { encoder.block(it) }
    }

    /** For each observation of a set, its taxon's position in the block, or `-1` when it has no genotype. */
    private class ObservationSet(val blockTaxon: IntArray) {
        val size: Int get() = blockTaxon.size
    }

    private class Encoder(
        val options: MarkerOptions,
        val sets: List<ObservationSet>,
        val genotyped: IntArray,
        val reference: ReferenceProbability?,
    ) {
        private val maxSet = sets.maxOfOrNull { it.size } ?: 0
        private val siteGenotypes = ByteArray(genotyped.size)
        private val dosages = DoubleArray(genotyped.size)
        private val classOf = IntArray(maxSet)
        private val classCounts = IntArray(GenotypeClasses.size)
        private val newIndex = IntArray(GenotypeClasses.size)
        private val firstGenotype = IntArray(GenotypeClasses.size)
        private val alleleCounts = IntArray(Allele.NUCLEOTIDE_STATES)
        private val alleleOrder = IntArray(Allele.NUCLEOTIDE_STATES)

        private var loadedSite = -1
        private var stats: AlleleStats? = null
        private var dosagesLoaded = false

        fun block(block: SiteBlock): MarkerBlock {
            loadedSite = -1
            return object : MarkerBlock {
                override val markers: IntRange get() = block.sites

                // A site block is a read-only view of genotypes, so encoders with their own buffers can share it.
                override fun fork(): MarkerBlock = Encoder(options, sets, genotyped, reference).block(block)

                override fun load(marker: Int, set: Int, into: MarkerColumns): Boolean {
                    require(marker in block.sites) { "Marker $marker is not in block ${block.sites}" }
                    require(set in sets.indices) { "Observation set $set of ${sets.size}" }
                    require(into.capacity >= sets[set].size) {
                        "MarkerColumns holds ${into.capacity} observations, set $set has ${sets[set].size}"
                    }
                    into.reset()
                    if (marker != loadedSite) {
                        block.copySiteGenotypes(marker, siteGenotypes)
                        loadedSite = marker
                        stats = null
                        dosagesLoaded = false
                    }
                    return if (options.encoding == MarkerEncoding.dosage) {
                        encodeDosage(block, marker, sets[set], into)
                    } else {
                        encodeCalls(marker, sets[set], into)
                    }
                }
            }
        }

        private fun encodeDosage(block: SiteBlock, marker: Int, set: ObservationSet, into: MarkerColumns): Boolean {
            if (!dosagesLoaded) {
                reference!!.column(block, marker, genotyped, dosages)
                dosagesLoaded = true
            }
            // NaN is missing. TASSEL 5's NaN check for reference probabilities never fires (a recorded departure).
            var kept = 0
            for (i in 0 until set.size) {
                val position = set.blockTaxon[i]
                if (position >= 0 && !dosages[position].isNaN()) kept++
            }
            if (kept == 0) return fail(into, SkipReason.NoObservations)
            into.columns(listOf(DOSAGE_LABEL), byteArrayOf(Allele.UNKNOWN), Allele.UNKNOWN)
            val column = into.column(0)
            var k = 0
            for (i in 0 until set.size) {
                val position = set.blockTaxon[i]
                if (position < 0 || dosages[position].isNaN()) continue
                into.observation[k] = i
                column[k] = dosages[position]
                k++
            }
            into.kept(k, minorObs = 0)
            return true
        }

        private fun encodeCalls(marker: Int, set: ObservationSet, into: MarkerColumns): Boolean {
            val alleleBased = options.encoding != MarkerEncoding.genotype
            classCounts.fill(0)
            var present = 0
            for (i in 0 until set.size) {
                val position = set.blockTaxon[i]
                var c = -1
                if (position >= 0) {
                    val genotype = siteGenotypes[position].toInt() and 0xFF
                    if (!alleleBased || GenotypeClasses.scorable[genotype]) c = GenotypeClasses.classOf[genotype]
                }
                classOf[i] = c
                if (c >= 0) {
                    classCounts[c]++
                    present++
                }
            }
            if (present == 0) return fail(into, SkipReason.NoObservations)

            var distinct = 0
            var heterozygous = 0
            var largest = 0
            var second = 0
            for (c in classCounts.indices) {
                val count = classCounts[c]
                if (count == 0) continue
                distinct++
                if (GenotypeClasses.heterozygous[c]) heterozygous++
                if (count > largest) {
                    second = largest
                    largest = count
                } else if (count > second) {
                    second = count
                }
            }
            val minorObs = if (distinct > 1) second else 0

            if (options.biallelicOnly && !(distinct == 2 || (distinct == 3 && heterozygous == 1))) {
                return fail(into, SkipReason.NotBiallelic, minorObs)
            }

            if (options.minMac > 0 || options.minMaf > 0.0) {
                val alleles = countAlleles(set)
                val total = alleleCounts.sum()
                val minor = if (alleles > 1) alleleCounts[alleleOrder[1]] else 0
                if (options.minMac > 0 && minor < options.minMac) {
                    return fail(into, SkipReason.MinorAlleleCount, minorObs)
                }
                if (options.minMaf > 0.0 && (total == 0 || minor.toDouble() / total < options.minMaf)) {
                    return fail(into, SkipReason.MinorAlleleFrequency, minorObs)
                }
            }

            if (options.minClassSize > 0) {
                var bigEnough = 0
                for (c in classCounts.indices) if (classCounts[c] >= options.minClassSize) bigEnough++
                if (bigEnough < 2) return fail(into, SkipReason.ClassSize, minorObs)
                if (bigEnough < distinct) {
                    for (i in 0 until set.size) {
                        val c = classOf[i]
                        if (c >= 0 && classCounts[c] < options.minClassSize) classOf[i] = -1
                    }
                    for (c in classCounts.indices) if (classCounts[c] < options.minClassSize) classCounts[c] = 0
                    distinct = bigEnough
                    present = classCounts.sum()
                }
            }
            if (distinct < 2) {
                // One class is still a model with no marker columns, and TASSEL 5's GLM reports it (marker df 0), so
                // the observations that have a call and the class stay readable beside the skip reason.
                var kept = 0
                for (i in 0 until set.size) if (classOf[i] >= 0) into.observation[kept++] = i
                into.kept(kept, minorObs)
                if (options.encoding == MarkerEncoding.genotype) encodeGenotype(set, into)
                into.flagSkipped(SkipReason.Monomorphic)
                return false
            }

            var k = 0
            for (i in 0 until set.size) {
                if (classOf[i] >= 0) into.observation[k++] = i
            }
            into.kept(k, minorObs)
            check(k == present)
            return when (options.encoding) {
                MarkerEncoding.genotype -> encodeGenotype(set, into)
                else -> encodeAdditive(set, into)
            }
        }

        private fun encodeGenotype(set: ObservationSet, into: MarkerColumns): Boolean {
            val order = if (options.classOrder == ClassOrder.Glm) GenotypeClasses.glmOrder else GenotypeClasses.mlmOrder
            val labels = ArrayList<String>(MAX_CLASSES_HINT)
            val alleleLabels = ArrayList<String>(MAX_CLASSES_HINT)
            val counts = IntArray(classCounts.count { it > 0 })
            // A class can hold several genotypes (every indel heterozygote is `0`), so its allele string is the first
            // one observed, with the smaller allele code first as TASSEL 5 sorts unphased genotypes.
            firstGenotype.fill(-1)
            for (k in 0 until into.size) {
                val c = classOf[into.observation[k]]
                if (firstGenotype[c] < 0) {
                    val position = set.blockTaxon[into.observation[k]]
                    firstGenotype[c] = siteGenotypes[position].toInt() and 0xFF
                }
            }
            var next = 0
            for (c in order) {
                if (classCounts[c] == 0) continue
                newIndex[c] = next
                labels.add(GenotypeClasses.label[c])
                alleleLabels.add(Nucleotide.genotypeString(GenotypeClasses.unphasedSorted(firstGenotype[c].toByte())))
                counts[next++] = classCounts[c]
            }
            into.classes(labels, alleleLabels, counts)
            for (k in 0 until into.size) into.classIndex[k] = newIndex[classOf[into.observation[k]]]

            // TASSEL 5's add_F inputs: every allele but the last of the frequency order over the genotyped taxa,
            // counting allele copies (TASSEL 5 codes every heterozygote as 1; the two agree at biallelic sites).
            val alleles = (stats ?: AlleleStats.of(siteGenotypes).also { stats = it }).alleles
            if (alleles.size < 2) return true
            val columns = alleles.subList(0, alleles.size - 1)
            into.columns(
                columns.map { Nucleotide.alleleString(it) },
                columns.toByteArray(),
                alleles.last(),
            )
            for ((j, allele) in columns.withIndex()) fillCopies(set, into, into.column(j), allele)
            return true
        }

        private fun encodeAdditive(set: ObservationSet, into: MarkerColumns): Boolean {
            val alleles = countAlleles(set)
            if (alleles < 2) return fail(into, SkipReason.Monomorphic, into.minorObs)
            val dominance = options.encoding == MarkerEncoding.additiveDominance
            if (alleles > 2 && (dominance || options.multiallelic == MultiallelicAdditive.Skip)) {
                return fail(into, SkipReason.Multiallelic, into.minorObs)
            }
            val major = alleleOrder[0].toByte()
            if (alleles > 2) {
                val effects = ByteArray(alleles - 1) { alleleOrder[it + 1].toByte() }
                into.columns(effects.map { Nucleotide.alleleString(it) }, effects, major)
                for (j in effects.indices) fillCopies(set, into, into.column(j), effects[j])
                return true
            }
            val minor = alleleOrder[1].toByte()
            val effect = when (val choice = options.effectAllele) {
                EffectAllele.Minor -> minor
                EffectAllele.Major -> major
                is EffectAllele.Named -> choice.allele
            }
            if (effect != major && effect != minor) return fail(into, SkipReason.EffectAlleleAbsent, into.minorObs)
            val comparison = if (effect == major) minor else major
            val labels = if (dominance) {
                listOf(Nucleotide.alleleString(effect), DOMINANCE_LABEL)
            } else {
                listOf(Nucleotide.alleleString(effect))
            }
            into.columns(
                labels,
                if (dominance) byteArrayOf(effect, Allele.UNKNOWN) else byteArrayOf(effect),
                comparison,
            )
            fillCopies(set, into, into.column(0), effect)
            if (dominance) {
                val heterozygote = into.column(1)
                for (k in 0 until into.size) {
                    val genotype = siteGenotypes[set.blockTaxon[into.observation[k]]]
                    heterozygote[k] = if (Genotype.isHeterozygous(genotype)) 1.0 else 0.0
                }
            }
            return true
        }

        private fun fillCopies(set: ObservationSet, into: MarkerColumns, column: DoubleArray, allele: Byte) {
            for (k in 0 until into.size) {
                val genotype = siteGenotypes[set.blockTaxon[into.observation[k]]]
                var copies = 0
                if (Genotype.firstAllele(genotype) == allele) copies++
                if (Genotype.secondAllele(genotype) == allele) copies++
                column[k] = copies.toDouble()
            }
        }

        /**
         * Counts allele copies over the observations with a class into [alleleCounts], and orders the alleles present
         * into [alleleOrder] by count, highest first, ties by code (`AlleleStats` order). Returns how many are present.
         */
        private fun countAlleles(set: ObservationSet): Int {
            alleleCounts.fill(0)
            for (i in 0 until set.size) {
                if (classOf[i] < 0) continue
                val genotype = siteGenotypes[set.blockTaxon[i]]
                val first = Genotype.firstAllele(genotype).toInt()
                val second = Genotype.secondAllele(genotype).toInt()
                if (first < Allele.NUCLEOTIDE_STATES) alleleCounts[first]++
                if (second < Allele.NUCLEOTIDE_STATES) alleleCounts[second]++
            }
            var present = 0
            for (allele in alleleCounts.indices) {
                if (alleleCounts[allele] == 0) continue
                var at = present++
                while (at > 0 && alleleCounts[alleleOrder[at - 1]] < alleleCounts[allele]) {
                    alleleOrder[at] = alleleOrder[at - 1]
                    at--
                }
                alleleOrder[at] = allele
            }
            return present
        }

        private fun fail(into: MarkerColumns, reason: SkipReason, minorObs: Int = 0): Boolean {
            into.skipped(reason, minorObs)
            return false
        }
    }

    companion object {
        private const val DOSAGE_LABEL = "dosage"
        private const val DOMINANCE_LABEL = "dom"
        private const val MAX_CLASSES_HINT = 8

        /** A source for the `genotype`, `additive`, and `additive-dominance` encodings. */
        fun of(data: GenotypePhenotype): GenotypeMarkerSource = GenotypeMarkerSource(data, null)

        /**
         * A source that also serves the `dosage` encoding from `ReferenceProbability`, which reads the whole
         * genotype table once to find each site's major allele. A missing genotype is a missing dosage (not the taxon
         * mean PCA imputes), as in TASSEL 5's GLM, so the taxon is dropped at that site.
         */
        suspend fun withDosage(
            data: GenotypePhenotype,
            progress: ProgressReporter = ProgressReporter.None,
        ): GenotypeMarkerSource =
            GenotypeMarkerSource(data, ReferenceProbability.of(data.genotypes, progress, imputeMissing = false))
    }
}

/**
 * Genotype classes as TASSEL 5's GLM sees them: the IUPAC string of the genotype. Class ids follow GLM order
 * (`ACGT` first, then the rest by character, as `ModelEffectUtils.getIntegerLevelsSortedByGenotype` sorts), so GLM
 * order is id order; [mlmOrder] is MLM's (`ModelEffectUtils.getIntegerLevels(byte[], ...)`).
 */
private object GenotypeClasses {
    /** Genotype byte to class id, or `-1` when an allele is unknown. */
    val classOf = IntArray(256)

    /** Genotype byte to whether both alleles are A, C, G, T, insert, or gap. */
    val scorable = BooleanArray(256)

    /** TASSEL 5 IUPAC string of each class. */
    val label: List<String>

    /** Whether each class is a heterozygote string (`R S Y W K M 0`), as TASSEL 5's filters test it. */
    val heterozygous: BooleanArray

    /** Class ids in GLM order. */
    val glmOrder: IntArray

    /** Class ids in MLM order. */
    val mlmOrder: IntArray

    val size: Int

    private const val NUCLEOTIDES = "ACGT"
    private const val HETEROZYGOTES = "RSYWKM0"

    init {
        val unknown = Allele.UNKNOWN.toInt()
        val members = HashMap<Char, MutableList<Byte>>()
        for (g in 0..255) {
            if ((g shr 4) == unknown || (g and 0xF) == unknown) continue
            members.getOrPut(Nucleotide.iupac(g.toByte())) { ArrayList() }.add(unphasedSorted(g.toByte()))
        }
        val chars = members.keys.sortedWith(
            compareBy<Char>({ if (it in NUCLEOTIDES) 0 else 1 }, { it }),
        )
        size = chars.size
        val id = chars.withIndex().associate { (index, char) -> char to index }
        for (g in 0..255) {
            classOf[g] =
                if ((g shr 4) == unknown || (g and 0xF) == unknown) -1 else id.getValue(Nucleotide.iupac(g.toByte()))
            scorable[g] = (g shr 4) < Allele.NUCLEOTIDE_STATES && (g and 0xF) < Allele.NUCLEOTIDE_STATES
        }
        label = chars.map { it.toString() }
        heterozygous = BooleanArray(size) { chars[it] in HETEROZYGOTES }
        // The representative is the homozygote TASSEL 5's MLM would list first if the class has one, else the
        // smallest unphased genotype byte (Java bytes are signed).
        val representative = chars.map { char ->
            members.getValue(char).minWith(compareBy<Byte>({ if (isListedHomozygote(it)) 0 else 1 }, { it.toInt() }))
        }
        glmOrder = IntArray(size) { it }
        mlmOrder = (0 until size)
            .sortedWith(
                compareBy<Int>(
                    { if (isListedHomozygote(representative[it])) 0 else 1 },
                    { representative[it].toInt() },
                ),
            )
            .toIntArray()
    }

    /** TASSEL 5's `homGeno`: AA, CC, GG, TT, ZZ, ++, --. */
    private fun isListedHomozygote(genotype: Byte): Boolean =
        genotype in byteArrayOf(0x00, 0x11, 0x22, 0x33, 0xEE.toByte(), 0x44, 0x55)

    /** `GenotypeTableUtils.getUnphasedSortedDiploidValue`: the smaller allele code first. */
    fun unphasedSorted(genotype: Byte): Byte {
        val first = Genotype.firstAllele(genotype)
        val second = Genotype.secondAllele(genotype)
        return if (first < second) genotype else Genotype.of(second, first)
    }
}
