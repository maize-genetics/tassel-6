package net.maizegenetics.tassel.distance

import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.coroutineScope
import net.maizegenetics.tassel.genotype.AlleleStats
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.SiteBlock
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report

/**
 * The `-maxThreads` TASSEL 5 ran with when the goldens were generated (see `testdata/golden/MANIFEST`).
 *
 * TASSEL 5 splits sites into blocks of `max(numberOfSites / maxThreads, 1000)`, each summing into its
 * own `float` array, and then adds those arrays together. The block size therefore fixes the order the
 * single-precision sums are made in, and so the last bit or two of every value. Taking the thread count
 * from the host would make the output depend on the machine, so this kernel uses the 4 the goldens were
 * made with. Extra cores come from disjoint taxa bands inside a block, which do not reorder any addition.
 */
private const val GOLDEN_MAX_THREADS = 4

/** Floor on the block size, from TASSEL 5's `max(numberOfSites / maxThreads, 1000)`. */
private const val MIN_SITES_PER_BLOCK = 1000

/** Pseudo-sites TASSEL packs into one allele-count word: five three-bit counts in a 15-bit short. */
private const val PSEUDO_SITES_PER_WORD = 5

/** Allele-count words TASSEL sums together before adding to the running total: three, so 15 pseudo-sites. */
private const val WORDS_PER_GROUP = 3

/** Terms per pseudo-site in the term table: one per three-bit allele-count code. */
private const val TERMS_PER_PSEUDO_SITE = 8

/** All five pseudo-site slots of a word set to the "not comparable" code. */
private const val NOTHING_COMPARABLE = 0x7FFF

/** Progress reports per block, following TASSEL 5's `(fence - start) / 10`. */
private const val REPORTS_PER_BLOCK = 10

/** Floor on the sites per progress report, following TASSEL 5. */
private const val MIN_SITES_PER_REPORT = 100

/** Groups accumulated before a parallel join. Small enough that the count tables stay modest. */
private const val BATCH_GROUPS = 16

private const val TERMS_PER_WORD = PSEUDO_SITES_PER_WORD * TERMS_PER_PSEUDO_SITE

/**
 * Bit-exact Centered IBS, the transcription of TASSEL 5's `EndelmanDistanceMatrix` that
 * [centeredIbs] replaced as the default (ADR-0020).
 *
 * Sums are accumulated in [Float] exactly where TASSEL 5 accumulates them. [workers] partitions
 * each block's triangle into private slices; each cell still adds groups in group order, and the
 * parent copies the slices back after the block. `workers = 1` is the sequential path.
 *
 * @see net.maizegenetics.analysis.distance.EndelmanDistanceMatrix
 */
internal suspend fun centeredIbsTassel5Order(
    table: GenotypeTable,
    maxAlleles: Int = DEFAULT_MAX_ALLELES,
    progress: ProgressReporter = ProgressReporter.None,
    workers: Int = kinshipWorkerCount(),
): DistanceMatrix {
    require(maxAlleles in 2..DEFAULT_MAX_ALLELES) {
        "centeredIbs: maxAlleles must be between 2 and $DEFAULT_MAX_ALLELES inclusive: $maxAlleles"
    }
    val numberOfTaxa = table.numberOfTaxa
    var sitesDone = 0
    val blockSize = maxOf(table.numberOfSites / GOLDEN_MAX_THREADS, MIN_SITES_PER_BLOCK)
    val blocks = buildList {
        for (block in table.siteBlocks(blockSize = blockSize)) {
            add(
                accumulate(block, maxAlleles, workers) { sites ->
                    sitesDone += sites
                    progress.report(sitesDone, table.numberOfSites)
                },
            )
        }
    }
    if (blocks.isEmpty()) {
        return finishCenteredIbs(table, FloatArray(DistanceMatrix.triangleSize(numberOfTaxa)), 0.0)
    }
    // TASSEL 5 reduces its parallel stream right to left, so the blocks are added as a + (b + (c + d)).
    val total = blocks.reduceRight { block, rest -> block.also { it.add(rest) } }
    return finishCenteredIbs(table, total.distances, total.sumPi * 2.0)
}

/**
 * Partial sums over a run of sites: the unscaled pairwise sums in the same upper-triangle order
 * TASSEL 5 fills them (`(0,0), (0,1), ... (0,n-1), (1,1), ...`), and the allele-frequency sum.
 */
private class Sums(numberOfTaxa: Int) {
    val distances = FloatArray(DistanceMatrix.triangleSize(numberOfTaxa))
    var sumPi: Double = 0.0

    fun add(other: Sums) {
        for (i in distances.indices) distances[i] += other.distances[i]
        sumPi += other.sumPi
    }
}

/**
 * One group of allele-count words, reused across a batch so a parallel window does not allocate
 * per site. [size] is how many groups [counts] and [terms] currently hold.
 */
private class GroupBatch(numberOfTaxa: Int, val capacity: Int) {
    val counts = Array(capacity) { Array(WORDS_PER_GROUP) { ShortArray(numberOfTaxa) } }
    val terms = Array(capacity) { Array(WORDS_PER_GROUP) { FloatArray(TERMS_PER_WORD) } }
    var size: Int = 0

    fun clear() {
        size = 0
    }
}

/**
 * Sums one block of sites into a triangle of its own, reporting the sites covered as it goes.
 * Genotypes are copied into one buffer for the block; allele statistics are counted from that buffer.
 */
private suspend fun accumulate(
    block: SiteBlock,
    maxAlleles: Int,
    workers: Int,
    reportSites: suspend (Int) -> Unit,
): Sums {
    val numberOfTaxa = block.numberOfTaxa
    val sums = Sums(numberOfTaxa)
    val fence = block.sites.last + 1
    val sitesPerReport = maxOf((fence - block.sites.first) / REPORTS_PER_BLOCK, MIN_SITES_PER_REPORT)
    val genotypes = ByteArray(numberOfTaxa)
    val bands = pairBands(numberOfTaxa, workers.coerceAtLeast(1))
    val parallel = bands.size > 1
    val slices = if (parallel) bands.map { FloatArray(pairCount(numberOfTaxa, it)) } else null
    val batch = if (parallel) GroupBatch(numberOfTaxa, BATCH_GROUPS) else null
    val counts = Array(WORDS_PER_GROUP) { ShortArray(numberOfTaxa) }
    val terms = Array(WORDS_PER_GROUP) { FloatArray(TERMS_PER_WORD) }

    var site = block.sites.first
    while (site < fence) {
        val reportFence = minOf(site + sitesPerReport, fence)
        val reportStart = site
        if (batch == null || slices == null) {
            while (site < reportFence) {
                val consumed = fillGroup(block, site, fence, maxAlleles, counts, terms, genotypes, sums)
                site += consumed
                addGroupRows(sums.distances, numberOfTaxa, 0 until numberOfTaxa, counts, terms)
            }
        } else {
            while (site < reportFence) {
                batch.clear()
                while (site < reportFence && batch.size < batch.capacity) {
                    val slot = batch.size
                    val consumed = fillGroup(
                        block,
                        site,
                        fence,
                        maxAlleles,
                        batch.counts[slot],
                        batch.terms[slot],
                        genotypes,
                        sums,
                    )
                    site += consumed
                    batch.size++
                }
                addBatch(slices, bands, batch, numberOfTaxa)
            }
        }
        reportSites(site - reportStart)
    }
    if (slices != null) {
        var destination = 0
        for (slice in slices) {
            slice.copyInto(sums.distances, destination)
            destination += slice.size
        }
        check(destination == sums.distances.size) {
            "taxa bands covered $destination pairs, expected ${sums.distances.size}"
        }
    }
    return sums
}

/** Fills the three words of one group starting at [site] and returns the sites consumed. */
private fun fillGroup(
    block: SiteBlock,
    site: Int,
    fence: Int,
    maxAlleles: Int,
    counts: Array<ShortArray>,
    terms: Array<FloatArray>,
    genotypes: ByteArray,
    sums: Sums,
): Int {
    var consumed = 0
    for (word in 0 until WORDS_PER_GROUP) {
        consumed += fillWord(
            block,
            site + consumed,
            fence,
            maxAlleles,
            counts[word],
            terms[word],
            genotypes,
            sums,
        )
    }
    return consumed
}

/**
 * Fills one allele-count word and its term table from the sites at [startSite] and after, and adds
 * their allele frequencies to [sums]. Returns the number of real sites consumed.
 *
 * A site contributes one pseudo-site per allele evaluated, and a site's alleles are never split
 * across words, so a site that would not fit is left for the next word. [counts] holds, per taxon,
 * five three-bit codes saying how many copies of that pseudo-site's allele the taxon carries: 1 for
 * none, 2 for one, 4 for two, and 7 for not comparable. [terms] holds, per pseudo-site and per pair
 * of codes OR-ed together, the term that pair contributes.
 */
private fun fillWord(
    block: SiteBlock,
    startSite: Int,
    fence: Int,
    maxAlleles: Int,
    counts: ShortArray,
    terms: FloatArray,
    genotypes: ByteArray,
    sums: Sums,
): Int {
    counts.fill(NOTHING_COMPARABLE.toShort())
    terms.fill(0.0f)
    var pseudoSite = 0
    var site = startSite
    var consumed = 0
    while (pseudoSite < PSEUDO_SITES_PER_WORD && site < fence) {
        block.copySiteGenotypes(site, genotypes)
        val stats = AlleleStats.of(genotypes)
        val evaluated = minOf(stats.alleles.size - 1, maxAlleles - 1)
        if (evaluated + pseudoSite > PSEUDO_SITES_PER_WORD) return consumed
        if (evaluated > 0) {
            val gametes = stats.alleleGameteCount
            for (alleleIndex in 0 until evaluated) {
                val allele = stats.alleles[alleleIndex]
                val frequency = stats.alleleCounts[alleleIndex].toFloat() / gametes.toFloat()
                val twiceFrequency = frequency * 2.0f
                sums.sumPi += frequency * (1.0 - frequency)
                addPseudoSite(pseudoSite, allele, twiceFrequency, genotypes, counts, terms)
                pseudoSite++
            }
        }
        site++
        consumed++
    }
    return consumed
}

/** Writes one pseudo-site's six terms and every taxon's allele count into slot [pseudoSite]. */
private fun addPseudoSite(
    pseudoSite: Int,
    allele: Byte,
    twiceFrequency: Float,
    genotypes: ByteArray,
    counts: ShortArray,
    terms: FloatArray,
) {
    val none = 0.0f - twiceFrequency
    val one = 1.0f - twiceFrequency
    val two = 2.0f - twiceFrequency
    val base = pseudoSite * TERMS_PER_PSEUDO_SITE
    terms[base + 1] = none * none
    terms[base + 3] = none * one
    terms[base + 5] = none * two
    terms[base + 2] = one * one
    terms[base + 6] = one * two
    terms[base + 4] = two * two

    val shift = countShift(pseudoSite)
    val keep = (0x7 shl shift).inv() and NOTHING_COMPARABLE
    val codes = CODE_BY_ALLELE_AND_GENOTYPE[allele.toInt() and 0x7]
    for (taxon in counts.indices) {
        val code = codes[genotypes[taxon].toInt() and 0xFF].toInt()
        counts[taxon] = (counts[taxon].toInt() and (keep or (code shl shift))).toShort()
    }
}

/** Adds one batch of groups into each band's private slice. The slices outlive the batch. */
private suspend fun addBatch(slices: List<FloatArray>, bands: List<IntRange>, batch: GroupBatch, numberOfTaxa: Int) {
    if (batch.size == 0) return
    coroutineScope {
        bands.mapIndexed { band, rows ->
            async {
                val dest = slices[band]
                for (group in 0 until batch.size) {
                    addGroupRows(dest, numberOfTaxa, rows, batch.counts[group], batch.terms[group])
                }
            }
        }.awaitAll()
    }
}

/** Adds one group of 15 pseudo-sites to the pairs in [rows], in TASSEL 5's summation order. */
private fun addGroupRows(
    distances: FloatArray,
    numberOfTaxa: Int,
    rows: IntRange,
    counts: Array<ShortArray>,
    terms: Array<FloatArray>,
) {
    val counts0 = counts[0]
    val counts1 = counts[1]
    val counts2 = counts[2]
    var index = 0
    for (first in rows) {
        val word0 = counts0[first].toInt()
        val word1 = counts1[first].toInt()
        val word2 = counts2[first].toInt()
        if (word0 == NOTHING_COMPARABLE && word1 == NOTHING_COMPARABLE && word2 == NOTHING_COMPARABLE) {
            // Every pseudo-site of this taxon is unknown, so every term it contributes would be zero.
            index += numberOfTaxa - first
            continue
        }
        for (second in first until numberOfTaxa) {
            distances[index] += (
                wordSum(terms[0], word0 or counts0[second].toInt()) +
                    wordSum(terms[1], word1 or counts1[second].toInt())
                ) + wordSum(terms[2], word2 or counts2[second].toInt())
            index++
        }
    }
}

/**
 * Sums the five pseudo-site terms a pair's OR-ed allele-count [word] selects, in slot order. The five
 * additions are written out rather than looped so the order and the single-precision rounding are
 * TASSEL 5's, down to the sign of a zero.
 */
private fun wordSum(terms: FloatArray, word: Int): Float {
    var sum = terms[termIndex(0, word)]
    sum += terms[termIndex(1, word)]
    sum += terms[termIndex(2, word)]
    sum += terms[termIndex(3, word)]
    sum += terms[termIndex(4, word)]
    return sum
}

/** Where pseudo-site [slot]'s term for the allele-count code held in [word] sits in a term table. */
private fun termIndex(slot: Int, word: Int): Int = slot * TERMS_PER_PSEUDO_SITE + ((word shr countShift(slot)) and 0x7)

/** Bit position of pseudo-site [slot]'s three-bit allele-count code; slot 0 is the highest. */
private fun countShift(slot: Int): Int = (PSEUDO_SITES_PER_WORD - slot - 1) * 3
