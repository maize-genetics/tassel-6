package net.maizegenetics.tassel.distance

import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.SiteBlock
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.table.doubleToString

/** TASSEL 5's `KinshipPlugin` default for the number of allele states per site to evaluate. */
const val DEFAULT_MAX_ALLELES: Int = 6

/**
 * The `-maxThreads` TASSEL 5 ran with when the goldens were generated (see `testdata/golden/MANIFEST`).
 *
 * TASSEL 5 splits sites into blocks of `max(numberOfSites / maxThreads, 1000)`, each summing into its
 * own `float` array, and then adds those arrays together. The block size therefore fixes the order the
 * single-precision sums are made in, and so the last bit or two of every value. Taking the thread count
 * from the host would make the output depend on the machine, so this port uses the 4 the goldens were
 * made with and stays single-threaded.
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

/**
 * Centered IBS (Endelman) kinship for every pair of taxa in [table]. Missing genotypes are ignored.
 *
 * The matrix carries TASSEL 5's two annotations: `Matrix_Type=Centered_IBS` and `Centered_IBS.SumPk`,
 * the sum of `2 * sum(p * (1 - p))` over every allele evaluated, which is what the raw sums are
 * divided by. [maxAlleles] is how many allele states per site to evaluate, 2 to 6 inclusive.
 *
 * Sums are accumulated in [Float] exactly where TASSEL 5 accumulates them, which is the precision
 * ceiling of its output and the reason the golden tolerance is 8e-6. Sites are read one block at a
 * time; the working set is one block of genotypes plus at most [GOLDEN_MAX_THREADS] + 1 triangles of
 * partial sums.
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
    var sitesDone = 0
    val blockSize = maxOf(table.numberOfSites / GOLDEN_MAX_THREADS, MIN_SITES_PER_BLOCK)
    // A plain loop rather than Sequence.map: reporting progress suspends (ADR-0019), and a lazy
    // sequence's transform cannot.
    val blocks = buildList {
        for (block in table.siteBlocks(blockSize = blockSize)) {
            add(
                accumulate(block, maxAlleles) { sites ->
                    sitesDone += sites
                    progress.report(sitesDone, table.numberOfSites)
                },
            )
        }
    }
    // TASSEL 5 reduces its parallel stream right to left, so the blocks are added as a + (b + (c + d)).
    // Any other order shifts the float sums by a few ulps, which the golden text would show.
    val total = blocks.reduceRight { block, rest -> block.also { it.add(rest) } }

    val sumPk = total.sumPi * 2.0
    val scaled = FloatArray(DistanceMatrix.triangleSize(numberOfTaxa))
    var index = 0
    for (first in 0 until numberOfTaxa) {
        for (second in first until numberOfTaxa) {
            scaled[second * (second + 1) / 2 + first] = (total.distances[index] / sumPk).toFloat()
            index++
        }
    }
    return DistanceMatrix(
        taxa = table.taxa,
        lowerTriangle = scaled,
        annotations = mapOf(
            DistanceMatrix.MATRIX_TYPE to KinshipMethod.Centered_IBS.name,
            DistanceMatrix.CENTERED_IBS_SUMPK to doubleToString(sumPk),
        ),
    )
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
 * Sums one block of sites into a triangle of its own, reporting the sites covered as it goes. Each
 * site's genotypes are read once.
 */
private suspend fun accumulate(block: SiteBlock, maxAlleles: Int, reportSites: suspend (Int) -> Unit): Sums {
    val numberOfTaxa = block.numberOfTaxa
    val sums = Sums(numberOfTaxa)
    val fence = block.sites.last + 1
    val sitesPerReport = maxOf((fence - block.sites.first) / REPORTS_PER_BLOCK, MIN_SITES_PER_REPORT)
    val counts = Array(WORDS_PER_GROUP) { ShortArray(numberOfTaxa) }
    val terms = Array(WORDS_PER_GROUP) { FloatArray(PSEUDO_SITES_PER_WORD * TERMS_PER_PSEUDO_SITE) }

    var site = block.sites.first
    while (site < fence) {
        val reportFence = minOf(site + sitesPerReport, fence)
        val reportStart = site
        while (site < reportFence) {
            // A group may run past reportFence, as it does in TASSEL 5: only the block fence stops it.
            var sites = 0
            for (word in 0 until WORDS_PER_GROUP) {
                sites += fillWord(block, site + sites, fence, maxAlleles, counts[word], terms[word], sums)
            }
            site += sites
            addGroup(sums.distances, numberOfTaxa, counts, terms)
        }
        reportSites(site - reportStart)
    }
    return sums
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
    sums: Sums,
): Int {
    counts.fill(NOTHING_COMPARABLE.toShort())
    terms.fill(0.0f)
    var pseudoSite = 0
    var site = startSite
    var consumed = 0
    while (pseudoSite < PSEUDO_SITES_PER_WORD && site < fence) {
        val stats = block.siteStats(site)
        val evaluated = minOf(stats.alleles.size - 1, maxAlleles - 1)
        if (evaluated + pseudoSite > PSEUDO_SITES_PER_WORD) return consumed
        if (evaluated > 0) {
            val genotypes = block.siteGenotypes(site)
            val gametes = stats.alleleGameteCount
            for (a in 0 until evaluated) {
                val allele = stats.alleles[a]
                val frequency = stats.alleleCounts[a].toFloat() / gametes.toFloat()
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
    val alleleBits = (allele.toInt() and 0x7) shl 6
    for (taxon in counts.indices) {
        val genotype = genotypes[taxon].toInt()
        val code = COUNT_CODES[alleleBits or ((genotype and 0x70) shr 1) or (genotype and 0x7)].toInt()
        counts[taxon] = (counts[taxon].toInt() and (keep or (code shl shift))).toShort()
    }
}

/** Adds one group of 15 pseudo-sites to the running pairwise sums, in TASSEL 5's summation order. */
private fun addGroup(distances: FloatArray, numberOfTaxa: Int, counts: Array<ShortArray>, terms: Array<FloatArray>) {
    val counts0 = counts[0]
    val counts1 = counts[1]
    val counts2 = counts[2]
    var index = 0
    for (first in 0 until numberOfTaxa) {
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

/**
 * How many copies of an allele each diploid genotype carries, indexed by
 * `allele shl 6 or firstAllele shl 3 or secondAllele` with the three low bits of each allele code.
 * The values are the codes [addPseudoSite] packs into an allele-count word: 7 when the genotype or
 * the allele is unknown, 1 for no copies, 2 for one, and 4 for two. OR-ing two codes together gives
 * the pair a term-table index, and OR-ing anything with 7 lands on the term table's zero.
 */
private val COUNT_CODES = ByteArray(512) { index ->
    val allele = index shr 6
    val firstAllele = (index shr 3) and 0x7
    val secondAllele = index and 0x7
    val unknown = 0x7
    when {
        allele == unknown || (firstAllele == unknown && secondAllele == unknown) -> 7
        firstAllele == allele && secondAllele == allele -> 4
        firstAllele == allele || secondAllele == allele -> 2
        else -> 1
    }.toByte()
}
