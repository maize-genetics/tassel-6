package net.maizegenetics.tassel.analysis

import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.SiteBlock
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import kotlin.math.floor

/**
 * TASSEL 5's byte encoding for a site score in `[0, 1]`.
 *
 * A reference probability is not kept as the `float` it is computed from: `ReferenceProbabilityBuilder`
 * puts it through [encode], which quantizes it to one of 255 steps, and reading it back through
 * [decode] gives that step rounded to three decimals. **That byte, not the `float`, is the precision
 * ceiling of PCA on genotypes**, and reproducing it is what makes the eigenvalues match TASSEL 5 to
 * 1e-13 instead of 1e-3.
 *
 * The three values an observed genotype can take (0.0, 0.5, and 1.0) are fixed points of the round
 * trip, so only an imputed value actually moves.
 *
 * @see net.maizegenetics.dna.snp.score.SiteScoreUtil
 * @see net.maizegenetics.dna.snp.score.ReferenceProbabilityBuilder
 */
internal object SiteScore {
    /** The byte that stands for a missing score. */
    const val NAN: Int = 255

    // SiteScoreUtil.BYTE_TO_FLOAT, which TASSEL 5 fills once in a static initializer.
    private val BYTE_TO_FLOAT = FloatArray(256) { if (it == NAN) Float.NaN else step(it) }

    /** [value] as the byte TASSEL 5 stores, in `0..255`; [NAN] for `NaN`. */
    fun encode(value: Float): Int {
        if (value.isNaN()) return NAN
        require(value >= 0.0f && value <= 1.0f) { "A site score must be between 0.0 and 1.0: $value" }
        // Java's Math.round(float) is floor(x + 0.5) evaluated exactly, which a Double gives here.
        return floor((254.0f * value).toDouble() + 0.5).toInt()
    }

    /** The score stored in byte [code], which is `0..255`. */
    fun decode(code: Int): Float = BYTE_TO_FLOAT[code]

    /** [value] after the round trip through [encode] and [decode]. */
    fun quantize(value: Float): Float = decode(encode(value))

    /**
     * `SiteScoreUtil.decode`: `code / 254f`, then that float's own decimal spelling rounded half up
     * to three places. `floor(x * 1000 + 0.5)` is that rounding exactly, because a `Float` scaled by
     * 1000 needs 31 significant bits and so is exact in a `Double`.
     */
    private fun step(code: Int): Float = (floor((code.toFloat() / 254.0f).toDouble() * 1000.0 + 0.5) / 1000.0).toFloat()
}

/**
 * The numeric matrix PCA runs on: for every taxon and site, the probability that an allele drawn at
 * random there is the site's major allele. Homozygous major is 1.0, homozygous minor 0.0, and
 * heterozygous 0.5; a missing genotype takes the mean of its **taxon**, quantized by [SiteScore].
 *
 * Imputing by taxon rather than by site is deliberate, and is TASSEL 5's behavior rather than its
 * documented intent: `ImputationPlugin` hands `ImputationByMean.impute` a sites-by-taxa matrix, and
 * that method averages each *column*. Imputing by site instead moves the first eigenvalue by 7.6%,
 * so it is not a detail that could be left either way.
 *
 * Columns are produced one site at a time from the site blocks, so no taxa-by-sites matrix is held.
 *
 * @see net.maizegenetics.analysis.numericaltransform.NumericalGenotypePlugin.setAlternateMinorAllelesToMinor
 * @see net.maizegenetics.dna.snp.GenotypeTableUtils.convertGenotypeToFloatProbability
 * @see net.maizegenetics.analysis.numericaltransform.ImputationByMean
 */
internal class ReferenceProbability private constructor(
    /** The major allele of each site, by site index. */
    private val majorAlleles: ByteArray,
    /** The value a missing genotype takes, by taxon index. */
    private val imputed: FloatArray,
) {
    /** Fills [into] with site [site]'s probabilities, in taxon order. */
    fun column(block: SiteBlock, site: Int, into: DoubleArray) {
        val major = majorAlleles[site]
        for (taxon in into.indices) {
            val genotype = block.genotype(taxon, site)
            into[taxon] = if (genotype == Genotype.UNKNOWN) {
                imputed[taxon].toDouble()
            } else {
                probability(genotype, major).toDouble()
            }
        }
    }

    companion object {
        /**
         * Reads [genotypes] once to find each site's major allele and each taxon's mean, which is
         * the first of the three passes PCA makes over the data.
         */
        fun of(genotypes: GenotypeTable, progress: ProgressReporter = ProgressReporter.None): ReferenceProbability {
            val numberOfTaxa = genotypes.numberOfTaxa
            val numberOfSites = genotypes.numberOfSites
            val majorAlleles = ByteArray(numberOfSites)
            val sums = DoubleArray(numberOfTaxa)
            val counts = DoubleArray(numberOfTaxa)
            var done = 0
            genotypes.siteBlocks().forEach { block ->
                for (site in block.sites) {
                    val major = block.siteStats(site).majorAllele
                    majorAlleles[site] = major
                    for (taxon in 0 until numberOfTaxa) {
                        val genotype = block.genotype(taxon, site)
                        if (genotype == Genotype.UNKNOWN) continue
                        sums[taxon] += probability(genotype, major).toDouble()
                        counts[taxon] += 1.0
                    }
                }
                done += block.sites.count()
                progress.report(done, numberOfSites, MESSAGE)
            }
            // ImputationByMean writes the column mean as a float, which ReferenceProbabilityBuilder
            // then quantizes. A taxon with no called genotype keeps TASSEL 5's 0.0 / 0.0 = NaN.
            val imputed = FloatArray(numberOfTaxa) { SiteScore.quantize((sums[it] / counts[it]).toFloat()) }
            return ReferenceProbability(majorAlleles, imputed)
        }

        /** Progress message for the pass [of] makes. */
        const val MESSAGE: String = "PCA: reference probabilities"

        /**
         * 0.5 per allele copy that equals [major]. An allele of a partly missing genotype such as
         * `A/N` is simply not the major allele, which is why only a wholly missing genotype is
         * imputed; TASSEL 5 does the same, by adding to a cell it has already set to `NaN`.
         */
        private fun probability(genotype: Byte, major: Byte): Float {
            var probability = 0.0f
            if (Genotype.firstAllele(genotype) == major) probability += 0.5f
            if (Genotype.secondAllele(genotype) == major) probability += 0.5f
            return probability
        }
    }
}
