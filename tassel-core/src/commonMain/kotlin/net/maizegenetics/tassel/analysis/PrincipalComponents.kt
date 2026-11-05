package net.maizegenetics.tassel.analysis

import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.math.SymmetricEigenDecomposition
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.table.PhenotypeTable
import net.maizegenetics.tassel.table.ResultTable
import net.maizegenetics.tassel.table.resultTable
import kotlin.math.max
import kotlin.math.sqrt

/** TASSEL 5's default number of principal components. */
const val DEFAULT_COMPONENTS: Int = 5

/**
 * What [principalComponents] produces: the principal components themselves, and the two optional
 * tables TASSEL 5 returns beside them.
 *
 * `Operation` has a single output type, so the three results are wrapped here rather than returned
 * as a collection, following the `GenotypeSummary` pattern.
 */
class PrincipalComponentsResult internal constructor(
    /**
     * One covariate per component, over all taxa, as TASSEL 5 stores them: a `Phenotype` of
     * `Float` values. Always produced.
     */
    val principalComponents: PhenotypeTable,
    /** Every eigenvalue with its share of the total variance, or `null` if it was not requested. */
    val eigenvalues: ResultTable?,
    /** One row per site, with a column per component, or `null` if it was not requested. */
    val eigenvectors: ResultTable?,
) {
    /** The tables produced, in TASSEL 5's order. */
    val tables: List<ResultTable> get() = listOfNotNull(eigenvalues, eigenvectors)
}

/**
 * Principal components of [genotypes], by eigendecomposition of the taxa Gram matrix.
 *
 * ## What the data matrix is
 *
 * PCA does not run on genotype bytes. TASSEL 5 first turns them into `float` reference
 * probabilities, mean-imputes the missing ones, and stores the result as quantized bytes; see
 * [ReferenceProbability] and [SiteScore], which carry that arithmetic. Only then does `PrinComp`
 * center each site column, scale the matrix by `1 / sqrt(nTaxa - 1)`, and decompose it.
 *
 * ## Why the Gram matrix
 *
 * With `A` the centered, scaled taxa-by-sites matrix and `A = U S V'`, the eigenvalues are `S^2`,
 * the site-space eigenvectors TASSEL prints are `V = A' U S^-1`, and the principal components are
 * `sqrt(nTaxa - 1) * U * S`. Eigendecomposing `G = A A'` therefore gives everything from a matrix
 * that is taxa by taxa, never sites by sites: 281 x 281 for the tutorial data rather than
 * 3,093 x 3,093. The data are read in three passes over the site blocks, one to find the major
 * alleles and taxon means, one to accumulate `G`, and one to build `V`.
 *
 * @param numberOfComponents components to return, capped at `min(nTaxa, nSites)`.
 * @param useCovariance decompose the covariance matrix (TASSEL 5's default), or the correlation
 *   matrix when false, which additionally divides each centered column by its standard deviation.
 * @see net.maizegenetics.analysis.data.PrincipalComponentsPlugin
 * @see net.maizegenetics.stats.PCA.PrinComp
 */
suspend fun principalComponents(
    genotypes: GenotypeTable,
    numberOfComponents: Int = DEFAULT_COMPONENTS,
    useCovariance: Boolean = true,
    reportEigenvalues: Boolean = true,
    reportEigenvectors: Boolean = true,
    progress: ProgressReporter = ProgressReporter.None,
): PrincipalComponentsResult {
    val numberOfTaxa = genotypes.numberOfTaxa
    val numberOfSites = genotypes.numberOfSites
    require(numberOfTaxa > 1) { "PCA needs at least two taxa, found $numberOfTaxa" }
    require(numberOfSites > 0) { "PCA needs at least one site" }
    require(numberOfComponents > 0) { "PCA needs at least one component, asked for $numberOfComponents" }

    // TASSEL 5's SVD is compact, so it yields min(nTaxa, nSites) singular values; the Gram matrix
    // has nTaxa eigenvalues, and any beyond that count are zero.
    val values = minOf(numberOfTaxa, numberOfSites)
    val components = minOf(numberOfComponents, values)

    val probabilities = ReferenceProbability.of(genotypes, progress)
    val gram = gramMatrix(genotypes, probabilities, useCovariance, progress)
    val eigen = SymmetricEigenDecomposition.of(numberOfTaxa) { row, column -> gram[row * numberOfTaxa + column] }

    // Round-off can push the smallest eigenvalues of a positive-semidefinite matrix below zero.
    val eigenvalues = DoubleArray(values) { eigen.eigenvalue(it) }
    val singular = DoubleArray(components) { sqrt(max(eigenvalues[it], 0.0)) }

    return PrincipalComponentsResult(
        principalComponents = principalComponentTable(genotypes, eigen, singular),
        eigenvalues = if (reportEigenvalues) eigenvalueTable(eigenvalues) else null,
        eigenvectors = if (reportEigenvectors) {
            eigenvectorTable(genotypes, probabilities, eigen, singular, useCovariance, progress)
        } else {
            null
        },
    )
}

/**
 * `G = A A'`, accumulated one site column at a time, with only the lower triangle filled because
 * that is all [SymmetricEigenDecomposition] reads. Row-major, so `G[row, column]` is at
 * `row * nTaxa + column`.
 */
private suspend fun gramMatrix(
    genotypes: GenotypeTable,
    probabilities: ReferenceProbability,
    useCovariance: Boolean,
    progress: ProgressReporter,
): DoubleArray {
    val numberOfTaxa = genotypes.numberOfTaxa
    val gram = DoubleArray(numberOfTaxa * numberOfTaxa)
    val column = DoubleArray(numberOfTaxa)
    forEachColumn(genotypes, probabilities, useCovariance, column, progress, GRAM_MESSAGE) {
        for (row in 0 until numberOfTaxa) {
            val value = column[row]
            val offset = row * numberOfTaxa
            for (other in 0..row) gram[offset + other] += value * column[other]
        }
    }
    return gram
}

/** `V = A' U S^-1`: one row per site, holding the first [singular]`.size` eigenvectors. */
private suspend fun eigenvectorTable(
    genotypes: GenotypeTable,
    probabilities: ReferenceProbability,
    eigen: SymmetricEigenDecomposition,
    singular: DoubleArray,
    useCovariance: Boolean,
    progress: ProgressReporter,
): ResultTable {
    val numberOfTaxa = genotypes.numberOfTaxa
    val column = DoubleArray(numberOfTaxa)
    val names = List(singular.size) { "Eigenvector${it + 1}" }
    return resultTable(EIGENVECTOR_TITLE, listOf(TRAIT) + names) {
        forEachColumn(genotypes, probabilities, useCovariance, column, progress, EIGENVECTOR_MESSAGE) { site ->
            val cells = ArrayList<Any?>(singular.size + 1)
            cells += genotypes.positions[site].name
            for (component in singular.indices) {
                var projection = 0.0
                for (taxon in 0 until numberOfTaxa) projection += column[taxon] * eigen[taxon, component]
                cells += if (singular[component] == 0.0) 0.0 else projection / singular[component]
            }
            row(cells)
        }
    }
}

/**
 * Every eigenvalue with its share of the total variance, as `PrincipalComponentsPlugin` builds it:
 * the `PC` column counts from zero and is text, and the proportions divide by the running total.
 */
private fun eigenvalueTable(eigenvalues: DoubleArray): ResultTable {
    val cumulative = DoubleArray(eigenvalues.size)
    var running = 0.0
    eigenvalues.forEachIndexed { index, eigenvalue ->
        running += eigenvalue
        cumulative[index] = running
    }
    val total = cumulative[cumulative.size - 1]
    return resultTable(EIGENVALUE_TITLE, EIGENVALUE_COLUMNS) {
        eigenvalues.forEachIndexed { index, eigenvalue ->
            row(index.toString(), eigenvalue, eigenvalue / total, cumulative[index] / total)
        }
    }
}

/** The requested components as TASSEL 5 stores them: `Float` covariates in a `Phenotype`. */
private fun principalComponentTable(
    genotypes: GenotypeTable,
    eigen: SymmetricEigenDecomposition,
    singular: DoubleArray,
): PhenotypeTable {
    val numberOfTaxa = genotypes.numberOfTaxa
    val scale = sqrt((numberOfTaxa - 1).toDouble())
    val attributes = ArrayList<PhenotypeTable.Attribute>(singular.size + 1)
    attributes += PhenotypeTable.taxa(genotypes.taxa.names)
    for (component in singular.indices) {
        val factor = scale * singular[component]
        attributes += PhenotypeTable.Attribute(
            "PC${component + 1}",
            FloatArray(numberOfTaxa) { (factor * eigen[it, component]).toFloat() },
        )
    }
    return PhenotypeTable(attributes)
}

/**
 * Runs [body] once per site with [column] holding that site's column of `A`: the reference
 * probabilities centered on their mean, divided by their standard deviation when the correlation
 * matrix was asked for, and scaled by `1 / sqrt(nTaxa - 1)`.
 *
 * @see net.maizegenetics.stats.PCA.PrinComp.centerCols
 * @see net.maizegenetics.stats.PCA.PrinComp.scaleCenteredMatrix
 */
private suspend inline fun forEachColumn(
    genotypes: GenotypeTable,
    probabilities: ReferenceProbability,
    useCovariance: Boolean,
    column: DoubleArray,
    progress: ProgressReporter,
    message: String,
    body: (site: Int) -> Unit,
) {
    val numberOfTaxa = column.size
    val multiplier = 1.0 / sqrt((numberOfTaxa - 1).toDouble())
    var done = 0
    genotypes.siteBlocks().forEach { block ->
        for (site in block.sites) {
            probabilities.column(block, site, column)
            var sum = 0.0
            for (taxon in 0 until numberOfTaxa) sum += column[taxon]
            val mean = sum / numberOfTaxa
            for (taxon in 0 until numberOfTaxa) column[taxon] -= mean
            if (!useCovariance) {
                var sumOfSquares = 0.0
                for (taxon in 0 until numberOfTaxa) sumOfSquares += column[taxon] * column[taxon]
                val deviation = sqrt(sumOfSquares / (numberOfTaxa - 1))
                for (taxon in 0 until numberOfTaxa) column[taxon] /= deviation
            }
            for (taxon in 0 until numberOfTaxa) column[taxon] *= multiplier
            body(site)
        }
        done += block.sites.count()
        progress.report(done, genotypes.numberOfSites, message)
    }
}

private const val GRAM_MESSAGE = "PCA: taxa covariance"
private const val EIGENVECTOR_MESSAGE = "PCA: eigenvectors"
private const val EIGENVALUE_TITLE = "Proportion of Variance Explained"
private const val EIGENVECTOR_TITLE = "Eigenvectors"
private const val TRAIT = "Trait"
private val EIGENVALUE_COLUMNS = listOf("PC", "eigenvalue", "proportion of total", "cumulative proportion")
