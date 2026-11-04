package net.maizegenetics.tassel.ops

import net.maizegenetics.tassel.analysis.DEFAULT_COMPONENTS
import net.maizegenetics.tassel.analysis.PrincipalComponentsResult
import net.maizegenetics.tassel.analysis.principalComponents
import net.maizegenetics.tassel.genotype.GenotypeTable

/**
 * How many principal components to return. The constant names are TASSEL 5's, because they are
 * what `-limitBy` accepts on the command line.
 *
 * @see net.maizegenetics.analysis.data.PrincipalComponentsPlugin.PCA_LIMIT
 */
@Suppress("ktlint:standard:enum-entry-name-case")
enum class PcaLimit { number_of_components, min_eigenvalue, total_variance }

/**
 * Principal components of a genotype table.
 *
 * Every TASSEL 5 parameter is declared, so the command line and the legacy translator keep working,
 * but only `limitBy = number_of_components` is ported; see the `run` block for why the other two
 * have no behavior to match.
 *
 * @see net.maizegenetics.analysis.data.PrincipalComponentsPlugin
 */
val pcaOperation = operation<GenotypeTable, PrincipalComponentsResult>(
    id = "pca",
    legacyPlugin = "PrincipalComponentsPlugin",
    help = "Principal components of the taxa in a genotype table",
) {
    val covariance by boolean(
        default = true,
        legacy = "covariance",
        help = "If the box is checked, then the analysis will do an eigenvalue decomposition of the covariance " +
            "matrix. If the box is unchecked, it will use a correlation matrix. Using the covariance matrix is " +
            "recommended for genotypes while the correlation matrix is often used for phenotypes.",
    )
    val limitBy by enum(
        default = PcaLimit.number_of_components,
        legacy = "limitBy",
        help = "This parameter determines the type of value that will be used to limit the number of principal " +
            "components (axes) returned.",
    )
    val numberOfComponents by int(
        default = DEFAULT_COMPONENTS,
        legacy = "ncomponents",
        help = "The analysis will return this many principal components up to the number of taxa.",
    )
    val minimumEigenvalue by double(
        default = 0.0,
        legacy = "minEigenval",
        help = "All principal components with an eigenvalue greater than or equal to this value will be returned.",
    )
    val totalVariance by double(
        default = 0.5,
        range = 0.0..1.0,
        legacy = "totalVar",
        help = "The first principal components that together explain this proportion of the total variance will " +
            "be returned.",
    )
    val reportEigenvalues by boolean(
        default = true,
        legacy = "reportEigenvalues",
        help = "Returns a list of eigenvalues sorted high to low.",
    )
    val reportEigenvectors by boolean(
        default = true,
        legacy = "reportEigenvectors",
        help = "Returns the eigenvectors calculated from a Singular Value Decomposition of the data. The resulting " +
            "table can be quite large if the number of variants and taxa are big.",
    )
    run { input, progress ->
        val limit = limitBy()
        // TASSEL 5 picks the count for the other two with Arrays.binarySearch over an array its SVD
        // does not sort, so neither has defined behavior to port. #19 scoped them out.
        require(limit == PcaLimit.number_of_components) {
            "pca: limitBy = ${limit.name} is not ported; minEigenval (${minimumEigenvalue()}) and " +
                "totalVar (${totalVariance()}) are declared for the legacy translator only. " +
                "Use limitBy = ${PcaLimit.number_of_components.name}."
        }
        principalComponents(
            genotypes = input,
            numberOfComponents = numberOfComponents(),
            useCovariance = covariance(),
            reportEigenvalues = reportEigenvalues(),
            reportEigenvectors = reportEigenvectors(),
            progress = progress,
        )
    }
}
