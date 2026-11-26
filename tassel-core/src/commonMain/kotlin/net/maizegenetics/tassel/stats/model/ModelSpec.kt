package net.maizegenetics.tassel.stats.model

/**
 * How a factor's levels become design-matrix columns. A factor with `k` levels present in the fitted rows gets
 * `k - 1` columns; levels are in the order of [net.maizegenetics.tassel.phenotype.FactorAttribute.levels] (sorted,
 * ADR-0023), and levels no fitted observation uses are left out, as R's `lm` leaves them out.
 *
 * @see net.maizegenetics.stats.linearmodels.FactorModelEffect
 */
@Suppress("ktlint:standard:enum-entry-name-case")
enum class Contrast {
    /**
     * TASSEL 5's convention: the LAST level is the reference and its effect is zero. Column `i` is the indicator of
     * level `i`, for `i < k - 1`. (R's default drops the first level; the sums of squares do not depend on the choice,
     * and R can be asked for the same coding with `contr.treatment(levels, base = k)`.)
     */
    treatment,

    /** Sum-to-zero coding: column `i` is the indicator of level `i` minus the indicator of the last level (R's `contr.sum`). */
    sum,
}

/** One term of a [ModelSpec]. The hierarchy is sealed, but open to new subclasses later (ADR-0025). */
sealed class ModelTerm {
    /** The term's name: unique within a spec, and the label of its row in a table of sums of squares. */
    abstract val label: String

    /** The intercept: one column of ones. */
    data object Intercept : ModelTerm() {
        override val label: String get() = "(Intercept)"
    }

    /** A numeric attribute entered as it is, one column. Missing values are not allowed in the fitted rows. */
    data class Covariate(val name: String) : ModelTerm() {
        init {
            require(name.isNotEmpty()) { "A covariate needs a name" }
        }

        override val label: String get() = name
    }

    /** A categorical attribute, coded by [contrast] into one column per level but one. */
    data class Factor(val name: String, val contrast: Contrast = Contrast.treatment) : ModelTerm() {
        init {
            require(name.isNotEmpty()) { "A factor needs a name" }
        }

        override val label: String get() = name
    }

    /** Placeholder for a genotype-marker term (tassel-6#189 defines the encodings). A design cannot be built from it yet. */
    data class MarkerTerm(val name: String = "Marker") : ModelTerm() {
        override val label: String get() = name
    }

    /** Placeholder for a random effect with a kinship covariance (ADR-0028, tassel-6#196). Not a fixed-effect column. */
    data class RandomEffect(val name: String) : ModelTerm() {
        override val label: String get() = name
    }
}

/**
 * A linear model: a [response] attribute and an ordered list of [terms]. The order matters, because sequential sums
 * of squares add terms in it and a collinear term is aliased with the terms before it (R's convention).
 *
 * A formula string such as `Height ~ Marker + PC1` is an optional front end for later (tassel-6#208).
 */
data class ModelSpec(val response: String, val terms: List<ModelTerm>) {
    init {
        require(response.isNotEmpty()) { "A model needs a response" }
        require(terms.isNotEmpty()) { "A model needs at least one term" }
        val labels = terms.map { it.label }
        require(labels.size == labels.distinct().size) { "Term labels repeat in $labels" }
    }

    constructor(response: String, vararg terms: ModelTerm) : this(response, terms.toList())

    /** The spec with only the terms at [indices], in that order. */
    fun withTerms(indices: List<Int>): ModelSpec = ModelSpec(response, indices.map { terms[it] })
}
