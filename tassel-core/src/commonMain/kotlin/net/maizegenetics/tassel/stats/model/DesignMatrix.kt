package net.maizegenetics.tassel.stats.model

import net.maizegenetics.tassel.math.DenseMatrix
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable

/** One column of a [DesignMatrix]: its [name] and the index in [ModelSpec.terms] of the term it belongs to. */
data class DesignColumn(val name: String, val termIndex: Int)

/**
 * The response and fixed-effect design of a [ModelSpec] over some observations of a [PhenotypeTable], with a map from
 * terms to columns so that sums of squares and tests can be computed per term.
 *
 * Columns follow the term order. A covariate gives one column and the intercept one; a factor gives `k - 1` columns
 * for the `k` levels present among [observations] (see [Contrast]), so a level emptied by row selection does not
 * create an all-zero column. Rows are the observations in the order given. Nothing here is estimated: [rank] and the
 * aliased-column report come from the same decomposition a fit uses, and a collinear design reports instead of
 * throwing (ADR-0025).
 */
class DesignMatrix internal constructor(
    /** The model this design was built for. */
    val spec: ModelSpec,
    /** The phenotype row of each design row. */
    val observations: IntArray,
    /** `observations.size x columns`. */
    val matrix: DenseMatrix,
    /** The response, one value per row. */
    val response: DoubleArray,
    /** Name and term of each column. */
    val columns: List<DesignColumn>,
    /**
     * The relative tolerance below which a column counts as a combination of the columns before it, or null for the
     * default of [net.maizegenetics.tassel.math.HouseholderQr.of]. Phenotype values are single precision (ADR-0023), so
     * a dependency among them holds only to about `1e-7`.
     */
    val rankTolerance: Double? = null,
) {
    init {
        require(matrix.columns == columns.size) { "${columns.size} column descriptions for ${matrix.columns} columns" }
        require(matrix.rows == response.size && matrix.rows == observations.size) {
            "matrix has ${matrix.rows} rows, response ${response.size}, observations ${observations.size}"
        }
    }

    /** Number of rows. */
    val rows: Int get() = matrix.rows

    /** Number of columns, before aliasing. */
    val columnCount: Int get() = matrix.columns

    /** The column indices that belong to the term at [termIndex], ascending (empty for a one-level factor). */
    fun termColumns(termIndex: Int): IntArray {
        require(termIndex in spec.terms.indices) { "term $termIndex is outside ${spec.terms.indices}" }
        return columns.indices.filter { columns[it].termIndex == termIndex }.toIntArray()
    }

    internal val decomposition: AliasAwareQr by lazy { AliasAwareQr.of(matrix, rankTolerance) }

    /** Numerical rank of the design. */
    val rank: Int get() = decomposition.rank

    /** Indices of the aliased columns: numerically combinations of the columns before them. Empty when full rank. */
    val aliasedColumns: IntArray get() = decomposition.aliased

    /** The names of [aliasedColumns]. */
    val aliasedColumnNames: List<String> get() = aliasedColumns.map { columns[it].name }

    /** Labels of the terms with columns, every one of which is aliased. */
    val aliasedTerms: List<String>
        get() = spec.terms.indices.filter { term ->
            val own = termColumns(term)
            own.isNotEmpty() && own.all { it in aliasedColumns }
        }.map { spec.terms[it].label }

    /** The design restricted to the terms at [termIndices] (indices into [ModelSpec.terms]), in that order. */
    fun subset(termIndices: List<Int>): DesignMatrix {
        val chosen = termIndices.flatMap { term -> termColumns(term).toList().map { it to term } }
        val newColumns = chosen.map { (column, term) -> DesignColumn(columns[column].name, termIndices.indexOf(term)) }
        val sub = DenseMatrix(rows, chosen.size)
        chosen.forEachIndexed { target, (source, _) ->
            matrix.data.copyInto(sub.data, target * rows, source * rows, source * rows + rows)
        }
        return DesignMatrix(spec.withTerms(termIndices), observations, sub, response, newColumns, rankTolerance)
    }

    /** The same design with another matrix and response of the same shape, for a transformed (whitened) model. */
    internal fun withData(newMatrix: DenseMatrix, newResponse: DoubleArray): DesignMatrix =
        DesignMatrix(spec, observations, newMatrix, newResponse, columns, rankTolerance)

    companion object {
        /**
         * The observations of [phenotype] with no missing value in the response or any covariate or factor of [spec]
         * (ADR-0023: NaN and code -1 are missing), ascending.
         */
        fun completeObservations(spec: ModelSpec, phenotype: PhenotypeTable): IntArray {
            val used = used(spec, phenotype)
            return (0 until phenotype.numberOfObservations).filter { row ->
                used.none { it.isMissing(row) }
            }.toIntArray()
        }

        /**
         * Builds the design of [spec] for [observations], rows of [phenotype] with no missing value in any attribute
         * the model uses (use [completeObservations] to pick them). Marker and random-effect terms are not
         * supported yet.
         *
         * @param rankTolerance see [DesignMatrix.rankTolerance]
         * @throws IllegalArgumentException for an unknown or wrongly typed attribute, an observation outside the
         *   table, or a missing value in a selected row
         */
        fun build(
            spec: ModelSpec,
            phenotype: PhenotypeTable,
            observations: IntArray = completeObservations(spec, phenotype),
            rankTolerance: Double? = null,
        ): DesignMatrix {
            observations.forEach {
                require(it in 0 until phenotype.numberOfObservations) {
                    "Observation $it is outside 0 until ${phenotype.numberOfObservations}"
                }
            }
            val n = observations.size
            val responseAttribute = numeric(phenotype, spec.response, "response")
            val response = DoubleArray(n) { checkedValue(responseAttribute, observations[it]) }

            val parts = mutableListOf<Pair<String, DoubleArray>>()
            val columns = mutableListOf<DesignColumn>()
            spec.terms.forEachIndexed { termIndex, term ->
                val built: List<Pair<String, DoubleArray>> = when (term) {
                    ModelTerm.Intercept -> listOf(term.label to DoubleArray(n) { 1.0 })
                    is ModelTerm.Covariate -> {
                        val attribute = numeric(phenotype, term.name, "covariate")
                        listOf(term.label to DoubleArray(n) { checkedValue(attribute, observations[it]) })
                    }
                    is ModelTerm.Factor -> factorColumns(term, factor(phenotype, term.name), observations)
                    is ModelTerm.MarkerTerm ->
                        throw UnsupportedOperationException("Marker terms are not supported yet (tassel-6#189)")
                    is ModelTerm.RandomEffect ->
                        throw UnsupportedOperationException("Random effects are not fixed-effect columns (ADR-0028)")
                }
                for ((name, values) in built) {
                    parts += name to values
                    columns += DesignColumn(name, termIndex)
                }
            }
            val matrix = DenseMatrix(n, parts.size)
            parts.forEachIndexed { column, (_, values) -> values.copyInto(matrix.data, column * n) }
            return DesignMatrix(spec, observations.copyOf(), matrix, response, columns, rankTolerance)
        }

        private fun factorColumns(
            term: ModelTerm.Factor,
            attribute: FactorAttribute,
            observations: IntArray,
        ): List<Pair<String, DoubleArray>> {
            val codes = IntArray(observations.size) {
                attribute.code(observations[it]).also { code ->
                    require(code >= 0) { "Factor '${term.name}' is missing at observation ${observations[it]}" }
                }
            }
            val present = codes.distinct().sorted()
            val last = present.lastOrNull() ?: return emptyList()
            return present.dropLast(1).map { level ->
                "${term.name}=${attribute.levels[level]}" to DoubleArray(codes.size) { row ->
                    val own = if (codes[row] == level) 1.0 else 0.0
                    when (term.contrast) {
                        Contrast.treatment -> own
                        Contrast.sum -> own - if (codes[row] == last) 1.0 else 0.0
                    }
                }
            }
        }

        private fun checkedValue(attribute: NumericAttribute, observation: Int): Double {
            require(!attribute.isMissing(observation)) {
                "'${attribute.name}' is missing at observation $observation"
            }
            return attribute.value(observation)
        }

        private fun numeric(phenotype: PhenotypeTable, name: String, role: String): NumericAttribute {
            require(phenotype.indexOf(name) >= 0) { "The $role '$name' is not an attribute of '${phenotype.name}'" }
            val attribute = phenotype.attribute(name)
            require(attribute is NumericAttribute) { "The $role '$name' must be numeric, not ${attribute.type}" }
            return attribute
        }

        private fun factor(phenotype: PhenotypeTable, name: String): FactorAttribute {
            require(phenotype.indexOf(name) >= 0) { "The factor '$name' is not an attribute of '${phenotype.name}'" }
            val attribute = phenotype.attribute(name)
            require(attribute is FactorAttribute) { "'$name' must be a factor, not ${attribute.type}" }
            return attribute
        }

        private fun used(spec: ModelSpec, phenotype: PhenotypeTable) = buildList {
            add(numeric(phenotype, spec.response, "response"))
            for (term in spec.terms) {
                when (term) {
                    is ModelTerm.Covariate -> add(numeric(phenotype, term.name, "covariate"))
                    is ModelTerm.Factor -> add(factor(phenotype, term.name))
                    else -> Unit
                }
            }
        }
    }
}
