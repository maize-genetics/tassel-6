package net.maizegenetics.tassel.stats.model

/** One row of an [AnovaTable]: a term's degrees of freedom, sum of squares, mean square, and F test. */
class AnovaRow(val term: String, val df: Int, val sumOfSquares: Double, val test: TestResult) {
    /** Mean square, `NaN` when the term has no degrees of freedom. */
    val meanSquare: Double get() = if (df > 0) sumOfSquares / df else Double.NaN

    val f: Double get() = test.statistic
    val p: Double get() = test.p
    val log10P: Double get() = test.log10P
}

/** Per-term sums of squares and F tests against the full model's residual mean square, plus the residual row. */
class AnovaTable(val rows: List<AnovaRow>, val residualDf: Int, val residualSumOfSquares: Double) {
    /** The row of [term], or null. */
    fun row(term: String): AnovaRow? = rows.firstOrNull { it.term == term }
}

/**
 * Sequential (type I) and marginal (drop-one) sums of squares, from nested fits of the same design. A term that is
 * aliased, or has no columns, gets zero degrees of freedom and an inestimable test. The intercept has no row.
 *
 * Sequential sums of squares add the terms in model order: row `i` is `RSS(terms before i) - RSS(through i)` on
 * `rank` difference degrees of freedom, as R's `anova()`. Marginal sums of squares drop one term at a time from the
 * full model: `RSS(without i) - RSS(full)`, as R's `drop1()`.
 *
 * @see net.maizegenetics.stats.linearmodels.LinearModelUtils
 */
object Anova {
    /** Sequential (type I) table for [design]. */
    fun sequential(design: DesignMatrix, refine: Boolean = false): AnovaTable {
        val terms = design.spec.terms
        val full = LinearModelFit.of(design, refine)
        var previousRank = 0
        var previousRss = design.response.fold(0.0) { sum, v -> sum + v * v }
        val rows = mutableListOf<AnovaRow>()
        for (index in terms.indices) {
            val fit = if (index ==
                terms.lastIndex
            ) {
                full
            } else {
                LinearModelFit.of(design.subset((0..index).toList()), refine)
            }
            if (terms[index] != ModelTerm.Intercept) {
                val df = fit.rank - previousRank
                val ss = (previousRss - fit.residualSumOfSquares).coerceAtLeast(0.0)
                rows += row(terms[index].label, df, ss, full)
            }
            previousRank = fit.rank
            previousRss = fit.residualSumOfSquares
        }
        return AnovaTable(rows, full.residualDf, full.residualSumOfSquares)
    }

    /** Marginal (drop-one) table for [design]: every term but the intercept, each tested after all the others. */
    fun marginal(design: DesignMatrix, refine: Boolean = false): AnovaTable {
        val terms = design.spec.terms
        val full = LinearModelFit.of(design, refine)
        val rows = terms.indices.filter { terms[it] != ModelTerm.Intercept }.map { index ->
            val reduced = LinearModelFit.of(design.subset(terms.indices.filter { it != index }), refine)
            row(
                terms[index].label,
                full.rank - reduced.rank,
                (reduced.residualSumOfSquares - full.residualSumOfSquares).coerceAtLeast(0.0),
                full,
            )
        }
        return AnovaTable(rows, full.residualDf, full.residualSumOfSquares)
    }

    private fun row(label: String, df: Int, ss: Double, full: LinearModelFit): AnovaRow {
        val f = if (df > 0 && full.residualDf > 0) (ss / df) / full.residualVariance else Double.NaN
        return AnovaRow(label, df, ss, TestResult.fromF(f, df, full.residualDf))
    }
}
