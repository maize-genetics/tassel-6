package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.table.DoubleFormat
import net.maizegenetics.tassel.table.doubleToString

/** How a `Double` cell is written. */
internal enum class NumberForm {
    /** TASSEL 5's pipeline `-export`: [DoubleFormat], five digits. */
    Tassel,

    /** TASSEL 5's `saveToFile`, and the modern reports: Java's `Double.toString`, every digit. */
    Full,
    ;

    fun format(value: Double): String = when (this) {
        Tassel -> DoubleFormat.format(value)
        Full -> doubleToString(value)
    }
}

/** Which TASSEL 5 GLM class's layout to write. */
internal enum class GlmLayout {
    /** `DiscreteSitesFELM`: marker classes, with optional add/dom columns. */
    Genotype,

    /** `ReferenceProbabilityFELM`: one dosage, no add/dom tests. */
    ReferenceProbability,
}

private val STATS_START = listOf("Trait", "Marker", "Chr", "Pos", "marker_F", "p")

/**
 * Writes TASSEL 5's GLM statistics report: a header, then one row per (trait, marker), trait by trait, with `\n` line
 * endings including after the last row, as `TableReportUtils` writes it.
 *
 * For [GlmLayout.Genotype] the row's main columns come from the marker's omnibus record, `add_F` and `add_p` from its
 * additive record, and `dom_F` and `dom_p` from its dominance record; a test the table lacks is `NaN`. Columns:
 * `Trait, Marker, Chr, Pos, marker_F, p, [perm_p,] marker_Rsq, add_F, add_p, dom_F, dom_p, marker_df, marker_MS,
 * error_df, error_MS, model_df, model_MS, minorObs[, addEffect, domEffect]`. For [GlmLayout.ReferenceProbability]:
 * `Trait, Marker, Chr, Pos, marker_F, p, [perm_p,] marker_Rsq, marker_df, marker_MS, error_df, error_MS, model_df,
 * model_MS`.
 *
 * @param permutations Write `perm_p` after `p`.
 * @param appendAddDom Write `addEffect` and `domEffect` last (genotype layout only).
 * @see net.maizegenetics.analysis.association.AbstractFixedEffectLM
 * @see net.maizegenetics.analysis.association.DiscreteSitesFELM
 * @see net.maizegenetics.analysis.association.ReferenceProbabilityFELM
 */
internal fun AssociationTable.writeGlmStatistics(
    out: Appendable,
    layout: GlmLayout = GlmLayout.Genotype,
    form: NumberForm = NumberForm.Tassel,
    permutations: Boolean = false,
    appendAddDom: Boolean = false,
) {
    val genotype = layout == GlmLayout.Genotype
    val header = buildList {
        addAll(STATS_START)
        if (permutations) add("perm_p")
        add("marker_Rsq")
        if (genotype) addAll(listOf("add_F", "add_p", "dom_F", "dom_p"))
        addAll(listOf("marker_df", "marker_MS", "error_df", "error_MS", "model_df", "model_MS"))
        if (genotype) add("minorObs")
        if (genotype && appendAddDom) addAll(listOf("addEffect", "domEffect"))
    }
    out.appendLine(header, '\t')
    val order = recordOrder()
    var i = 0
    while (i < order.size) {
        val main = get(order[i])
        var additive: AssociationRecord? = null
        var dominance: AssociationRecord? = null
        var j = i + 1
        while (j < order.size) {
            val next = get(order[j])
            if (next.trait != main.trait || next.markerIndex != main.markerIndex) break
            when (next.test) {
                AssociationTestKind.Additive -> additive = next
                AssociationTestKind.Dominance -> dominance = next
                else -> {}
            }
            j++
        }
        val cells = buildList {
            add(main.trait)
            add(main.marker)
            add(main.chromosome)
            add(main.position.toString())
            add(form.format(main.statistic))
            add(form.format(main.p))
            if (permutations) add(form.format(main.permutationP))
            add(form.format(main.rSquared))
            if (genotype) {
                add(form.format(additive?.statistic ?: Double.NaN))
                add(form.format(additive?.p ?: Double.NaN))
                add(form.format(dominance?.statistic ?: Double.NaN))
                add(form.format(dominance?.p ?: Double.NaN))
            }
            add(form.format(main.numeratorDf))
            add(form.format(main.markerMeanSquare))
            add(form.format(main.denominatorDf))
            add(form.format(main.errorMeanSquare))
            add(form.format(main.modelDf))
            add(form.format(main.modelMeanSquare))
            if (genotype) add(main.minorObs.toString())
            if (genotype && appendAddDom) {
                add(form.format(additive?.effect ?: Double.NaN))
                add(form.format(dominance?.effect ?: Double.NaN))
            }
        }
        out.appendLine(cells, '\t')
        i = j
    }
}

/**
 * Writes TASSEL 5's GLM genotype-effects report, one row per effect, trait by trait: `Trait, Marker, Chr, Pos, Obs,
 * Allele, Estimate`, or `Trait, Marker, Chr, Pos, Estimate` for [GlmLayout.ReferenceProbability].
 *
 * @see net.maizegenetics.analysis.association.AbstractFixedEffectLM
 */
internal fun AssociationTable.writeGlmEffects(
    out: Appendable,
    layout: GlmLayout = GlmLayout.Genotype,
    form: NumberForm = NumberForm.Tassel,
) {
    val genotype = layout == GlmLayout.Genotype
    out.appendLine(
        buildList {
            addAll(listOf("Trait", "Marker", "Chr", "Pos"))
            if (genotype) addAll(listOf("Obs", "Allele"))
            add("Estimate")
        },
        '\t',
    )
    for (row in effectOrder()) {
        val e = effect(row)
        out.appendLine(
            buildList {
                add(e.trait)
                add(e.marker)
                add(e.chromosome)
                add(e.position.toString())
                if (genotype) {
                    add(e.observations.toString())
                    add(e.allele)
                }
                add(form.format(e.estimate))
            },
            '\t',
        )
    }
}

/** First line of a modern statistics report. */
internal const val ASSOCIATION_SCHEMA = "#schema=association/1"

/** First line of a modern effects report. */
internal const val ASSOCIATION_EFFECTS_SCHEMA = "#schema=association-effects/1"

internal val MODERN_COLUMNS = listOf(
    "trait", "marker_index", "marker", "chr", "pos", "model", "test", "encoding", "effect_allele",
    "comparison_allele", "observations", "taxa", "allele_count", "allele_frequency", "minor_allele_count",
    "minor_obs", "effect", "standard_error", "ci_lower", "ci_upper", "statistic", "df1", "df2", "p", "log10_p",
    "r_squared", "marker_ms", "error_ms", "model_df", "model_ms", "perm_p", "genetic_variance",
    "residual_variance", "heritability", "variance_estimation", "status", "skip_reason", "message",
)

internal val MODERN_EFFECT_COLUMNS = listOf(
    "trait", "marker_index", "marker", "chr", "pos", "allele", "observations", "estimate", "standard_error",
)

/**
 * Writes the modern report: [ASSOCIATION_SCHEMA], the provenance as `##key=value` lines, a header, and one tab
 * separated row per record with every field of [AssociationRecord] in [MODERN_COLUMNS] order. Numbers are written with
 * every digit, a count or double that does not apply is empty, and `log10_p` is explicit. Rows keep append order.
 */
internal fun AssociationTable.writeModernStatistics(out: Appendable) {
    out.append(ASSOCIATION_SCHEMA).append('\n')
    provenance.headerLines().forEach { out.append(it).append('\n') }
    out.appendLine(MODERN_COLUMNS, '\t')
    for (row in 0 until size) {
        val r = get(row)
        out.appendLine(
            listOf(
                r.trait, r.markerIndex.toString(), r.marker, r.chromosome, r.position.toString(),
                r.model.name.lowercase(), r.test.name.lowercase(), r.encoding.name, r.effectAllele,
                r.comparisonAllele, count(r.observations), count(r.taxa), count(r.alleleCount),
                number(r.alleleFrequency), count(r.minorAlleleCount), count(r.minorObs), number(r.effect),
                number(r.standardError), number(r.ciLower), number(r.ciUpper), number(r.statistic),
                number(r.numeratorDf), number(r.denominatorDf), number(r.p), number(r.log10P),
                number(r.rSquared), number(r.markerMeanSquare), number(r.errorMeanSquare), number(r.modelDf),
                number(r.modelMeanSquare), number(r.permutationP), number(r.geneticVariance),
                number(r.residualVariance), number(r.heritability), r.varianceEstimation.name.lowercase(),
                r.status.name.lowercase(), r.skipReason?.name ?: "", clean(r.message),
            ),
            '\t',
        )
    }
}

/** Writes the effects in the modern form: [ASSOCIATION_EFFECTS_SCHEMA], provenance, a header, and every field. */
internal fun AssociationTable.writeModernEffects(out: Appendable) {
    out.append(ASSOCIATION_EFFECTS_SCHEMA).append('\n')
    provenance.headerLines().forEach { out.append(it).append('\n') }
    out.appendLine(MODERN_EFFECT_COLUMNS, '\t')
    for (row in 0 until effectCount) {
        val e = effect(row)
        out.appendLine(
            listOf(
                e.trait, e.markerIndex.toString(), e.marker, e.chromosome, e.position.toString(), e.allele,
                count(e.observations), number(e.estimate), number(e.standardError),
            ),
            '\t',
        )
    }
}

private fun count(value: Int): String = if (value < 0) "" else value.toString()

private fun number(value: Double): String = if (value.isNaN()) "" else doubleToString(value)

private fun clean(text: String): String = text.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ')

private fun Appendable.appendLine(cells: List<String>, delimiter: Char) {
    cells.forEachIndexed { index, cell ->
        if (index != 0) append(delimiter)
        append(cell)
    }
    append('\n')
}
