package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.math.DenseMatrix
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.phenotype.TaxaAttribute
import net.maizegenetics.tassel.stats.model.Contrast
import net.maizegenetics.tassel.stats.model.DesignColumn
import net.maizegenetics.tassel.stats.model.DesignMatrix
import net.maizegenetics.tassel.stats.model.LinearModelFit
import net.maizegenetics.tassel.stats.model.ModelSpec
import net.maizegenetics.tassel.stats.model.ModelTerm

/**
 * Settings of [phenotypeLm].
 *
 * @property atCovariateMeans Evaluate each BLUE at the mean of every covariate over the trait's observations. Off by
 * default: TASSEL 5 computes BLUEs at covariate value 0, so the covariates contribute nothing. New in TASSEL 6.
 */
internal data class PhenotypeLmOptions(val atCovariateMeans: Boolean = false)

/**
 * One row of TASSEL 5's phenotype-only report: the test of the taxon factor for a trait, as the last effect of the
 * model after the mean, the factors, and the covariates.
 *
 * With unreplicated taxa the model fits every observation exactly, so [errorDf] is 0 and [f], [p], and [errorMs] are
 * `NaN`. TASSEL 5 prints the sign of its sweep's rounding noise there (`F = -0`, `errorMS = Infinity`).
 */
internal class PhenotypeAnovaRow(
    val trait: String,
    val f: Double,
    val p: Double,
    val taxaDf: Double,
    val taxaMs: Double,
    val errorDf: Double,
    val errorMs: Double,
    val modelDf: Double,
    val modelMs: Double,
)

/** The result of [phenotypeLm]: [blues] as a phenotype of one row per taxon, and one [anova] row per trait. */
internal class PhenotypeLmResult(val blues: PhenotypeTable, val anova: List<PhenotypeAnovaRow>)

/** The columns of TASSEL 5's phenotype-only report. */
internal val PHENOTYPE_ANOVA_COLUMNS =
    listOf("Trait", "F", "p", "taxaDF", "taxaMS", "errorDF", "errorMS", "modelDF", "modelMS")

/**
 * Writes the report as TASSEL 5's `PhenotypeLM` report: a header and one row per trait, with `\n` line endings. Every
 * number is a `Double`, so degrees of freedom are `29.0` in [NumberForm.Full] and `29` in [NumberForm.Tassel].
 */
internal fun List<PhenotypeAnovaRow>.writePhenotypeAnova(out: Appendable, form: NumberForm = NumberForm.Tassel) {
    out.append(PHENOTYPE_ANOVA_COLUMNS.joinToString("\t")).append('\n')
    for (row in this) {
        val cells = listOf(
            row.trait,
            form.format(row.f),
            form.format(row.p),
            form.format(row.taxaDf),
            form.format(row.taxaMs),
            form.format(row.errorDf),
            form.format(row.errorMs),
            form.format(row.modelDf),
            form.format(row.modelMs),
        )
        out.append(cells.joinToString("\t")).append('\n')
    }
}

/**
 * `glm` on a phenotype without genotypes (TASSEL 5's `phenoOnly`): for each trait, fits the fixed-effect model
 * `trait ~ 1 + factors + covariates + taxon`, tests the taxon effect, and reports the best linear unbiased estimate
 * (BLUE) of each taxon. Taxa may be replicated.
 *
 * Per trait, the observations are those where the trait, every covariate, and every factor are present. The model is
 * the mean, the factors, the covariates, and then the taxa as a factor whose last taxon (in order of first appearance)
 * is the zero reference. The taxon test is the sequential sum of squares of that last effect against the residual, so
 * unreplicated taxa give error df 0 and `NaN` (see [PhenotypeAnovaRow]). `modelDf` and `modelMs` exclude the mean.
 *
 * A BLUE is the taxon's effect plus the overall mean, where the overall mean is the intercept plus, for each factor, the
 * sum of its level effects over the number of its levels (the average level). Covariates contribute nothing, as in
 * TASSEL 5 ([PhenotypeLmOptions.atCovariateMeans] changes that). A taxon with no observation of a trait has a missing
 * BLUE; a taxon aliased with the other effects is estimated as zero effect, as TASSEL 5's sweep does. BLUEs are stored
 * as `Float`s (ADR-0023), and so are TASSEL 5's. Rows follow the taxa in order of first appearance.
 *
 * The sums of squares come from `Q'y` of one Householder QR per trait and nothing is inverted (ADR-0025). The BLUE
 * phenotype has a taxa attribute and one `data` attribute per trait, so it can be joined to genotypes and analyzed by
 * `glm`.
 *
 * @see net.maizegenetics.analysis.association.PhenotypeLM
 * @see net.maizegenetics.analysis.association.FixedEffectLMPlugin
 */
internal fun phenotypeLm(
    phenotype: PhenotypeTable,
    spec: AssociationSpec = AssociationSpec(),
    options: PhenotypeLmOptions = PhenotypeLmOptions(),
): PhenotypeLmResult {
    val roles = resolveRoles(phenotype, spec)
    val taxonNames = phenotype.taxa.map { it.name }
    val taxonIndex = HashMap<String, Int>(taxonNames.size * 2).also { map ->
        taxonNames.forEachIndexed { index, name -> map[name] = index }
    }
    val covariates = roles.covariates.map { phenotype.attribute(it) as NumericAttribute }
    val factors = roles.factors.map { phenotype.attribute(it) as FactorAttribute }

    val blueAttributes = ArrayList<NumericAttribute>(roles.traits.size)
    val rows = ArrayList<PhenotypeAnovaRow>(roles.traits.size)
    for (name in roles.traits) {
        val trait = phenotype.attribute(name) as NumericAttribute
        val observations = (0 until phenotype.numberOfObservations).filter { row ->
            !trait.isMissing(row) && covariates.none { it.isMissing(row) } && factors.none { it.isMissing(row) }
        }.toIntArray()
        val fit = fitTrait(phenotype, name, roles, observations, options)
        val values = FloatArray(taxonNames.size) { Float.NaN }
        for ((taxon, blue) in fit.blues) values[taxonIndex.getValue(taxon)] = blue.toFloat()
        blueAttributes += NumericAttribute(name, AttributeType.data, values)
        rows += fit.row
    }
    val taxa = TaxaAttribute(phenotype.taxaAttribute.name, taxonNames)
    return PhenotypeLmResult(
        PhenotypeTable(listOf(taxa) + blueAttributes, "BLUEs_${phenotype.name}"),
        rows,
    )
}

private class TraitFit(val row: PhenotypeAnovaRow, val blues: Map<String, Double>)

private fun fitTrait(
    phenotype: PhenotypeTable,
    trait: String,
    roles: AssociationRoles,
    observations: IntArray,
    options: PhenotypeLmOptions,
): TraitFit {
    val baseSpec = ModelSpec(
        trait,
        buildList {
            add(ModelTerm.Intercept)
            roles.factors.forEach { add(ModelTerm.Factor(it, Contrast.treatment)) }
            roles.covariates.forEach { add(ModelTerm.Covariate(it)) }
        },
    )
    if (observations.isEmpty()) {
        val none = Double.NaN
        return TraitFit(PhenotypeAnovaRow(trait, none, none, 0.0, none, 0.0, none, 0.0, none), emptyMap())
    }
    val base = DesignMatrix.build(baseSpec, phenotype, observations, GLM_RANK_TOLERANCE)

    // One indicator column per taxon but the last in order of first appearance, as a restricted factor (TASSEL 5's
    // `taxaEffect`); that last taxon is the zero reference.
    val taxonOfObservation = Array(observations.size) { phenotype.taxon(observations[it]) }
    val present = LinkedHashSet<String>().also { it.addAll(taxonOfObservation) }.toList()
    val n = observations.size
    val p0 = base.columnCount
    val taxonColumns = present.size - 1
    val matrix = DenseMatrix(n, p0 + maxOf(taxonColumns, 0))
    base.matrix.data.copyInto(matrix.data, 0, 0, n * p0)
    for (j in 0 until taxonColumns) {
        val level = present[j]
        for (row in 0 until n) if (taxonOfObservation[row] == level) matrix.data[(p0 + j) * n + row] = 1.0
    }
    val termIndex = baseSpec.terms.size
    val design = DesignMatrix(
        ModelSpec(trait, baseSpec.terms + ModelTerm.Factor(TAXON_TERM)),
        base.observations,
        matrix,
        base.response,
        base.columns + (0 until taxonColumns).map { DesignColumn("$TAXON_TERM=${present[it]}", termIndex) },
        GLM_RANK_TOLERANCE,
    )

    val sums = SequentialSums(design)
    val rank = sums.rank
    val taxaRank = rank - sums.rankThrough(p0)
    val taxaSs = sums.ss(rank - taxaRank, rank)
    val errorDf = n - rank
    val errorSs = sums.ss(rank, n)
    val modelDf = rank - 1
    val modelSs = sums.ss(1, rank)
    val test = Fstat.of(taxaSs, taxaRank, errorSs, errorDf)
    val row = PhenotypeAnovaRow(
        trait = trait,
        f = test.f,
        p = test.p,
        taxaDf = taxaRank.toDouble(),
        taxaMs = taxaSs / taxaRank,
        errorDf = errorDf.toDouble(),
        errorMs = errorSs / errorDf,
        modelDf = modelDf.toDouble(),
        modelMs = modelSs / modelDf,
    )

    // BLUE = taxon effect + overall mean; an aliased column is a zero effect.
    val beta = LinearModelFit.of(design).coefficients.map { if (it.isNaN()) 0.0 else it }
    var overall = beta[0]
    roles.factors.indices.forEach { f ->
        val columns = design.termColumns(1 + f)
        val levels = columns.size + 1
        overall += columns.sumOf { beta[it] } / levels
    }
    if (options.atCovariateMeans) {
        roles.covariates.forEachIndexed { c, name ->
            val column = design.termColumns(1 + roles.factors.size + c).singleOrNull() ?: return@forEachIndexed
            val attribute = phenotype.attribute(name) as NumericAttribute
            overall += beta[column] * (observations.sumOf { attribute.value(it) } / n)
        }
    }
    val blues = LinkedHashMap<String, Double>(present.size * 2)
    present.forEachIndexed { index, taxon ->
        blues[taxon] = if (index < taxonColumns) overall + beta[p0 + index] else overall
    }
    return TraitFit(row, blues)
}

private const val TAXON_TERM = "(taxon)"
