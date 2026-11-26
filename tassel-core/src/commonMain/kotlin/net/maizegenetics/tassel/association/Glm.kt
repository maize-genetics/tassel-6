package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.math.DenseMatrix
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.ops.report
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.stats.StudentTDistribution
import net.maizegenetics.tassel.stats.fTestLog10PValue
import net.maizegenetics.tassel.stats.fTestPValue
import net.maizegenetics.tassel.stats.model.DesignColumn
import net.maizegenetics.tassel.stats.model.DesignMatrix
import net.maizegenetics.tassel.stats.model.LinearModelFit
import net.maizegenetics.tassel.stats.model.ModelSpec
import net.maizegenetics.tassel.stats.model.ModelTerm
import kotlin.math.abs

/**
 * Settings of a GLM run beyond the model: marker encoding and filters, `appendAddDom`, `maxP`, and the site block size.
 *
 * @property appendAddDom TASSEL 5's `appendAddDom`: put the additive and dominance effects in the `effect` of the
 * additive and dominance records of a biallelic marker (written as `addEffect` and `domEffect` by the compatible report).
 * @property maxP Keep a marker iff `maxP == 1.0 || p <= maxP`, as TASSEL 5's GLM does.
 */
internal data class GlmOptions(
    val marker: MarkerOptions = MarkerOptions(),
    val appendAddDom: Boolean = false,
    val maxP: Double = 1.0,
    val blockSize: Int = MarkerSource.DEFAULT_MARKER_BLOCK_SIZE,
) {
    init {
        require(maxP in 0.0..1.0) { "maxP must be between 0 and 1: $maxP" }
        require(blockSize > 0) { "blockSize must be positive: $blockSize" }
    }
}

/**
 * How close to a combination of earlier columns a GLM column must be to count as aliased: `sqrt(1e-13)`, TASSEL 5's
 * sweep tolerance (`SweepFast.TOL`, a bound on the residual sum of squares) as a bound on the residual norm. TASSEL 5
 * reads phenotype values as `float`, and so does this port (ADR-0023), so a covariate that is the exact sum of two others
 * in the file is the sum only to about `1e-7` once stored, and still has to be recognized as dependent.
 */
internal const val GLM_RANK_TOLERANCE = 3.1622776601683794E-7

/** What a replicated GLM tests the marker against, named in the `message` of its records. */
internal const val TAXA_WITHIN_MARKER = "error term: taxa within marker"

/**
 * The GLM reference kernel: for every trait and every marker of [source], fits TASSEL 5's fixed-effect model and sends
 * one chunk per site block to [sink], records in marker order with each marker's traits together, omnibus test first
 * (ADR-0027). The model is the mean, the factors, the covariates, and then the marker, so the marker's sum of squares
 * is sequential: what it adds after everything else.
 *
 * Every marker is fitted from scratch on its own observations (those with a call and a trait), from a Householder QR
 * of the design, and nothing is inverted (ADR-0025). The sums of squares are read from `Q'y` and so are never a
 * difference of two residual sums. This is the reference path the fast path (tassel-6#192) is checked against: it does
 * not reuse work between markers.
 *
 * - `genotype` is TASSEL 5's omnibus test of the genotype classes, with `add_F` and `dom_F` for biallelic sites
 *   (`DiscreteSitesFELM`). A monomorphic marker is reported, as TASSEL 5 does, with no marker degrees of freedom and
 *   `NaN` statistics, as a `Skipped(Monomorphic)` record.
 * - `dosage` is `ReferenceProbabilityFELM`: one numeric column per marker.
 * - `additive` is the 1-df allele-count test with a named effect allele, and `additiveDominance` adds the
 *   heterozygote indicator; both report a signed effect with its standard error and confidence interval.
 * - With replicated taxa, which TASSEL 5 refuses, the model adds taxa nested within the marker and the marker is tested
 *   against that term's mean square ([TAXA_WITHIN_MARKER]), per decision Q1 and ADR-0026. Allele effects come from the
 *   model without the taxa term.
 *
 * Reports progress after each site block, which is also where an engine cancels the run.
 *
 * @see net.maizegenetics.analysis.association.AbstractFixedEffectLM
 * @see net.maizegenetics.analysis.association.DiscreteSitesFELM
 * @see net.maizegenetics.analysis.association.ReferenceProbabilityFELM
 */
internal suspend fun glm(
    data: AssociationData,
    source: MarkerSource,
    options: GlmOptions,
    sink: AssociationSink,
    progress: ProgressReporter = ProgressReporter.None,
) {
    val fitter = GlmFitter(data, options)
    val capacity = data.observationSets.maxOfOrNull { it.size } ?: 0
    val columns = MarkerColumns(capacity)
    val total = source.numberOfMarkers.toLong()
    var done = 0L
    var sequence = 0
    progress.report(0L, total, "glm")
    for (block in source.scan(data.observationSets, options.marker, options.blockSize)) {
        val records = ArrayList<AssociationRecord>()
        val effects = ArrayList<AssociationEffect>()
        for (marker in block.markers) {
            val position = source.position(marker)
            for (trait in data.traits.indices) {
                val loaded = block.load(marker, trait, columns)
                fitter.analyze(trait, marker, position, columns, loaded, records, effects)
            }
        }
        sink.accept(AssociationChunk(sequence++, records, effects))
        done += block.markers.last - block.markers.first + 1
        progress.report(done, total, "glm")
    }
}

/**
 * Runs [glm] into an in-memory table with the `maxP` filter, the counts of what was tested or skipped, and provenance.
 * [parameters] are recorded as given.
 */
internal suspend fun glmTable(
    data: AssociationData,
    source: MarkerSource,
    options: GlmOptions = GlmOptions(),
    parameters: Map<String, String> = emptyMap(),
    progress: ProgressReporter = ProgressReporter.None,
): AssociationTable {
    val collector = AssociationCollector(options.maxP, MaxPRule.Glm, data.traits)
    glm(data, source, options, collector, progress)
    val excluded = LinkedHashMap<String, Int>()
    data.traits.forEachIndexed { index, trait ->
        for ((reason, count) in data.excludedObservations[index]) excluded["$trait.$reason"] = count
    }
    val provenance = AnalysisProvenance(
        operation = "glm",
        parameters = parameters,
        inputs = listOf(
            AnalysisInput(
                name = data.data.phenotype.name,
                taxa = data.data.genotypeTaxonIndices().size,
                markers = source.numberOfMarkers,
                observations = data.numberOfObservations,
            ),
        ),
        model = describeModel(data, options),
        excludedObservations = excluded,
    )
    return collector.finish(provenance)
}

private fun describeModel(data: AssociationData, options: GlmOptions): String = buildString {
    append("trait ~ 1")
    data.factors.forEach { append(" + factor(").append(it).append(')') }
    data.covariates.forEach { append(" + ").append(it) }
    append(" + marker[").append(options.marker.encoding.name).append(']')
    append("; traits=").append(data.traits.joinToString(","))
}

/** Fits one trait's model at one marker at a time. Holds each trait's base design so it is built once. */
private class GlmFitter(private val data: AssociationData, private val options: GlmOptions) {
    private class Trait(val spec: ModelSpec, val observations: IntArray, val replicated: Boolean) {
        var base: DesignMatrix? = null
    }

    private val traits = data.traits.indices.map {
        Trait(data.baseSpec(it), data.observations(it), data.hasReplicates(it))
    }

    fun analyze(
        traitIndex: Int,
        marker: Int,
        position: Position,
        columns: MarkerColumns,
        loaded: Boolean,
        records: MutableList<AssociationRecord>,
        effects: MutableList<AssociationEffect>,
    ) {
        val trait = traits[traitIndex]
        val identity =
            Identity(data.traits[traitIndex], marker, position.name, position.chromosome.name, position.position)
        val monomorphicGenotype = !loaded &&
            columns.skip == SkipReason.Monomorphic &&
            options.marker.encoding == MarkerEncoding.genotype &&
            columns.size > 0
        if (!loaded && !monomorphicGenotype) {
            records += skipped(identity, columns)
            return
        }
        val n = columns.size
        val base = baseDesign(trait, columns)
        val taxonOf = if (trait.replicated) {
            IntArray(n) { data.taxonIds[trait.observations[columns.observation[it]]] }
        } else {
            null
        }
        val message = if (trait.replicated) TAXA_WITHIN_MARKER else ""
        val context = Context(identity, columns, base, taxonOf, message)
        when (options.marker.encoding) {
            MarkerEncoding.genotype -> genotype(context, monomorphicGenotype, records, effects)
            MarkerEncoding.dosage -> dosage(context, records, effects)
            MarkerEncoding.additive -> additive(context, records)
            MarkerEncoding.additiveDominance -> additiveDominance(context, records)
        }
    }

    private class Identity(
        val trait: String,
        val marker: Int,
        val name: String,
        val chromosome: String,
        val position: Int,
    )

    private class Context(
        val id: Identity,
        val columns: MarkerColumns,
        val base: DesignMatrix,
        val taxonOf: IntArray?,
        val message: String,
    ) {
        val n: Int get() = columns.size
    }

    private fun skipped(id: Identity, columns: MarkerColumns) = AssociationRecord(
        trait = id.trait,
        markerIndex = id.marker,
        marker = id.name,
        chromosome = id.chromosome,
        position = id.position,
        encoding = options.marker.encoding,
        minorObs = columns.minorObs,
        status = AssociationStatus.Skipped,
        skipReason = columns.skip,
    )

    private fun baseDesign(trait: Trait, columns: MarkerColumns): DesignMatrix {
        val n = columns.size
        if (n == trait.observations.size) {
            return trait.base
                ?: DesignMatrix.build(trait.spec, data.data.phenotype, trait.observations, GLM_RANK_TOLERANCE)
                    .also { trait.base = it }
        }
        return DesignMatrix.build(
            trait.spec,
            data.data.phenotype,
            IntArray(n) { trait.observations[columns.observation[it]] },
            GLM_RANK_TOLERANCE,
        )
    }

    // --- genotype -------------------------------------------------------------------------------------------------

    private fun genotype(
        c: Context,
        monomorphic: Boolean,
        records: MutableList<AssociationRecord>,
        effects: MutableList<AssociationEffect>,
    ) {
        val k = c.columns.numberOfClasses
        val markerColumns = (0 until k - 1).map { cls ->
            c.columns.classLabels[cls] to DoubleArray(c.n) { if (c.columns.classIndex[it] == cls) 1.0 else 0.0 }
        }
        val full = fitBlock(c.base, markerColumns, c.taxonOf)
        val p0 = c.base.columnCount
        val fit = LinearModelFit.of(full.design)
        val test = full.test()

        // add_F and dom_F: the base model plus a copy count per allele but the last (DiscreteSitesFELM).
        val withAlleles = c.columns.columnCount > 0
        var additive: BlockResult? = null
        var addTest = Fstat.NONE
        var domTest = Fstat.NONE
        if (withAlleles) {
            val addColumns = (0 until c.columns.columnCount).map { j ->
                c.columns.columnLabels[j] to c.columns.column(j).copyOf(c.n)
            }
            additive = fitBlock(c.base, addColumns, c.taxonOf)
            // TASSEL 5 turns a non-finite add or dom statistic, and its p-value, into NaN.
            addTest = additive.test().strict()
            domTest = Fstat.of(full.ss - additive.ss, full.df - additive.df, full.errorSs, full.errorDf).strict()
        }

        var addEffect = Double.NaN
        var domEffect = Double.NaN
        if (options.appendAddDom && withAlleles && c.columns.columnCount == 1) {
            // Q6: TASSEL 5 indexes the marker's coefficients by the number of effects, not columns, which is wrong
            // when a factor has other than two levels; these are the marker's own coefficients.
            val b0 = coefficient(fit, p0)
            val b1 = coefficient(fit, p0 + 1)
            if (k == 2) {
                addEffect = abs(b0 / 2)
            } else if (k == 3) {
                addEffect = abs((b0 - b1) / 2)
                domEffect = -(b0 + b1) / 2
            }
        }

        val main = AssociationRecord(
            trait = c.id.trait,
            markerIndex = c.id.marker,
            marker = c.id.name,
            chromosome = c.id.chromosome,
            position = c.id.position,
            test = AssociationTestKind.Genotype,
            encoding = MarkerEncoding.genotype,
            observations = c.n,
            taxa = taxa(c),
            minorObs = c.columns.minorObs,
            statistic = test.f,
            numeratorDf = full.df.toDouble(),
            denominatorDf = full.errorDf.toDouble(),
            p = test.p,
            log10P = test.log10P,
            rSquared = full.ss / full.totalSs,
            markerMeanSquare = full.ss / full.df,
            errorMeanSquare = full.errorSs / full.errorDf,
            modelDf = full.modelDf.toDouble(),
            modelMeanSquare = full.modelSs / full.modelDf,
            status = if (monomorphic) AssociationStatus.Skipped else AssociationStatus.Ok,
            skipReason = if (monomorphic) SkipReason.Monomorphic else null,
            message = c.message,
        )
        records += main
        if (additive != null) {
            records += main.copy(
                test = AssociationTestKind.Additive,
                statistic = addTest.f,
                numeratorDf = additive.df.toDouble(),
                denominatorDf = additive.errorDf.toDouble(),
                p = addTest.p,
                log10P = addTest.log10P,
                effect = addEffect,
                rSquared = Double.NaN,
                markerMeanSquare = Double.NaN,
                errorMeanSquare = Double.NaN,
                modelDf = Double.NaN,
                modelMeanSquare = Double.NaN,
            )
            records += main.copy(
                test = AssociationTestKind.Dominance,
                statistic = domTest.f,
                numeratorDf = (full.df - additive.df).toDouble(),
                denominatorDf = full.errorDf.toDouble(),
                p = domTest.p,
                log10P = domTest.log10P,
                effect = domEffect,
                rSquared = Double.NaN,
                markerMeanSquare = Double.NaN,
                errorMeanSquare = Double.NaN,
                modelDf = Double.NaN,
                modelMeanSquare = Double.NaN,
            )
        }
        if (withAlleles) {
            for (cls in 0 until k) {
                val last = cls == k - 1
                effects += AssociationEffect(
                    trait = c.id.trait,
                    markerIndex = c.id.marker,
                    marker = c.id.name,
                    chromosome = c.id.chromosome,
                    position = c.id.position,
                    allele = c.columns.classLabels[cls],
                    observations = c.columns.classCounts[cls],
                    estimate = if (last) 0.0 else coefficient(fit, p0 + cls),
                    standardError = if (last) Double.NaN else fit.standardErrors[p0 + cls],
                )
            }
        }
    }

    /** TASSEL 5 writes 0 for a coefficient it swept out; an aliased column here is `NaN`. */
    private fun coefficient(fit: LinearModelFit, column: Int): Double =
        if (column < fit.coefficients.size) fit.coefficients[column].let { if (it.isNaN()) 0.0 else it } else 0.0

    // --- dosage ---------------------------------------------------------------------------------------------------

    private fun dosage(c: Context, records: MutableList<AssociationRecord>, effects: MutableList<AssociationEffect>) {
        val dosage = listOf(c.columns.columnLabels[0] to c.columns.column(0).copyOf(c.n))
        val block = fitBlock(c.base, dosage, c.taxonOf)
        val fit = LinearModelFit.of(block.design)
        val p0 = c.base.columnCount
        val test = block.test()
        records += AssociationRecord(
            trait = c.id.trait,
            markerIndex = c.id.marker,
            marker = c.id.name,
            chromosome = c.id.chromosome,
            position = c.id.position,
            test = AssociationTestKind.Dosage,
            encoding = MarkerEncoding.dosage,
            observations = c.n,
            taxa = taxa(c),
            effect = fit.coefficients[p0],
            standardError = fit.standardErrors[p0],
            statistic = test.f,
            numeratorDf = block.df.toDouble(),
            denominatorDf = block.errorDf.toDouble(),
            p = test.p,
            log10P = test.log10P,
            rSquared = block.ss / block.totalSs,
            markerMeanSquare = block.ss / block.df,
            errorMeanSquare = block.errorSs / block.errorDf,
            modelDf = block.modelDf.toDouble(),
            modelMeanSquare = block.modelSs / block.modelDf,
            message = c.message,
        )
        effects += AssociationEffect(
            trait = c.id.trait,
            markerIndex = c.id.marker,
            marker = c.id.name,
            chromosome = c.id.chromosome,
            position = c.id.position,
            estimate = coefficient(fit, p0),
            standardError = fit.standardErrors[p0],
        )
    }

    // --- additive and additive-dominance ----------------------------------------------------------------------------

    private fun additive(c: Context, records: MutableList<AssociationRecord>) {
        val labels = c.columns.columnLabels
        val columns = labels.indices.map { labels[it] to c.columns.column(it).copyOf(c.n) }
        val block = fitBlock(c.base, columns, c.taxonOf)
        val fit = LinearModelFit.of(block.design)
        val single = labels.size == 1
        records += encodedRecord(
            c,
            block,
            block.test(),
            AssociationTestKind.Additive,
            fit,
            effectColumn = if (single) c.base.columnCount else -1,
            effectAllele = labels.joinToString(","),
            singleAllele = single,
        )
    }

    private fun additiveDominance(c: Context, records: MutableList<AssociationRecord>) {
        val labels = c.columns.columnLabels
        val addColumns = listOf(labels[0] to c.columns.column(0).copyOf(c.n))
        val domColumns = listOf(labels[1] to c.columns.column(1).copyOf(c.n))
        val full = fitBlock(c.base, addColumns + domColumns, c.taxonOf)
        val additive = fitBlock(c.base, addColumns, c.taxonOf)
        val fit = LinearModelFit.of(full.design)
        val p0 = c.base.columnCount
        // As TASSEL 5: add against the additive model's own error, dom against the full model's.
        val joint = encodedRecord(
            c,
            full,
            full.test(),
            AssociationTestKind.Genotype,
            fit,
            -1,
            labels.joinToString(","),
            true,
        )
        records += joint
        records += encodedRecord(c, additive, additive.test(), AssociationTestKind.Additive, fit, p0, labels[0], true)
        val domTest = Fstat.of(full.ss - additive.ss, full.df - additive.df, full.errorSs, full.errorDf)
        records += encodedRecord(c, full, domTest, AssociationTestKind.Dominance, fit, p0 + 1, labels[1], true)
            .copy(numeratorDf = (full.df - additive.df).toDouble(), comparisonAllele = "")
    }

    private fun encodedRecord(
        c: Context,
        block: BlockResult,
        test: Fstat,
        kind: AssociationTestKind,
        fit: LinearModelFit,
        effectColumn: Int,
        effectAllele: String,
        singleAllele: Boolean,
    ): AssociationRecord {
        val copies = if (singleAllele) c.columns.column(0).sumOfFirst(c.n) else Double.NaN
        val alleles = 2.0 * c.n
        val comparison = c.columns.comparisonAllele
        val effect = if (effectColumn >= 0) fit.coefficients[effectColumn] else Double.NaN
        val standardError = if (effectColumn >= 0) fit.standardErrors[effectColumn] else Double.NaN
        var lower = Double.NaN
        var upper = Double.NaN
        if (!standardError.isNaN() && fit.residualDf > 0) {
            val t = StudentTDistribution(fit.residualDf.toDouble()).quantile(0.975)
            lower = effect - t * standardError
            upper = effect + t * standardError
        }
        return AssociationRecord(
            trait = c.id.trait,
            markerIndex = c.id.marker,
            marker = c.id.name,
            chromosome = c.id.chromosome,
            position = c.id.position,
            test = kind,
            encoding = options.marker.encoding,
            effectAllele = effectAllele,
            comparisonAllele = if (comparison == Allele.UNKNOWN) "" else Nucleotide.alleleString(comparison),
            observations = c.n,
            taxa = taxa(c),
            alleleCount = if (singleAllele) copies.toInt() else -1,
            alleleFrequency = if (singleAllele) copies / alleles else Double.NaN,
            minorAlleleCount = if (singleAllele) minOf(copies, alleles - copies).toInt() else -1,
            minorObs = c.columns.minorObs,
            effect = effect,
            standardError = standardError,
            ciLower = lower,
            ciUpper = upper,
            statistic = test.f,
            numeratorDf = block.df.toDouble(),
            denominatorDf = block.errorDf.toDouble(),
            p = test.p,
            log10P = test.log10P,
            rSquared = block.ss / block.totalSs,
            markerMeanSquare = block.ss / block.df,
            errorMeanSquare = block.errorSs / block.errorDf,
            modelDf = block.modelDf.toDouble(),
            modelMeanSquare = block.modelSs / block.modelDf,
            message = c.message,
        )
    }

    private fun DoubleArray.sumOfFirst(count: Int): Double {
        var sum = 0.0
        for (i in 0 until count) sum += this[i]
        return sum
    }

    private fun taxa(c: Context): Int {
        val taxonOf = c.taxonOf ?: return c.n
        return taxonOf.toSet().size
    }

    // --- fitting ----------------------------------------------------------------------------------------------------

    /** Sums of squares and degrees of freedom of the block of columns added after the base model, and its error term. */
    private class BlockResult(
        val design: DesignMatrix,
        val ss: Double,
        val df: Int,
        val errorSs: Double,
        val errorDf: Int,
        val modelSs: Double,
        val modelDf: Int,
        val totalSs: Double,
    ) {
        fun test(): Fstat = Fstat.of(ss, df, errorSs, errorDf)
    }

    /**
     * Fits `[base | block]`, plus taxa nested in the block when [taxonOf] is given, and reads the block's sequential sum
     * of squares from `Q'y`. The error term is the residual, or with replicated taxa the taxa term added after the block.
     */
    private fun fitBlock(base: DesignMatrix, block: List<Pair<String, DoubleArray>>, taxonOf: IntArray?): BlockResult {
        val p0 = base.columnCount
        val withBlock = extend(base, ModelTerm.MarkerTerm(BLOCK_LABEL), block)
        val p1 = p0 + block.size
        val n = base.rows
        if (taxonOf == null) {
            val sums = Sums(withBlock)
            val r0 = sums.rankThrough(p0)
            val r1 = sums.rankThrough(p1)
            return BlockResult(
                withBlock,
                sums.ss(r0, r1),
                r1 - r0,
                sums.ss(sums.rank, n),
                n - sums.rank,
                sums.ss(1, sums.rank),
                sums.rank - 1,
                sums.ss(1, n),
            )
        }
        val withTaxa = extend(withBlock, ModelTerm.Factor(TAXON_LABEL), taxonColumns(taxonOf))
        val sums = Sums(withTaxa)
        val r0 = sums.rankThrough(p0)
        val r1 = sums.rankThrough(p1)
        val r2 = sums.rank
        return BlockResult(
            withBlock,
            sums.ss(r0, r1),
            r1 - r0,
            sums.ss(r1, r2),
            r2 - r1,
            sums.ss(1, r2),
            r2 - 1,
            sums.ss(1, n),
        )
    }

    /** `base` with [columns] appended as one term. Each array holds at least `base.rows` values. */
    private fun extend(base: DesignMatrix, term: ModelTerm, columns: List<Pair<String, DoubleArray>>): DesignMatrix {
        val n = base.rows
        val p = base.columnCount
        val matrix = DenseMatrix(n, p + columns.size)
        base.matrix.data.copyInto(matrix.data, 0, 0, n * p)
        columns.forEachIndexed { j, (_, values) -> values.copyInto(matrix.data, (p + j) * n, 0, n) }
        val termIndex = base.spec.terms.size
        return DesignMatrix(
            ModelSpec(base.spec.response, base.spec.terms + term),
            base.observations,
            matrix,
            base.response,
            base.columns + columns.map { DesignColumn(it.first, termIndex) },
            base.rankTolerance,
        )
    }

    /** One indicator column per taxon but the last present, as a restricted factor (TASSEL 5's `taxaEffect`). */
    private fun taxonColumns(taxonOf: IntArray): List<Pair<String, DoubleArray>> {
        val present = taxonOf.distinct().sorted()
        return present.dropLast(1).map { id ->
            "$TAXON_LABEL=$id" to DoubleArray(taxonOf.size) { if (taxonOf[it] == id) 1.0 else 0.0 }
        }
    }

    /** `Q'y` of a design's QR: the sequential sums of squares are sums of its squared entries. */
    private class Sums(design: DesignMatrix) {
        private val decomposition = design.decomposition
        private val qty = design.response.copyOf().also { decomposition.qr?.applyQt(it) }
        val rank: Int get() = decomposition.rank

        /** The rank of the columns before [columnEnd]: how many of them are not aliased. */
        fun rankThrough(columnEnd: Int): Int = decomposition.kept.count { it < columnEnd }

        /** The sum of squares added by the rank positions [from] until [to]. */
        fun ss(from: Int, to: Int): Double {
            var sum = 0.0
            for (i in from until minOf(to, qty.size)) sum += qty[i] * qty[i]
            return sum
        }
    }

    /** An F statistic and its p-values, computed as TASSEL 5 does: `((SS / df) / errorSS) * errorDf`. */
    private class Fstat(val f: Double, val p: Double, val log10P: Double) {
        fun strict(): Fstat = if (f.isFinite()) this else NONE

        companion object {
            val NONE = Fstat(Double.NaN, Double.NaN, Double.NaN)

            fun of(ss: Double, df: Int, errorSs: Double, errorDf: Int): Fstat {
                val f = ss / df / errorSs * errorDf
                if (!f.isFinite()) {
                    // TASSEL 5 turns a non-finite statistic into NaN for add and dom, and its p into NaN throughout;
                    // an infinite F with a positive error df is a perfect fit, which has p = 0.
                    if (f == Double.POSITIVE_INFINITY &&
                        df > 0 &&
                        errorDf > 0
                    ) {
                        return Fstat(f, 0.0, Double.NEGATIVE_INFINITY)
                    }
                    return Fstat(f, Double.NaN, Double.NaN)
                }
                val p = try {
                    fTestPValue(f, df.toDouble(), errorDf.toDouble())
                } catch (e: IllegalArgumentException) {
                    Double.NaN
                }
                val log10P = try {
                    fTestLog10PValue(f, df.toDouble(), errorDf.toDouble())
                } catch (e: IllegalArgumentException) {
                    Double.NaN
                }
                return Fstat(f, p, log10P)
            }
        }
    }

    private companion object {
        const val BLOCK_LABEL = "(marker)"
        const val TAXON_LABEL = "(taxon)"
    }
}
