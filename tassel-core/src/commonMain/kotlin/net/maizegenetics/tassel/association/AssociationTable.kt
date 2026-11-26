package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.table.ResultTable

/**
 * Counts of what was tested, always kept in full even when a filter keeps only some rows (ADR-0027). Every record
 * counts once, so an additive or dominance test counts beside its omnibus test.
 */
internal data class AssociationSummary(
    val records: Long = 0,
    val tested: Long = 0,
    val skipped: Map<SkipReason, Long> = emptyMap(),
    val failed: Long = 0,
) {
    val skippedTotal: Long get() = skipped.values.sum()
}

/**
 * Results of an association scan in primitive columns: about 8 bytes per numeric cell and no object per row, so a
 * million markers cost a few hundred megabytes rather than gigabytes (ADR-0027). Strings are stored once per distinct
 * value.
 *
 * Read-only; make one with [AssociationTableBuilder]. Rows are in the order they were appended, which for a scan is
 * (trait, marker) order; reports write them trait by trait.
 */
internal class AssociationTable internal constructor(
    private val columns: RecordColumns,
    private val effectColumns: EffectColumns,
    private val strings: StringDictionary,
    val provenance: AnalysisProvenance,
    val summary: AssociationSummary,
    private val reservedTraits: List<String> = emptyList(),
) {
    /** Number of records. */
    val size: Int get() = columns.size

    /** Number of effect rows. */
    val effectCount: Int get() = effectColumns.size

    /** Record [index], rebuilt from the columns. */
    operator fun get(index: Int): AssociationRecord {
        require(index in 0 until size) { "Record $index of $size" }
        return columns.record(index, strings)
    }

    /** Effect row [index], rebuilt from the columns. */
    fun effect(index: Int): AssociationEffect {
        require(index in 0 until effectCount) { "Effect $index of $effectCount" }
        return effectColumns.effect(index, strings)
    }

    val records: Sequence<AssociationRecord> get() = (0 until size).asSequence().map { get(it) }

    val effects: Sequence<AssociationEffect> get() = (0 until effectCount).asSequence().map { effect(it) }

    /** Traits in the order reserved by the analysis, then in order of first appearance. */
    val traits: List<String> by lazy {
        val seen = LinkedHashSet<Int>()
        for (trait in reservedTraits) seen.add(strings.code(trait))
        for (row in 0 until size) seen.add(columns.traitCode(row))
        for (row in 0 until effectCount) seen.add(effectColumns.traitCode(row))
        seen.map { strings[it] }
    }

    /** Record rows trait by trait (TASSEL 5's report order), each trait in append order. */
    fun recordOrder(): IntArray = traitMajor(size) { columns.traitCode(it) }

    /** Effect rows trait by trait, each trait in append order. */
    fun effectOrder(): IntArray = traitMajor(effectCount) { effectColumns.traitCode(it) }

    private fun traitMajor(count: Int, traitCode: (Int) -> Int): IntArray {
        val rank = HashMap<Int, Int>()
        for (trait in traits) rank[strings.code(trait)] = rank.size
        val start = IntArray(rank.size + 1)
        for (row in 0 until count) start[rank.getValue(traitCode(row)) + 1]++
        for (i in 1 until start.size) start[i] += start[i - 1]
        val order = IntArray(count)
        for (row in 0 until count) order[start[rank.getValue(traitCode(row))]++] = row
        return order
    }

    /**
     * What the table holds in memory, from the sizes of its arrays: exact for the numeric columns, an estimate for the
     * strings (see [StringDictionary.estimatedBytes]).
     */
    val estimatedBytes: Long
        get() = columns.bytes + effectColumns.bytes + strings.estimatedBytes

    /** A small view for display: the main columns of the first [maxRows] records, as a [ResultTable]. */
    fun toResultTable(title: String = "Association", maxRows: Int = 1000): ResultTable {
        val names = listOf("Trait", "Marker", "Chr", "Pos", "Test", "Status", "Statistic", "df1", "df2", "p", "log10_p")
        val rows = (0 until minOf(size, maxRows)).map { row ->
            val r = get(row)
            listOf(
                r.trait, r.marker, r.chromosome, r.position, r.test.name, r.status.name,
                r.statistic, r.numeratorDf, r.denominatorDf, r.p, r.log10P,
            )
        }
        return ResultTable(title, names, rows)
    }

    internal companion object {
        /** Builds a table from [records] and [effects] in the order given. */
        fun of(
            records: List<AssociationRecord>,
            effects: List<AssociationEffect> = emptyList(),
            provenance: AnalysisProvenance,
            summary: AssociationSummary = AssociationSummary(),
        ): AssociationTable {
            val builder = AssociationTableBuilder()
            records.forEach(builder::append)
            effects.forEach(builder::append)
            return builder.build(provenance, summary)
        }
    }
}

/** Collects records and effects into an [AssociationTable]. Not thread safe: feed it from one coroutine at a time. */
internal class AssociationTableBuilder {
    private val strings = StringDictionary()
    private val columns = RecordColumns()
    private val effectColumns = EffectColumns()
    private var built = false
    private var reserved: List<String> = emptyList()

    val size: Int get() = columns.size

    /** Lists [traits] first, in this order, in [AssociationTable.traits], whether or not they get a row. */
    fun reserveTraits(traits: List<String>) {
        check(columns.size == 0) { "Traits must be reserved before the first record" }
        reserved = traits.distinct()
        reserved.forEach { strings.code(it) }
    }

    fun append(record: AssociationRecord) {
        check(!built) { "The table has been built" }
        columns.add(record, strings)
    }

    fun append(effect: AssociationEffect) {
        check(!built) { "The table has been built" }
        effectColumns.add(effect, strings)
    }

    /** Seals the columns and returns the table. The builder cannot be used again. */
    fun build(provenance: AnalysisProvenance, summary: AssociationSummary = AssociationSummary()): AssociationTable {
        check(!built) { "The table has been built" }
        built = true
        columns.seal()
        effectColumns.seal()
        return AssociationTable(columns, effectColumns, strings, provenance, summary, reserved)
    }
}

internal class RecordColumns {
    private val trait = IntColumn()
    private val markerIndex = IntColumn()
    private val marker = IntColumn()
    private val chromosome = IntColumn()
    private val position = IntColumn()
    private val model = ByteColumn()
    private val test = ByteColumn()
    private val encoding = ByteColumn()
    private val effectAllele = IntColumn()
    private val comparisonAllele = IntColumn()
    private val observations = IntColumn()
    private val taxa = IntColumn()
    private val alleleCount = IntColumn()
    private val alleleFrequency = DoubleColumn()
    private val minorAlleleCount = IntColumn()
    private val minorObs = IntColumn()
    private val effect = DoubleColumn()
    private val standardError = DoubleColumn()
    private val ciLower = DoubleColumn()
    private val ciUpper = DoubleColumn()
    private val statistic = DoubleColumn()
    private val numeratorDf = DoubleColumn()
    private val denominatorDf = DoubleColumn()
    private val p = DoubleColumn()
    private val log10P = DoubleColumn()
    private val rSquared = DoubleColumn()
    private val markerMeanSquare = DoubleColumn()
    private val errorMeanSquare = DoubleColumn()
    private val modelDf = DoubleColumn()
    private val modelMeanSquare = DoubleColumn()
    private val permutationP = DoubleColumn()
    private val geneticVariance = DoubleColumn()
    private val residualVariance = DoubleColumn()
    private val heritability = DoubleColumn()
    private val varianceEstimation = ByteColumn()
    private val status = ByteColumn()
    private val skipReason = ByteColumn()
    private val message = IntColumn()

    private val doubles = listOf(
        alleleFrequency, effect, standardError, ciLower, ciUpper, statistic, numeratorDf, denominatorDf, p, log10P,
        rSquared, markerMeanSquare, errorMeanSquare, modelDf, modelMeanSquare, permutationP, geneticVariance,
        residualVariance, heritability,
    )
    private val ints = listOf(
        trait, markerIndex, marker, chromosome, position, effectAllele, comparisonAllele, observations, taxa,
        alleleCount, minorAlleleCount, minorObs, message,
    )
    private val byteColumns = listOf(model, test, encoding, varianceEstimation, status, skipReason)

    val size: Int get() = trait.size

    fun traitCode(row: Int): Int = trait[row]

    fun add(r: AssociationRecord, strings: StringDictionary) {
        trait.add(strings.code(r.trait))
        markerIndex.add(r.markerIndex)
        marker.add(strings.code(r.marker))
        chromosome.add(strings.code(r.chromosome))
        position.add(r.position)
        model.add(r.model.ordinal)
        test.add(r.test.ordinal)
        encoding.add(r.encoding.ordinal)
        effectAllele.add(strings.code(r.effectAllele))
        comparisonAllele.add(strings.code(r.comparisonAllele))
        observations.add(r.observations)
        taxa.add(r.taxa)
        alleleCount.add(r.alleleCount)
        alleleFrequency.add(r.alleleFrequency)
        minorAlleleCount.add(r.minorAlleleCount)
        minorObs.add(r.minorObs)
        effect.add(r.effect)
        standardError.add(r.standardError)
        ciLower.add(r.ciLower)
        ciUpper.add(r.ciUpper)
        statistic.add(r.statistic)
        numeratorDf.add(r.numeratorDf)
        denominatorDf.add(r.denominatorDf)
        p.add(r.p)
        log10P.add(r.log10P)
        rSquared.add(r.rSquared)
        markerMeanSquare.add(r.markerMeanSquare)
        errorMeanSquare.add(r.errorMeanSquare)
        modelDf.add(r.modelDf)
        modelMeanSquare.add(r.modelMeanSquare)
        permutationP.add(r.permutationP)
        geneticVariance.add(r.geneticVariance)
        residualVariance.add(r.residualVariance)
        heritability.add(r.heritability)
        varianceEstimation.add(r.varianceEstimation.ordinal)
        status.add(r.status.ordinal)
        skipReason.add(r.skipReason?.ordinal ?: NO_REASON)
        message.add(strings.code(r.message))
    }

    fun record(row: Int, strings: StringDictionary) = AssociationRecord(
        trait = strings[trait[row]],
        markerIndex = markerIndex[row],
        marker = strings[marker[row]],
        chromosome = strings[chromosome[row]],
        position = position[row],
        model = AssociationModel.entries[model[row]],
        test = AssociationTestKind.entries[test[row]],
        encoding = MarkerEncoding.entries[encoding[row]],
        effectAllele = strings[effectAllele[row]],
        comparisonAllele = strings[comparisonAllele[row]],
        observations = observations[row],
        taxa = taxa[row],
        alleleCount = alleleCount[row],
        alleleFrequency = alleleFrequency[row],
        minorAlleleCount = minorAlleleCount[row],
        minorObs = minorObs[row],
        effect = effect[row],
        standardError = standardError[row],
        ciLower = ciLower[row],
        ciUpper = ciUpper[row],
        statistic = statistic[row],
        numeratorDf = numeratorDf[row],
        denominatorDf = denominatorDf[row],
        p = p[row],
        log10P = log10P[row],
        rSquared = rSquared[row],
        markerMeanSquare = markerMeanSquare[row],
        errorMeanSquare = errorMeanSquare[row],
        modelDf = modelDf[row],
        modelMeanSquare = modelMeanSquare[row],
        permutationP = permutationP[row],
        geneticVariance = geneticVariance[row],
        residualVariance = residualVariance[row],
        heritability = heritability[row],
        varianceEstimation = VarianceEstimation.entries[varianceEstimation[row]],
        status = AssociationStatus.entries[status[row]],
        skipReason = skipReason[row].let { if (it == NO_REASON) null else SkipReason.entries[it] },
        message = strings[message[row]],
    )

    fun seal() {
        doubles.forEach { it.seal() }
        ints.forEach { it.seal() }
        byteColumns.forEach { it.seal() }
    }

    val bytes: Long
        get() = doubles.sumOf { it.bytes } + ints.sumOf { it.bytes } + byteColumns.sumOf { it.bytes }

    private companion object {
        const val NO_REASON = -1
    }
}

internal class EffectColumns {
    private val trait = IntColumn()
    private val markerIndex = IntColumn()
    private val marker = IntColumn()
    private val chromosome = IntColumn()
    private val position = IntColumn()
    private val allele = IntColumn()
    private val observations = IntColumn()
    private val estimate = DoubleColumn()
    private val standardError = DoubleColumn()

    val size: Int get() = trait.size

    fun traitCode(row: Int): Int = trait[row]

    fun add(e: AssociationEffect, strings: StringDictionary) {
        trait.add(strings.code(e.trait))
        markerIndex.add(e.markerIndex)
        marker.add(strings.code(e.marker))
        chromosome.add(strings.code(e.chromosome))
        position.add(e.position)
        allele.add(strings.code(e.allele))
        observations.add(e.observations)
        estimate.add(e.estimate)
        standardError.add(e.standardError)
    }

    fun effect(row: Int, strings: StringDictionary) = AssociationEffect(
        trait = strings[trait[row]],
        markerIndex = markerIndex[row],
        marker = strings[marker[row]],
        chromosome = strings[chromosome[row]],
        position = position[row],
        allele = strings[allele[row]],
        observations = observations[row],
        estimate = estimate[row],
        standardError = standardError[row],
    )

    fun seal() {
        listOf(trait, markerIndex, marker, chromosome, position, allele, observations).forEach { it.seal() }
        estimate.seal()
        standardError.seal()
    }

    val bytes: Long
        get() = listOf(trait, markerIndex, marker, chromosome, position, allele, observations).sumOf { it.bytes } +
            estimate.bytes + standardError.bytes
}
