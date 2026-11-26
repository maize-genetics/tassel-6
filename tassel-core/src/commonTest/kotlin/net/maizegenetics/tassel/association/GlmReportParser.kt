package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.table.readsBackAs
import kotlin.math.log10
import kotlin.math.nextDown
import kotlin.math.nextUp

/**
 * Reads a TASSEL 5 GLM report written at full precision back into records, the way a kernel would produce them, so
 * the writers can be checked against TASSEL 5's own files. Test code only.
 */
internal class ParsedGlm(val records: List<AssociationRecord>, val effects: List<AssociationEffect>)

internal class ParsedLayout(val layout: GlmLayout, val permutations: Boolean, val appendAddDom: Boolean)

internal fun detectLayout(statisticsHeader: String): ParsedLayout {
    val names = statisticsHeader.split('\t')
    return ParsedLayout(
        layout = if ("add_F" in names) GlmLayout.Genotype else GlmLayout.ReferenceProbability,
        permutations = "perm_p" in names,
        appendAddDom = "addEffect" in names,
    )
}

internal fun parseGlm(statistics: String, effects: String): ParsedGlm {
    val lines = statistics.trimEnd('\n').split('\n')
    val header = lines.first().split('\t')
    val genotype = "add_F" in header
    fun column(name: String) = header.indexOf(name)

    val records = ArrayList<AssociationRecord>()
    val markerIndexes = HashMap<Pair<String, String>, Int>()
    for ((row, line) in lines.drop(1).withIndex()) {
        val cells = line.split('\t')
        fun double(name: String): Double = exactDouble(cells[column(name)])
        val trait = cells[0]
        val marker = cells[1]
        markerIndexes[trait to marker] = row
        val p = double("p")
        fun record(test: AssociationTestKind, statistic: Double, pValue: Double, effect: Double) = AssociationRecord(
            trait = trait,
            markerIndex = row,
            marker = marker,
            chromosome = cells[2],
            position = cells[3].toInt(),
            test = test,
            observations = -1,
            minorObs = if (genotype) cells[column("minorObs")].toInt() else -1,
            effect = effect,
            statistic = statistic,
            numeratorDf = double("marker_df"),
            denominatorDf = double("error_df"),
            p = pValue,
            log10P = log10(pValue),
            rSquared = double("marker_Rsq"),
            markerMeanSquare = double("marker_MS"),
            errorMeanSquare = double("error_MS"),
            modelDf = double("model_df"),
            modelMeanSquare = double("model_MS"),
            permutationP = if ("perm_p" in header) double("perm_p") else Double.NaN,
        )
        records += record(
            if (genotype) AssociationTestKind.Genotype else AssociationTestKind.Dosage,
            double("marker_F"),
            p,
            Double.NaN,
        )
        if (genotype) {
            records += record(
                AssociationTestKind.Additive,
                double("add_F"),
                double("add_p"),
                if ("addEffect" in header) double("addEffect") else Double.NaN,
            )
            records += record(
                AssociationTestKind.Dominance,
                double("dom_F"),
                double("dom_p"),
                if ("domEffect" in header) double("domEffect") else Double.NaN,
            )
        }
    }

    val effectLines = effects.trimEnd('\n').split('\n')
    val effectHeader = effectLines.first().split('\t')
    val withClasses = "Allele" in effectHeader
    val parsedEffects = effectLines.drop(1).map { line ->
        val cells = line.split('\t')
        AssociationEffect(
            trait = cells[0],
            markerIndex = markerIndexes.getValue(cells[0] to cells[1]),
            marker = cells[1],
            chromosome = cells[2],
            position = cells[3].toInt(),
            allele = if (withClasses) cells[5] else "",
            observations = if (withClasses) cells[4].toInt() else -1,
            estimate = exactDouble(cells.last()),
        )
    }
    return ParsedGlm(records, parsedEffects)
}

/**
 * [text] as the nearest `Double`, on every target. Kotlin/Wasm's `String.toDouble` is off by one ulp for some 16 and
 * 17 digit strings, which would make a round trip through TASSEL 5's full-precision files differ in the last digit, so
 * the result is checked against the text exactly and moved to a neighbour if it does not read back.
 */
internal fun exactDouble(text: String): Double {
    val guess = text.toDouble()
    if (guess.isNaN() || guess.isInfinite() || guess == 0.0) return guess
    val unsigned = text.trimStart('-', '+')
    val mantissa = unsigned.substringBefore('E').substringBefore('e')
    val exponent = unsigned.substring(mantissa.length).drop(1).ifEmpty { "0" }.toInt()
    val point = mantissa.indexOf('.')
    val fractionDigits = if (point < 0) 0 else mantissa.length - point - 1
    val digits = mantissa.replace(".", "").toLong()
    for (candidate in listOf(guess, guess.nextDown(), guess.nextUp())) {
        if (readsBackAs(candidate, digits, exponent - fractionDigits)) return candidate
    }
    error("$text does not read back as $guess or a neighbour")
}
