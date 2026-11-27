package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.tasselVersion

/** One input to an analysis: its name and dimensions; a dimension that is unknown is `-1`. */
internal data class AnalysisInput(
    val name: String,
    val taxa: Int = -1,
    val markers: Int = -1,
    val observations: Int = -1,
)

/**
 * Where a result table came from (ADR-0027): enough to repeat the run and to know what was left out.
 *
 * @property operation The operation id, such as `glm`.
 * @property parameters Bound parameter values, in the order given.
 * @property model The model in canonical text.
 * @property kinship The kinship source and scaling, or empty.
 * @property seed The random seed (permutations), or null.
 * @property permutations The permutation count and per-trait experiment-wise thresholds, or null.
 * @property excludedObservations Observations excluded, by reason.
 * @property excludedMarkers Markers excluded, by reason.
 */
internal data class AnalysisProvenance(
    val operation: String,
    val version: String = tasselVersion(),
    val parameters: Map<String, String> = emptyMap(),
    val inputs: List<AnalysisInput> = emptyList(),
    val model: String = "",
    val kinship: String = "",
    val seed: Long? = null,
    val permutations: PermutationSummary? = null,
    val excludedObservations: Map<String, Int> = emptyMap(),
    val excludedMarkers: Map<String, Int> = emptyMap(),
) {
    /**
     * The `##key=value` header lines, as `DistanceMatrix` writes annotations. Line breaks and tabs in a value become
     * spaces so a value cannot break the header.
     */
    fun headerLines(): List<String> = buildList {
        add(line("tassel_version", version))
        add(line("operation", operation))
        for ((name, value) in parameters) add(line("parameter.$name", value))
        for (input in inputs) {
            add(
                line(
                    "input",
                    buildString {
                        append(input.name)
                        if (input.taxa >= 0) append(" taxa=").append(input.taxa)
                        if (input.markers >= 0) append(" markers=").append(input.markers)
                        if (input.observations >= 0) append(" observations=").append(input.observations)
                    },
                ),
            )
        }
        if (model.isNotEmpty()) add(line("model", model))
        if (kinship.isNotEmpty()) add(line("kinship", kinship))
        if (seed != null) add(line("seed", seed.toString()))
        if (permutations != null) {
            add(line("permutations", permutations.permutations.toString()))
            for ((trait, levels) in permutations.thresholds) {
                for ((alpha, threshold) in levels) {
                    add(
                        line("permutation_threshold.$trait.$alpha", threshold.toString()),
                    )
                }
            }
        }
        for ((reason, count) in excludedObservations) add(line("excluded_observations.$reason", count.toString()))
        for ((reason, count) in excludedMarkers) add(line("excluded_markers.$reason", count.toString()))
    }

    private fun line(key: String, value: String): String =
        "##$key=" + value.replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
}
