package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.join.GenotypePhenotype
import net.maizegenetics.tassel.join.JoinOptions
import net.maizegenetics.tassel.join.JoinType
import net.maizegenetics.tassel.join.alignGenotypePhenotype
import net.maizegenetics.tassel.join.joinPhenotypes
import net.maizegenetics.tassel.ops.ProgressReporter
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.stats.model.Contrast
import net.maizegenetics.tassel.stats.model.ModelSpec
import net.maizegenetics.tassel.stats.model.ModelTerm

/**
 * Which attributes of the phenotype an association analysis uses. A null list takes TASSEL 5's default: every `data`
 * attribute is a trait, every `covariate` attribute and every `factor` attribute is a fixed effect
 * (`AbstractFixedEffectLM` builds the model from `attributeListOfType`).
 *
 * A named list may also pick a numeric attribute of any type as a covariate, or a factor attribute as a factor. An empty
 * list means none.
 */
internal data class AssociationSpec(
    val traits: List<String>? = null,
    val covariates: List<String>? = null,
    val factors: List<String>? = null,
)

/**
 * What a GLM or MLM run needs from its inputs (ADR-0026): the genotypes aligned with the phenotype by taxon name, the
 * traits and fixed effects, and for each trait the base observation set, which is the observations where the trait,
 * every covariate, and every factor are present (and that have a genotype row).
 *
 * Genotypes are never copied. A marker's own missing calls are not handled here: they drop observations from that
 * marker's test only, in the [MarkerSource].
 *
 * @property data The aligned genotypes and phenotype; observations are its phenotype rows.
 * @property traits Trait names, in the order results follow.
 * @property covariates Covariate names, in model order.
 * @property factors Factor names, in model order.
 * @property excludedObservations For each trait, the observations left out of its base set, by reason: `trait`,
 * `covariate`, `factor` (the first missing value found, in that order), or `genotype` (no genotype row).
 * @see net.maizegenetics.analysis.association.AbstractFixedEffectLM
 */
internal class AssociationData private constructor(
    val data: GenotypePhenotype,
    val traits: List<String>,
    val covariates: List<String>,
    val factors: List<String>,
    private val sets: List<IntArray>,
    val excludedObservations: List<Map<String, Int>>,
) {
    /** The base observation set of each trait: observation indices, ascending. The arrays are shared; do not change them. */
    val observationSets: List<IntArray> get() = sets

    /** Base observations of the trait at [trait]. */
    fun observations(trait: Int): IntArray = sets[trait]

    /** The taxon of every observation as an id: observations of one taxon (replicates) share it. */
    val taxonIds: IntArray by lazy {
        val ids = HashMap<String, Int>()
        IntArray(data.numberOfObservations) { ids.getOrPut(data.phenotype.taxon(it)) { ids.size } }
    }

    /** Whether some taxon has more than one observation in the base set of [trait] (replicated taxa). */
    fun hasReplicates(trait: Int): Boolean {
        val seen = HashSet<Int>()
        return sets[trait].any { !seen.add(taxonIds[it]) }
    }

    /**
     * The fixed-effect model of [trait] without a marker, in TASSEL 5's order: the mean, then factors (the last level
     * present is the zero reference), then covariates.
     */
    fun baseSpec(trait: Int): ModelSpec = ModelSpec(
        traits[trait],
        buildList {
            add(ModelTerm.Intercept)
            factors.forEach { add(ModelTerm.Factor(it, Contrast.treatment)) }
            covariates.forEach { add(ModelTerm.Covariate(it)) }
        },
    )

    /** Total observations in the aligned phenotype, before any exclusion. */
    val numberOfObservations: Int get() = data.numberOfObservations

    companion object {
        /**
         * Reads [dataSet] and selects traits and fixed effects by [spec].
         *
         * The data set holds one `GenotypePhenotype`, or one genotype table and at least one phenotype table, which
         * are joined and aligned with intersect rules (`join-genotype-phenotype`). Anything else is an error that
         * names what was found.
         *
         * @throws IllegalArgumentException for an unusable data set, an unknown or wrongly typed attribute, an
         *   attribute named twice, or no trait.
         */
        suspend fun of(
            dataSet: DataSet,
            spec: AssociationSpec = AssociationSpec(),
            progress: ProgressReporter = ProgressReporter.None,
        ): AssociationData = of(align(dataSet, progress), spec)

        /** As [of] for data that is already aligned. */
        fun of(data: GenotypePhenotype, spec: AssociationSpec = AssociationSpec()): AssociationData {
            val phenotype = data.phenotype
            fun names(type: AttributeType) = phenotype.attributes(type).map { it.name }
            val traits = spec.traits ?: names(AttributeType.data)
            val covariates = spec.covariates ?: names(AttributeType.covariate)
            val factors = spec.factors ?: names(AttributeType.factor)
            require(traits.isNotEmpty()) { "There is no trait to analyze: the phenotype has no 'data' attribute" }

            val used = HashSet<String>()
            fun take(role: String, name: String) {
                require(phenotype.indexOf(name) >= 0) {
                    "The $role '$name' is not an attribute of '${phenotype.name}'"
                }
                require(used.add("$role:$name")) { "The $role '$name' is listed twice" }
            }
            traits.forEach { name ->
                take("trait", name)
                require(phenotype.attribute(name) is NumericAttribute) {
                    "The trait '$name' must be numeric, not ${phenotype.attribute(name).type}"
                }
            }
            covariates.forEach { name ->
                take("covariate", name)
                require(phenotype.attribute(name) is NumericAttribute) {
                    "The covariate '$name' must be numeric, not ${phenotype.attribute(name).type}"
                }
                require(name !in traits) { "'$name' cannot be both a trait and a covariate" }
            }
            factors.forEach { name ->
                take("factor", name)
                require(phenotype.attribute(name) is FactorAttribute) {
                    "The factor '$name' must be a factor attribute, not ${phenotype.attribute(name).type}"
                }
            }

            val covariateColumns = covariates.map { phenotype.attribute(it) as NumericAttribute }
            val factorColumns = factors.map { phenotype.attribute(it) as FactorAttribute }
            val sets = ArrayList<IntArray>(traits.size)
            val excluded = ArrayList<Map<String, Int>>(traits.size)
            for (name in traits) {
                val trait = phenotype.attribute(name) as NumericAttribute
                val kept = ArrayList<Int>()
                val reasons = LinkedHashMap<String, Int>()
                fun exclude(reason: String) {
                    reasons[reason] = (reasons[reason] ?: 0) + 1
                }
                for (observation in 0 until data.numberOfObservations) {
                    when {
                        trait.isMissing(observation) -> exclude("trait")
                        covariateColumns.any { it.isMissing(observation) } -> exclude("covariate")
                        factorColumns.any { it.isMissing(observation) } -> exclude("factor")
                        data.genotypeRow(observation) < 0 -> exclude("genotype")
                        else -> kept.add(observation)
                    }
                }
                sets.add(kept.toIntArray())
                excluded.add(reasons)
            }
            return AssociationData(data, traits, covariates, factors, sets, excluded)
        }

        private suspend fun align(dataSet: DataSet, progress: ProgressReporter): GenotypePhenotype {
            val pairs = dataSet.genotypePhenotypes()
            val genotypes = dataSet.genotypes()
            val phenotypes = dataSet.phenotypes()
            val others = dataSet.size - pairs.size - genotypes.size - phenotypes.size
            if (pairs.size == 1 && genotypes.isEmpty() && phenotypes.isEmpty() && others == 0) return pairs[0]
            require(pairs.isEmpty() && genotypes.size == 1 && phenotypes.isNotEmpty() && others == 0) {
                "An association analysis needs one genotype-phenotype pairing, or one genotype table and at " +
                    "least one phenotype table, and nothing else; got ${pairs.size} pairings, " +
                    "${genotypes.size} genotype tables, ${phenotypes.size} phenotype tables, and $others other items"
            }
            val phenotype = if (phenotypes.size == 1) {
                phenotypes[0]
            } else {
                joinPhenotypes(phenotypes, JoinOptions(type = JoinType.intersect), progress)
            }
            return alignGenotypePhenotype(genotypes[0], phenotype, JoinType.intersect)
        }
    }
}
