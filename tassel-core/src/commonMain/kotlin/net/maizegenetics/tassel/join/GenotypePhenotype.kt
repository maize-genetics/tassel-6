package net.maizegenetics.tassel.join

import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.ops.operation
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.phenotype.TaxaAttribute

/**
 * A phenotype aligned with a genotype table by taxon (ADR-0024). Phenotype rows drive it: [genotypeRow] gives, for each
 * observation, the taxon index in [genotypes] it was measured on, or `-1` when that taxon has no genotype (a union
 * only).
 *
 * The genotypes are held by reference and never copied. Read them through
 * `genotypes.siteBlocks(taxa = genotypeTaxonIndices())`, where block taxon `i` is then the `i`-th distinct genotyped
 * taxon in observation order. Replicates share one genotype row, so they are read once.
 *
 * @see net.maizegenetics.phenotype.GenotypePhenotype
 */
class GenotypePhenotype(val genotypes: GenotypeTable, val phenotype: PhenotypeTable, genotypeRow: IntArray) {
    private val rows: IntArray = genotypeRow.copyOf()

    init {
        require(rows.size == phenotype.numberOfObservations) {
            "genotypeRow has ${rows.size} entries for ${phenotype.numberOfObservations} observations"
        }
        rows.forEachIndexed { observation, row ->
            require(row in -1 until genotypes.numberOfTaxa) {
                "Observation $observation maps to genotype taxon $row, outside -1 until ${genotypes.numberOfTaxa}"
            }
            require(row < 0 || genotypes.taxa[row].name == phenotype.taxon(observation)) {
                "Observation $observation is taxon '${phenotype.taxon(observation)}' but genotype taxon $row is " +
                    "'${genotypes.taxa[row].name}'"
            }
        }
    }

    /** Number of observations (phenotype rows), with or without a genotype. */
    val numberOfObservations: Int get() = rows.size

    /** The genotype taxon index of [observation], or `-1` when it has none. */
    fun genotypeRow(observation: Int): Int = rows[observation]

    /** The genotype taxon index of every observation, `-1` for none. The caller owns the returned array. */
    fun genotypeRows(): IntArray = rows.copyOf()

    /** The distinct genotype taxon indices, in the order their observations first appear. */
    fun genotypeTaxonIndices(): IntArray {
        val seen = HashSet<Int>()
        return rows.filter { it >= 0 && seen.add(it) }.toIntArray()
    }

    /** The observations that have a genotype, in observation order. */
    fun observationsWithGenotype(): IntArray = rows.indices.filter { rows[it] >= 0 }.toIntArray()

    override fun toString(): String =
        "GenotypePhenotype($numberOfObservations observations, ${genotypeTaxonIndices().size} genotyped taxa)"
}

/**
 * Aligns [phenotype] with [genotypes] by exact, case-sensitive taxon name, as ADR-0024 defines.
 *
 * [JoinType.intersect] keeps the observations whose taxon is genotyped. [JoinType.union] keeps every observation and
 * appends one all-missing observation for each genotype-only taxon, in genotype order. Observation order is otherwise
 * the phenotype's. The genotypes are not read.
 *
 * @see net.maizegenetics.phenotype.GenotypePhenotypeBuilder
 */
fun alignGenotypePhenotype(
    genotypes: GenotypeTable,
    phenotype: PhenotypeTable,
    type: JoinType = JoinType.intersect,
): GenotypePhenotype {
    val taxa = genotypes.taxa
    val found = IntArray(phenotype.numberOfObservations) { taxa.indexOf(phenotype.taxon(it)) }
    if (type == JoinType.intersect) {
        val kept = found.indices.filter { found[it] >= 0 }.toIntArray()
        return GenotypePhenotype(genotypes, phenotype.selectRows(kept), IntArray(kept.size) { found[kept[it]] })
    }
    val present = BooleanArray(taxa.size).also { seen -> found.forEach { if (it >= 0) seen[it] = true } }
    val extra = taxa.indices.filter { !present[it] }
    if (extra.isEmpty()) return GenotypePhenotype(genotypes, phenotype, found)
    return GenotypePhenotype(
        genotypes,
        PhenotypeTable(phenotype.attributes.map { it.extended(extra.map { row -> taxa[row].name }) }, phenotype.name),
        found + IntArray(extra.size) { extra[it] },
    )
}

/** This column with one missing cell per name in [extraTaxa] appended; the taxa column gets the names. */
private fun PhenotypeAttribute.extended(extraTaxa: List<String>): PhenotypeAttribute = when (this) {
    is TaxaAttribute -> TaxaAttribute(name, names + extraTaxa)
    is FactorAttribute -> FactorAttribute.of(name, (0 until size).map { level(it) } + extraTaxa.map { null })
    is NumericAttribute -> NumericAttribute(name, type, toFloatArray() + FloatArray(extraTaxa.size) { Float.NaN })
}

/**
 * The spec for [alignGenotypePhenotype]. Declared here; `TasselOperations` registers it (#166).
 *
 * The input is a [DataSet] with exactly one genotype table and at least one phenotype table, in any order. Several
 * phenotypes are first joined with [joinPhenotypes] using the same `type` and the default collision rules. There is
 * no TASSEL 5 plugin to translate, so no legacy name.
 */
val joinGenotypePhenotypeOperation = operation<DataSet, GenotypePhenotype>(
    id = "join-genotype-phenotype",
    help = "Aligns a genotype table with one or more phenotype tables by taxon, keeping matched or all observations",
) {
    val type by enum(
        JoinType.intersect,
        help = "intersect keeps genotyped observations; union also adds genotype-only taxa",
    )

    run { input, progress ->
        val genotypes = input.genotypes()
        val phenotypes = input.phenotypes()
        require(genotypes.size == 1 && phenotypes.isNotEmpty() && genotypes.size + phenotypes.size == input.size) {
            "join-genotype-phenotype needs exactly one genotype table and at least one phenotype table, and " +
                "nothing else; got ${genotypes.size} genotype tables and ${phenotypes.size} phenotype tables " +
                "among ${input.size} items"
        }
        val phenotype = if (phenotypes.size == 1) {
            phenotypes[0]
        } else {
            joinPhenotypes(phenotypes, JoinOptions(type = type()), progress)
        }
        alignGenotypePhenotype(genotypes[0], phenotype, type())
    }
}
