package net.maizegenetics.tassel.association

import net.maizegenetics.tassel.dataset.DataSet
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.genotype.Nucleotide
import net.maizegenetics.tassel.phenotype.AttributeType
import net.maizegenetics.tassel.phenotype.FactorAttribute
import net.maizegenetics.tassel.phenotype.NumericAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeAttribute
import net.maizegenetics.tassel.phenotype.PhenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList
import kotlin.math.PI
import kotlin.math.cos
import kotlin.math.ln
import kotlin.math.sqrt
import kotlin.random.Random

/**
 * [taxa] x [markers] with hets, multiallelic and monomorphic sites, and missing calls: random ones at [missingRate],
 * none at every seventh marker, one shared pattern (so a pattern base is reused), and exactly the members of a rare
 * factor level (so dropping them removes the level and the kept rows' base loses rank). Three traits: `y1` and `y2`
 * are complete and share an observation set, `y3` misses some values; two covariates, one with a missing value, and a
 * factor.
 */
internal fun randomGlmDesign(seed: Int, taxa: Int = 90, markers: Int = 70, missingRate: Double = 0.06): DataSet {
    val random = Random(seed)
    val names = (0 until taxa).map { "t$it" }
    val rare = setOf(taxa - 2, taxa - 1)
    val shared = (0 until taxa).filter { random.nextDouble() < 0.25 }.toSet() - rare
    val nucleotides = listOf('A', 'C', 'G', 'T')
    val calls = Array(markers) { m ->
        val (a, b, c) = nucleotides.shuffled(random)
        val frequency = 0.1 + 0.8 * random.nextDouble()
        Array(taxa) { t ->
            val missing = when (m % 7) {
                0 -> false
                3 -> t in shared
                4 -> t in rare
                else -> random.nextDouble() < missingRate * (1 + m % 3)
            }
            fun allele() = if (random.nextDouble() < frequency) a else b
            when {
                missing -> "NN"
                m % 11 == 5 -> "$a$a"
                m % 9 == 2 -> listOf(a, b, c).let { "${it[random.nextInt(3)]}${it[random.nextInt(3)]}" }
                else -> "${allele()}${allele()}"
            }
        }
    }
    val genotypes = InMemoryGenotypeTable.of(
        TaxaList.of(names),
        PositionList((0 until markers).map { Position(Chromosome.of(1 + it / 40), (it + 1) * 100, "m$it") }),
    ) { taxon, site -> Nucleotide.genotype(calls[site][taxon]) }

    val effect = DoubleArray(taxa) { t ->
        calls[1][t].count { it == calls[1][0][0] } * 1.7 + calls[8][t].count { it == calls[8][0][0] } * 0.4
    }
    val c1 = FloatArray(taxa) { random.nextDouble(-2.0, 2.0).toFloat() }
    val c2 = FloatArray(taxa) { if (it == 7) Float.NaN else random.nextDouble(0.0, 10.0).toFloat() }
    val factorCodes = IntArray(taxa) { if (it in rare) 3 else random.nextInt(3) }
    fun trait(name: String, missing: Double, noise: Double) = NumericAttribute(
        name,
        AttributeType.data,
        FloatArray(taxa) {
            if (random.nextDouble() < missing) {
                Float.NaN
            } else {
                (100 + effect[it] + 0.8 * c1[it] + factorCodes[it] + noise * random.nextGaussian()).toFloat()
            }
        },
    )
    val attributes = listOf<PhenotypeAttribute>(
        PhenotypeTable.taxa(names),
        trait("y1", 0.0, 1.0),
        trait("y2", 0.0, 3.0),
        trait("y3", 0.1, 0.5),
        NumericAttribute("c1", AttributeType.covariate, c1),
        NumericAttribute("c2", AttributeType.covariate, c2),
        FactorAttribute("f", factorCodes, listOf("a", "b", "c", "z")),
    )
    return DataSet.of(genotypes, PhenotypeTable(attributes))
}

internal fun Random.nextGaussian(): Double {
    // Box-Muller; one draw is enough here.
    val u = 1.0 - nextDouble()
    val v = nextDouble()
    return sqrt(-2.0 * ln(u)) * cos(2.0 * PI * v)
}
