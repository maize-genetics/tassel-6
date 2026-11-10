package net.maizegenetics.tassel.spike.gram

import net.maizegenetics.tassel.distance.DEFAULT_MAX_ALLELES
import net.maizegenetics.tassel.genotype.Allele
import net.maizegenetics.tassel.genotype.AlleleStats
import net.maizegenetics.tassel.genotype.Genotype
import net.maizegenetics.tassel.genotype.GenotypeTable
import net.maizegenetics.tassel.genotype.InMemoryGenotypeTable
import net.maizegenetics.tassel.position.Chromosome
import net.maizegenetics.tassel.position.Position
import net.maizegenetics.tassel.position.PositionList
import net.maizegenetics.tassel.taxa.TaxaList

/**
 * Centered dosage `copies - twiceP`, or `0` when the genotype is not comparable at [allele].
 *
 * Same rule as `distance.centeredDosage`: allele codes are read in their low three bits, two
 * matching alleles contribute `2`, one contributes `1`, and a comparable genotype with no copies
 * contributes `0` before centering. Fully missing genotypes contribute `0` after centering.
 */
fun centeredDosageValue(genotype: Int, allele: Int, twiceP: Double): Double {
    val alleleBits = allele and 0x7
    val first = (genotype shr 4) and 0x7
    val second = genotype and 0x7
    val unknown = 0x7
    val copies = when {
        alleleBits == unknown || (first == unknown && second == unknown) -> -1
        first == alleleBits && second == alleleBits -> 2
        first == alleleBits || second == alleleBits -> 1
        else -> 0
    }
    return when (copies) {
        0 -> 0.0 - twiceP
        1 -> 1.0 - twiceP
        2 -> 2.0 - twiceP
        else -> 0.0
    }
}

/** Diploid `A/C` encoding of a dosage code. `0` is `C/C`, `1` is `A/C`, `2` is `A/A`, missing is `N/N`. */
fun encodeGenotype(code: Int): Byte = when (code) {
    0 -> Genotype.of(Allele.C, Allele.C)
    1 -> Genotype.of(Allele.A, Allele.C)
    2 -> Genotype.of(Allele.A, Allele.A)
    else -> Genotype.UNKNOWN
}

fun encodeGenotypes(codes: ByteArray): ByteArray = ByteArray(codes.size) { encodeGenotype(codes[it].toInt()) }

fun genotypeTable(codes: ByteArray, taxa: Int, sites: Int): GenotypeTable {
    val chromosome = Chromosome.of(1)
    val positions = PositionList(List(sites) { Position(chromosome, (it + 1) * 1_000) })
    val names = List(taxa) { "t${it + 1}" }
    return InMemoryGenotypeTable.of(TaxaList.of(names), positions, encodeGenotypes(codes))
}

class CenteredGram(val values: DoubleArray, val sumPk: Double)

/**
 * Walks centered-IBS columns in site order.
 *
 * A site contributes its major allele column when it has two or more alleles, matching
 * `centeredIbs` at the default `maxAlleles` (a biallelic site evaluates one allele). Frequency is
 * the Float division `centeredIbs` uses. [sumPk] is `2 * sum(p * (1 - p))` and is complete only
 * after [nextPanel] has returned `0`.
 */
class CenteredWalker(private val genotypes: ByteArray, private val taxa: Int, private val sites: Int) {
    private var site: Int = 0
    private val siteGenotypes = ByteArray(taxa)
    private var sumPi: Double = 0.0

    val sumPk: Double get() = sumPi * 2.0

    fun nextPanel(panel: DoubleArray, width: Int): Int {
        var filled = 0
        while (site < sites && filled < width) {
            val offset = site * taxa
            for (taxon in 0 until taxa) siteGenotypes[taxon] = genotypes[offset + taxon]
            val stats = AlleleStats.of(siteGenotypes)
            val evaluated = minOf(stats.alleles.size - 1, DEFAULT_MAX_ALLELES - 1)
            if (evaluated > width) error("panel width $width cannot hold $evaluated centered columns")
            if (evaluated > 0 && filled + evaluated > width) break
            if (evaluated > 0) {
                val gametes = stats.alleleGameteCount.toFloat()
                for (alleleIndex in 0 until evaluated) {
                    val frequency = stats.alleleCounts[alleleIndex].toFloat() / gametes
                    sumPi += frequency * (1.0 - frequency)
                    val twiceP = frequency.toDouble() * 2.0
                    val allele = stats.alleles[alleleIndex].toInt()
                    val column = filled
                    for (taxon in 0 until taxa) {
                        panel[taxon * width + column] =
                            centeredDosageValue(siteGenotypes[taxon].toInt(), allele, twiceP)
                    }
                    filled++
                }
            }
            site++
        }
        return filled
    }
}

/** f64 centered Gram at [panelColumns], plus `SumPk`. */
fun centeredGram(codes: ByteArray, taxa: Int, sites: Int, panelColumns: Int = PANEL_COLUMNS): CenteredGram {
    val walker = CenteredWalker(encodeGenotypes(codes), taxa, sites)
    val triangle = DoubleArray(pairCount(taxa))
    val panel = DoubleArray(taxa * panelColumns)
    while (true) {
        val columns = walker.nextPanel(panel, panelColumns)
        if (columns == 0) break
        addPanel(panel, taxa, columns, panelColumns, triangle)
    }
    return CenteredGram(triangle, walker.sumPk)
}

/** One tight f32 panel of centered columns, for upload. Returns the live column count. */
fun nextCenteredFloatPanel(walker: CenteredWalker, taxa: Int, width: Int): FloatArray? {
    val doubles = DoubleArray(taxa * width)
    val columns = walker.nextPanel(doubles, width)
    if (columns == 0) return null
    val panel = FloatArray(taxa * columns)
    for (taxon in 0 until taxa) {
        val source = taxon * width
        val dest = taxon * columns
        for (column in 0 until columns) {
            panel[dest + column] = doubles[source + column].toFloat()
        }
    }
    return panel
}
