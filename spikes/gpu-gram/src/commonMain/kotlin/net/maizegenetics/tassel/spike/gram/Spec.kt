package net.maizegenetics.tassel.spike.gram

/** Columns of one Gram panel. The f64 reference and the GPU dispatch both use this width. */
const val PANEL_COLUMNS: Int = 1024

/** Spec defaults for the generator. Issue #113 reproduces these numbers. */
const val PARK_MILLER_SEED: Int = 81

const val PARK_MILLER_MULTIPLIER: Int = 48_271

const val PARK_MILLER_MODULUS: Int = 2_147_483_647

const val MISSING_PERCENT: Int = 5

const val MAF_MIN_PERCENT: Int = 5

const val MAF_MAX_PERCENT: Int = 50

/** Dosage byte for a missing genotype. Called dosages are 0, 1, or 2. */
const val MISSING_CODE: Int = 3

/** Absolute tolerance on a scaled Float kinship value (ADR-0020). */
const val SCALED_TOLERANCE: Double = 8.0e-6

/** Denominator floor for relative error, so a zero reference does not explode. */
const val RELATIVE_FLOOR: Double = 1.0e-12

data class GramShape(val taxa: Int, val sites: Int) {
    val label: String = "${taxa}x$sites"
}

val GRAM_SHAPES: List<GramShape> = listOf(
    GramShape(281, 3_093),
    GramShape(1_000, 10_000),
    GramShape(1_000, 100_000),
    GramShape(5_000, 20_000),
)

fun pairCount(taxa: Int): Int {
    require(taxa >= 0) { "taxa must not be negative: $taxa" }
    val count = taxa.toLong() * (taxa + 1L) / 2L
    require(count <= Int.MAX_VALUE) { "triangle for $taxa taxa does not fit in an array" }
    return count.toInt()
}

/**
 * Upper-triangle index for pair `(row, column)` with `column >= row`:
 * `(0,0), (0,1), ... (0,n-1), (1,1), ...`.
 *
 * The WGSL helpers `triangle_start` and `triangle_index` follow this formula.
 */
fun triangleIndex(row: Int, column: Int, taxa: Int): Int {
    require(column >= row) { "column $column is above row $row" }
    return row * (2 * taxa - row + 1) / 2 + (column - row)
}
