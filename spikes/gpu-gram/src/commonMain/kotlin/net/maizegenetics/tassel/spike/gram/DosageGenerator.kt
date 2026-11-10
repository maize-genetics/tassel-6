package net.maizegenetics.tassel.spike.gram

/**
 * Park-Miller LCG already used by `scripts/bench/generate-wheat-vcf.sh`.
 *
 * `state = (state * 48271) % 2147483647`, kept in `1 .. 2147483646`. The product is taken in
 * 64-bit integer arithmetic, which matches the awk script: every intermediate value fits in the
 * exact integer range of an IEEE double, and in a signed 64-bit integer.
 */
class ParkMiller(seed: Int = PARK_MILLER_SEED) {
    private var state: Int = seed

    init {
        require(state in 1 until PARK_MILLER_MODULUS) { "seed must be in 1..2147483646: $seed" }
    }

    fun next(): Int {
        state = ((state.toLong() * PARK_MILLER_MULTIPLIER) % PARK_MILLER_MODULUS).toInt()
        return state
    }
}

/**
 * Site-major dosage codes for [taxa] by [sites].
 *
 * For each site, one draw sets the minor-allele percent to `5 + (next() % 46)` (an integer in
 * `5..50`). For each taxon at that site, the next draw is missing when `next() % 100 < 5`;
 * otherwise two further draws are Bernoulli trials at that percent and the code is their sum
 * (`0`, `1`, or `2` copies). Missing is stored as [MISSING_CODE].
 *
 * No floating-point value is used. The same arguments always produce the same bytes.
 */
fun generateDosages(taxa: Int, sites: Int, seed: Int = PARK_MILLER_SEED): ByteArray {
    require(taxa > 0 && sites > 0) { "taxa and sites must be positive: $taxa x $sites" }
    val random = ParkMiller(seed)
    val codes = ByteArray(sites * taxa)
    var index = 0
    for (site in 0 until sites) {
        val mafPercent = MAF_MIN_PERCENT + (random.next() % (MAF_MAX_PERCENT - MAF_MIN_PERCENT + 1))
        for (taxon in 0 until taxa) {
            val code = if (random.next() % 100 < MISSING_PERCENT) {
                MISSING_CODE
            } else {
                val first = if (random.next() % 100 < mafPercent) 1 else 0
                val second = if (random.next() % 100 < mafPercent) 1 else 0
                first + second
            }
            codes[index++] = code.toByte()
        }
    }
    return codes
}

/**
 * Wrapping unsigned 64-bit checksum: start at `0`, then `checksum = checksum * 131 + code` over
 * the site-major codes. Codes are `0`, `1`, `2`, or `3`.
 */
fun dosageChecksum(codes: ByteArray): ULong {
    var checksum = 0uL
    for (code in codes) {
        checksum = checksum * 131uL + code.toUByte().toULong()
    }
    return checksum
}
