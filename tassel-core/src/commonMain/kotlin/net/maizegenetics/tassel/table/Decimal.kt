package net.maizegenetics.tassel.table

/**
 * A non-negative whole number, in base 10^9 limbs, little endian: just enough to expand a `Double`
 * or `Float` exactly. [DoubleFormat] rounds the expansion; the phenotype reader compares a decimal
 * token against it.
 */
internal class Decimal(value: Long) {
    private var limbs = IntArray(4)
    private var size = 0

    init {
        var rest = value
        do {
            append((rest % BASE).toInt())
            rest /= BASE
        } while (rest > 0)
    }

    /** Multiplies by [base] (2 or 5) raised to [power]. */
    fun multiplyByPowerOf(base: Int, power: Int) {
        val chunk = if (base == 2) 30 else 13 // the largest power of each that fits an Int
        var remaining = power
        while (remaining > 0) {
            val step = minOf(chunk, remaining)
            var factor = 1
            repeat(step) { factor *= base }
            multiply(factor)
            remaining -= step
        }
    }

    /** The digits, with no leading zeros. */
    override fun toString(): String {
        val text = StringBuilder(size * LIMB_DIGITS)
        text.append(limbs[size - 1])
        for (index in size - 2 downTo 0) {
            val limb = limbs[index].toString()
            repeat(LIMB_DIGITS - limb.length) { text.append('0') }
            text.append(limb)
        }
        return text.toString()
    }

    private fun multiply(factor: Int) {
        var carry = 0L
        for (index in 0 until size) {
            // A limb is under 10^9 and a factor under 2^31, so this stays well inside a Long.
            val product = limbs[index].toLong() * factor + carry
            limbs[index] = (product % BASE).toInt()
            carry = product / BASE
        }
        while (carry > 0) {
            append((carry % BASE).toInt())
            carry /= BASE
        }
    }

    private fun append(limb: Int) {
        if (size == limbs.size) limbs = limbs.copyOf(limbs.size * 2)
        limbs[size++] = limb
    }

    private companion object {
        const val BASE = 1_000_000_000L
        const val LIMB_DIGITS = 9
    }
}
