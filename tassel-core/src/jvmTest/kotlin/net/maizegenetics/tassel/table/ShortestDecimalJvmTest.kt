package net.maizegenetics.tassel.table

import kotlin.random.Random
import kotlin.test.Test
import kotlin.test.assertEquals

/**
 * [floatToString] and [doubleToString] against the JVM's own `Float.toString` and `Double.toString`,
 * which are what TASSEL 5 wrote the goldens with. On the JVM the platform spelling is already
 * Java's, so this pins the re-layout as an exact identity over the whole value range, including the
 * subnormals that `ShortestDecimalTest` leaves to the platform.
 */
class ShortestDecimalJvmTest {
    @Test
    fun agreesWithJavaOverRandomFloats() {
        val random = Random(20261104)
        repeat(VALUES) {
            val value = Float.fromBits(random.nextInt())
            if (!value.isNaN() && !value.isInfinite()) {
                assertEquals(java.lang.Float.toString(value), floatToString(value))
            }
        }
    }

    @Test
    fun agreesWithJavaOverRandomDoubles() {
        val random = Random(20261104)
        repeat(VALUES) {
            val value = Double.fromBits(random.nextLong())
            if (!value.isNaN() && !value.isInfinite()) {
                assertEquals(java.lang.Double.toString(value), doubleToString(value))
            }
        }
    }

    @Test
    fun agreesWithJavaOnEveryKinshipGoldenValue() {
        var bits = 1
        repeat(VALUES) {
            // Values the size of a kinship coefficient, where the plain and scientific forms meet.
            bits = bits * 1664525 + 1013904223
            val value = (bits % 4_000_000) / 1_000_000.0f
            assertEquals(java.lang.Float.toString(value), floatToString(value))
            assertEquals(java.lang.Double.toString(value.toDouble()), doubleToString(value.toDouble()))
        }
    }

    private companion object {
        const val VALUES = 200_000
    }
}
