package net.maizegenetics.util;

/**
 * Drop-in replacement for TASSEL 5's {@code net.maizegenetics.util.DoubleFormat}, placed ahead of
 * {@code sTASSEL.jar} on the classpath by {@code scripts/generate-goldens.sh} to write the
 * full-precision ({@code *.full.txt}) goldens. TASSEL's version rounds to five significant digits;
 * this one writes {@link Double#toString(double)}. NaN and infinity are written as TASSEL writes them.
 */
public class DoubleFormat {

    private DoubleFormat() {
    }

    public static String format(double value) {
        if (Double.isNaN(value)) {
            return "NaN";
        }
        if (Double.isInfinite(value)) {
            return "Infinity";
        }
        return Double.toString(value);
    }

    public static String format(Double value) {
        if (value == null) {
            return "NaN";
        }
        return format(value.doubleValue());
    }
}
