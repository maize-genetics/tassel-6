import cern.jet.stat.Gamma;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Writes the Colt 1.2.0 {@code cern.jet.stat.Gamma} values that {@code scripts/generate-goldens.sh} records as
 * {@code golden/stats/special-functions.txt}. Every TASSEL 5 GLM and MLM p-value reaches Colt through
 * {@code LinearModelUtils.Ftest}, which is {@code Gamma.incompleteBeta(df2 / 2, df1 / 2, df2 / (df2 + df1 * F))}.
 *
 * <p>The file is tab-delimited with the header {@code function a b x value}. {@code incompleteBeta} rows use all three
 * arguments; {@code incompleteGamma} and {@code incompleteGammaComplement} use {@code a} and {@code x}; {@code gamma}
 * and {@code logGamma} use only {@code a}. An unused column is {@code -}. Arguments and values are written with
 * {@link Double#toString(double)}, and a call that throws is written as {@code error}. Rows whose value is
 * ill-conditioned are left out; see {@link #row}.
 *
 * <p>Usage: {@code java -cp colt-1.2.0.jar:. Distributions OUTPUT}
 */
public class Distributions {

    private static final double[] BETA_SHAPES = {
        1e-8, 0.001, 0.5, 1, 1.5, 2, 3.5, 10, 50, 100, 500, 1000, 1e5, 1e7
    };
    private static final double[] BETA_X = {
        1e-300, 1e-100, 1e-5, 1e-3, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999999, 1 - 1e-16, 1.0
    };
    // a, b, and x outside the domain; x = 0 and x = 1 are in BETA_X.
    private static final double[][] BETA_DOMAIN = {
        {0, 1, 0.5}, {-1, 1, 0.5}, {1, 0, 0.5}, {1, -2.5, 0.5}, {1, 1, -0.1}, {1, 1, 1.1}, {0, 0, 2}, {2, 3, 0}, {2, 3, 1}
    };
    // The F-test arguments df2 / 2, df1 / 2, df2 / (df2 + df1 * F) for ordinary and non-integer degrees of freedom.
    private static final double[][] F_TEST = {
        {1, 1, 1}, {1, 100, 3.9}, {5, 2, 10}, {10.5, 3.25, 2.7}, {278, 3, 0.4}, {2.5, 7.5, 15}, {1000, 12, 1.1},
        {0.5, 0.5, 0.25}, {50, 1, 1e-3}, {1, 40, 5e4}, {3, 3, 1e6}
    };
    private static final double[] GAMMA_SHAPES = {1e-8, 0.001, 0.5, 1, 1.5, 3, 7.5, 20, 100, 1000, 1e5};
    private static final double[] GAMMA_X = {
        1e-10, 0.001, 0.1, 0.5, 1, 2, 5, 10, 50, 100, 500, 1000, 1e5, 1e8
    };
    private static final double[] GAMMA_ARGS = {
        -100.5, -50.5, -35, -34.5, -33.5, -10.25, -3.5, -2, -1.5, -0.5, -1e-9, -1e-10, 0, 1e-10, 1e-9, 1e-5,
        0.001, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4.5, 10, 20.5, 33, 33.5, 34, 50, 100, 143, 143.5, 171, 171.5, 171.62,
        172, 200
    };
    private static final double[] LOG_GAMMA_ARGS = {
        -1e5 - 0.5, -100.5, -50.5, -35, -34.5, -34, -33.5, -10.25, -3.5, -3, -1.5, -0.5, -1e-9, 0, 1e-300, 1e-10,
        1e-5, 0.001, 0.5, 1, 1.5, 2, 2.5, 3, 4.5, 10, 12.5, 12.99, 13, 13.5, 100, 999, 1000, 1001, 1e5, 1e8,
        1e8 + 1, 1e10, 1e100, 1e300, 2.5e305, 2.6e305, Double.POSITIVE_INFINITY
    };

    public static void main(String[] args) throws IOException {
        if (args.length != 1) {
            throw new IllegalArgumentException("usage: Distributions OUTPUT");
        }
        Path out = Path.of(args[0]);
        Files.createDirectories(out.toAbsolutePath().getParent());
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(out, StandardCharsets.UTF_8))) {
            w.print("function\ta\tb\tx\tvalue\n");
            for (double a : BETA_SHAPES) {
                for (double b : BETA_SHAPES) {
                    for (double x : BETA_X) {
                        beta(w, a, b, x);
                    }
                }
            }
            for (double[] row : BETA_DOMAIN) {
                beta(w, row[0], row[1], row[2]);
            }
            for (double[] row : F_TEST) {
                double df1 = row[0] * 2, df2 = row[1] * 2, f = row[2];
                beta(w, df2 / 2, df1 / 2, df2 / (df2 + df1 * f));
            }
            for (double a : GAMMA_SHAPES) {
                for (double x : GAMMA_X) {
                    row(w, "incompleteGamma", a, null, x, (p, q, r) -> Gamma.incompleteGamma(p, r));
                    row(w, "incompleteGammaComplement", a, null, x, (p, q, r) -> Gamma.incompleteGammaComplement(p, r));
                }
            }
            for (double a : new double[] {0, -1}) {
                row(w, "incompleteGamma", a, null, 1.0, (p, q, r) -> Gamma.incompleteGamma(p, r));
                row(w, "incompleteGammaComplement", a, null, 1.0, (p, q, r) -> Gamma.incompleteGammaComplement(p, r));
            }
            for (double x : new double[] {0, -1}) {
                row(w, "incompleteGamma", 2.0, null, x, (p, q, r) -> Gamma.incompleteGamma(p, r));
                row(w, "incompleteGammaComplement", 2.0, null, x, (p, q, r) -> Gamma.incompleteGammaComplement(p, r));
            }
            for (double a : GAMMA_ARGS) {
                row(w, "gamma", a, null, null, (p, q, r) -> Gamma.gamma(p));
            }
            for (double a : LOG_GAMMA_ARGS) {
                row(w, "logGamma", a, null, null, (p, q, r) -> Gamma.logGamma(p));
            }
        }
    }

    private static void beta(PrintWriter w, double a, double b, double x) {
        row(w, "incompleteBeta", a, b, x, Gamma::incompleteBeta);
    }

    /** A Colt call on up to three arguments; the ones a function does not use are ignored. */
    private interface Call {
        double apply(double a, double b, double x);
    }

    /**
     * Writes one row, unless Colt's value there is ill-conditioned. Colt's results move by more than a double's
     * rounding with the platform's {@code exp}, {@code log}, and {@code pow} (they differ in the last bit between
     * the JVM on macOS and Linux, glibc, and V8), so no port can hold {@code rel 1e-13} on those rows. A row is
     * dropped when nudging any argument by two ulps either way moves the value by more than {@link #STABLE}
     * relative, or turns a value into an error. Rows that throw, and x = 0 and x = 1 for the beta function, are
     * kept as they are: they are exact.
     */
    private static void row(PrintWriter w, String function, Double a, Double b, Double x, Call call) {
        String value;
        try {
            double v = call.apply(orZero(a), orZero(b), orZero(x));
            if (!(x != null && function.equals("incompleteBeta") && (x == 0.0 || x == 1.0)) && !stable(v, a, b, x, call)) {
                return;
            }
            value = Double.toString(v);
        } catch (RuntimeException e) {
            value = "error";
        }
        w.print(function + "\t" + cell(a) + "\t" + cell(b) + "\t" + cell(x) + "\t" + value + "\n");
    }

    private static final double STABLE = 2e-14;

    private static boolean stable(double v, Double a, Double b, Double x, Call call) {
        for (int arg = 0; arg < 3; arg++) {
            Double[] in = {a, b, x};
            if (in[arg] == null) {
                continue;
            }
            for (int dir = -1; dir <= 1; dir += 2) {
                double n = in[arg];
                for (int k = 0; k < 2; k++) {
                    n = dir > 0 ? Math.nextUp(n) : Math.nextDown(n);
                }
                Double[] nudged = in.clone();
                nudged[arg] = n;
                double u;
                try {
                    u = call.apply(orZero(nudged[0]), orZero(nudged[1]), orZero(nudged[2]));
                } catch (RuntimeException e) {
                    return false;
                }
                double scale = Math.max(Math.abs(v), Math.abs(u));
                if (u != v && !(Math.abs(u - v) <= STABLE * scale)) {
                    return false;
                }
            }
        }
        return true;
    }

    private static double orZero(Double v) {
        return v == null ? 0.0 : v;
    }

    private static String cell(Double v) {
        return v == null ? "-" : Double.toString(v);
    }
}
