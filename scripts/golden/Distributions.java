import cern.jet.stat.Gamma;
import cern.jet.stat.Probability;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.distribution.TDistribution;

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
 * <p>It also writes {@code golden/stats/distributions.txt}: Colt's F, Student's t, and chi-square cdf and survival
 * values, in the same {@code function a b x value} layout ({@code a} and {@code b} are the degrees of freedom, {@code x}
 * the argument, and {@code b} is {@code -} when there is one). Colt has no F distribution, so the F rows are the
 * incomplete beta calls that TASSEL 5 and {@code Probability.beta} make: {@code fCdf} is
 * {@code beta(df1 / 2, df2 / 2, df1 x / (df2 + df1 x))} and {@code fSf} is the TASSEL 5 F-test p-value
 * {@code beta(df2 / 2, df1 / 2, df2 / (df2 + df1 x))}. {@code tCdf} is {@code Probability.studentT(df, t)} and
 * {@code tSf} is {@code studentT(df, -t)}, the accurate side of Colt's symmetric formula. {@code chiSquareCdf} and
 * {@code chiSquareSf} are {@code chiSquare} and {@code chiSquareComplemented}. A row is written only when its value is
 * in [1e-300, 1]: below that the double is subnormal or zero and no port can match it.
 *
 * <p>And it writes {@code golden/stats/distribution-quantiles.txt}: commons-math3 3.4.1 {@code
 * inverseCumulativeProbability}, a Brent solve to 1e-9 absolute, as {@code fQuantile}, {@code tQuantile}, and
 * {@code chiSquareQuantile} rows with {@code x} the probability. Quantiles above {@link #MAX_QUANTILE} are left out
 * because 1e-9 absolute is finer than commons-math resolves them.
 *
 * <p>Usage: {@code java -cp colt-1.2.0.jar:commons-math3-3.4.1.jar:. Distributions SPECIAL DISTRIBUTIONS QUANTILES}
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

    private static final double[] CHI_DF = {1, 1.5, 2, 3, 5, 10, 30, 100, 500, 1000, 1e4, 1e6};
    // Multiples of the degrees of freedom, so every df reaches both tails; CHI_X_ABS adds fixed points.
    private static final double[] CHI_X_SCALE = {
        1e-6, 1e-3, 0.01, 0.1, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 3, 5, 10, 30, 100, 300, 1000
    };
    private static final double[] CHI_X_ABS = {1e-100, 1e-10, 0.001, 0.5, 2, 20, 100, 700, 1400, 1e5};
    private static final double[] T_DF = {1, 1.5, 2, 3, 5, 10, 30, 100, 1000, 1e5, 1e7};
    private static final double[] T_X = {
        0, 0.1, 0.5, 1, 1.5, 2, 3, 5, 10, 30, 100, 1e3, 1e5, 1e10, 1e100, 1e150
    };
    private static final double[] F_DF1 = {1, 2, 3, 5, 10, 30, 100, 1000};
    private static final double[] F_DF2 = {1, 2, 3, 5, 10, 30, 100, 1000, 1e5};
    private static final double[] F_X = {1e-8, 1e-4, 0.01, 0.1, 0.5, 1, 2, 5, 10, 30, 100, 1e3, 1e5, 1e8, 1e20, 1e100};
    private static final double[] QUANTILE_DF = {1, 1.5, 2, 3, 5, 10, 30, 100, 1000};
    private static final double[] QUANTILE_P = {1e-6, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999};
    private static final double MAX_QUANTILE = 1e6;

    public static void main(String[] args) throws IOException {
        if (args.length != 3) {
            throw new IllegalArgumentException("usage: Distributions SPECIAL DISTRIBUTIONS QUANTILES");
        }
        writeDistributions(Path.of(args[1]));
        writeQuantiles(Path.of(args[2]));
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

    private static void writeDistributions(Path out) throws IOException {
        Files.createDirectories(out.toAbsolutePath().getParent());
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(out, StandardCharsets.UTF_8))) {
            w.print("function\ta\tb\tx\tvalue\n");
            for (double df : CHI_DF) {
                java.util.TreeSet<Double> xs = new java.util.TreeSet<>();
                for (double m : CHI_X_SCALE) {
                    xs.add(df * m);
                }
                for (double x : CHI_X_ABS) {
                    xs.add(x);
                }
                for (double x : xs) {
                    tail(w, "chiSquareCdf", df, null, x, Probability.chiSquare(df, x));
                    tail(w, "chiSquareSf", df, null, x, Probability.chiSquareComplemented(df, x));
                }
            }
            for (double df : T_DF) {
                for (double t : T_X) {
                    // Colt forms 1 - incompleteBeta(df / 2, 1/2, df / (df + t^2)). For |t| < 1 and a large df the
                    // argument rounds to within 1e-16 * df / t^2 of 1 and the value loses that much relative
                    // accuracy, more than 1e-13 once df / t^2 passes 1e3. Those rows are Colt's error, not an
                    // oracle, so they are left out.
                    if (t > 0 && t < 1 && df / (t * t) > 1e3) {
                        continue;
                    }
                    for (double signed : t == 0 ? new double[] {0} : new double[] {t, -t}) {
                        tail(w, "tCdf", df, null, signed, Probability.studentT(df, signed));
                        tail(w, "tSf", df, null, signed, Probability.studentT(df, -signed));
                    }
                }
            }
            for (double df1 : F_DF1) {
                for (double df2 : F_DF2) {
                    for (double x : F_X) {
                        double dx = df1 * x;
                        tail(w, "fCdf", df1, df2, x, Probability.beta(df1 / 2, df2 / 2, dx / (df2 + dx)));
                        tail(w, "fSf", df1, df2, x, Probability.beta(df2 / 2, df1 / 2, df2 / (df2 + dx)));
                    }
                }
            }
        }
    }

    /** One cdf or survival row, left out when the value is outside [1e-300, 1] (zero, subnormal, or above one). */
    private static void tail(PrintWriter w, String function, double a, Double b, double x, double v) {
        if (v >= 1e-300 && v <= 1.0) {
            w.print(function + "\t" + cell(a) + "\t" + cell(b) + "\t" + cell(x) + "\t" + Double.toString(v) + "\n");
        }
    }

    private static void writeQuantiles(Path out) throws IOException {
        Files.createDirectories(out.toAbsolutePath().getParent());
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(out, StandardCharsets.UTF_8))) {
            w.print("function\ta\tb\tx\tvalue\n");
            for (double df : QUANTILE_DF) {
                for (double p : QUANTILE_P) {
                    quantile(w, "chiSquareQuantile", df, null, p, new ChiSquaredDistribution(df));
                    quantile(w, "tQuantile", df, null, p, new TDistribution(df));
                    for (double df2 : QUANTILE_DF) {
                        quantile(w, "fQuantile", df, df2, p, new FDistribution(df, df2));
                    }
                }
            }
        }
    }

    private static void quantile(PrintWriter w, String function, double a, Double b, double p,
                                 AbstractRealDistribution d) {
        double v = d.inverseCumulativeProbability(p);
        // commons-math stops its Brent solve at 1e-9 absolute, so its answer can be a few 1e-9 off the true quantile.
        // Keep a row only when the quantile is within 1e-10 by commons-math's own cdf and density.
        double error = Math.abs(d.cumulativeProbability(v) - p) / d.density(v);
        if (Math.abs(v) <= MAX_QUANTILE && error <= 1e-10) {
            w.print(function + "\t" + cell(a) + "\t" + cell(b) + "\t" + cell(p) + "\t" + Double.toString(v) + "\n");
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
