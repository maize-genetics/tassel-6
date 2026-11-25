"""mpmath log-survival grids for F, Student's t, and chi-square (tassel-6#183, ADR-0029).

Usage: tails.py OUT_FILE

Writes a tab-separated table with columns dist, df1, df2, x, log_sf. `x` is a double printed
with repr(), so it round-trips exactly; `log_sf` is the natural log of the survival function
P(X > x) at exactly that double, computed with mpmath at 80 digits and printed to 50
significant digits. Nothing here uses a double-precision distribution function.

Two families of rows:
  * deep tails: for each df, the x at which log p first reaches each target from -1 down to
    -50,000 (located by bisection, then rounded to a double; the row stores the exact log_sf
    at that double, so the target is only a placement). Rows with x above 1e300 are not
    representable and are skipped, which is why heavy-tailed t with small df stops early.
  * the regions scripts/golden/Distributions.java drops: values below 1e-300 (the deep tails
    above), and small |t| with a large df, where Colt's 1 - incompleteBeta form loses digits.
"""
import sys

import mpmath
from mpmath import mp

mp.dps = 80
PRINT_DIGITS = 50

TARGETS = [-1, -3, -10, -30, -100, -300, -1000, -3000, -10000, -50000]
X_MAX = 1e300

CHI_DF = [1, 1.5, 2, 3, 5, 10, 30, 100, 500, 1000, 1e4, 1e6]
T_DF = [1, 1.5, 2, 3, 5, 10, 30, 100, 1000, 1e5, 1e7]
F_DF1 = [1, 2, 3, 5, 10, 30, 100, 1000]
F_DF2 = [1, 2, 3, 5, 10, 30, 100, 1000, 1e5]
# Small |t| with a large df, the Colt-lossy region.
T_SMALL = [0.001, 0.01, 0.1, 0.5, 0.9, 1.0]
T_SMALL_DF = [1e3, 1e5, 1e7]
# Small x, where the survival function is within rounding of 1.
CHI_SMALL_MULT = [1e-6, 1e-3, 1e-1]
F_SMALL = [1e-6, 1e-3, 1e-1]


def chi_log_sf(df, _unused, x):
    return mpmath.log(mpmath.gammainc(mpmath.mpf(df) / 2, mpmath.mpf(x) / 2, mpmath.inf, regularized=True))


def ibeta(a, b, x, xc):
    """Regularized incomplete beta I_x(a, b) by the modified-Lentz continued fraction (Numerical
    Recipes betacf, in arbitrary precision); xc = 1 - x is passed in exactly. It converges fast
    and without cancellation when x < (a + 1) / (a + b + 2). Returns (log I, I)."""
    log_prefix = (
        a * mpmath.log(x) + b * mpmath.log(xc)
        - (mpmath.loggamma(a) + mpmath.loggamma(b) - mpmath.loggamma(a + b))
    )
    eps = mpmath.mpf(10) ** (-(mp.dps + 5))
    tiny = mpmath.mpf(10) ** (-(mp.dps * 4))
    qab, qap, qam = a + b, a + 1, a - 1
    c = mpmath.mpf(1)
    d = 1 - qab * x / qap
    d = 1 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 10**6):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        d = 1 / (d if abs(d) > tiny else tiny)
        c = 1 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        d = 1 / (d if abs(d) > tiny else tiny)
        c = 1 + aa / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1) < eps:
            break
    else:
        raise RuntimeError("continued fraction did not converge")
    log_i = log_prefix + mpmath.log(h / a)
    return log_i, mpmath.exp(log_i)


def lower_or_complement(a, b, z, w):
    """log I_z(a, b) given both z and its exact complement w = 1 - z. Uses whichever side the
    continued fraction converges on, so a small result is never a difference of values near 1."""
    if z < (a + 1) / (a + b + 2):
        return ibeta(a, b, z, w)[0]
    return mpmath.log1p(-ibeta(b, a, w, z)[1])


def t_log_sf(df, _unused, x):
    df, x = mpmath.mpf(df), mpmath.mpf(x)
    z, w = df / (df + x * x), x * x / (df + x * x)
    tail = mpmath.exp(lower_or_complement(df / 2, mpmath.mpf(1) / 2, z, w)) / 2
    return mpmath.log(tail if x >= 0 else 1 - tail)


def f_log_sf(df1, df2, x):
    df1, df2, x = mpmath.mpf(df1), mpmath.mpf(df2), mpmath.mpf(x)
    z, w = df2 / (df2 + df1 * x), df1 * x / (df2 + df1 * x)
    return lower_or_complement(df2 / 2, df1 / 2, z, w)


def solve(log_sf, df1, df2, target):
    """The double x > 0 where log_sf crosses target, or None when it is above X_MAX."""
    saved = mp.dps
    mp.dps = 30
    try:
        lo, hi = -690.0, 690.0  # bisection over ln x
        if log_sf(df1, df2, mpmath.exp(hi)) > target:
            return None
        for _ in range(40):
            mid = (lo + hi) / 2
            if log_sf(df1, df2, mpmath.exp(mid)) > target:
                lo = mid
            else:
                hi = mid
        return float(mpmath.exp((lo + hi) / 2))
    finally:
        mp.dps = saved


def main(out_file):
    rows = []

    def add(dist, log_sf, df1, df2, x):
        value = log_sf(df1, df2, x)
        rows.append((dist, df1, df2, x, mpmath.nstr(value, PRINT_DIGITS, min_fixed=1, max_fixed=0)))

    for df in CHI_DF:
        for target in TARGETS:
            x = solve(chi_log_sf, df, None, target)
            if x is not None and x <= X_MAX:
                add("chisq", chi_log_sf, df, None, x)
        for m in CHI_SMALL_MULT:
            add("chisq", chi_log_sf, df, None, df * m)
    for df in T_DF:
        for target in TARGETS:
            x = solve(t_log_sf, df, None, target)
            if x is not None and x <= X_MAX:
                add("t", t_log_sf, df, None, x)
    for df in T_SMALL_DF:
        for t in T_SMALL:
            add("t", t_log_sf, df, None, t)
            add("t", t_log_sf, df, None, -t)
    for df1 in F_DF1:
        for df2 in F_DF2:
            for target in TARGETS:
                x = solve(f_log_sf, df1, df2, target)
                if x is not None and x <= X_MAX:
                    add("f", f_log_sf, df1, df2, x)
            for x in F_SMALL:
                add("f", f_log_sf, df1, df2, x)

    with open(out_file, "w", encoding="ascii", newline="\n") as out:
        out.write("dist\tdf1\tdf2\tx\tlog_sf\n")
        for dist, df1, df2, x, value in rows:
            out.write(f"{dist}\t{df1!r}\t{'-' if df2 is None else repr(df2)}\t{x!r}\t{value}\n")


if __name__ == "__main__":
    main(sys.argv[1])
