# R lm() / anova() / drop1() reference fixtures (tassel-6#183, ADR-0029).
#
# Usage: Rscript lm-fixtures.R OUT_DIR
#
# First a smoke fixture (a covariate and a three-level factor) that proves the pinned R container,
# the file layout, and the test-side loader; then the model designs tassel-6#186 needs (contrasts,
# collinear, missing-row, unbalanced, and weighted cases). Every number is written with %.17g so
# the files are byte-identical across runs of the same image.
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

fmt <- function(x) sprintf("%.17g", x)
write_tsv <- function(frame, name) {
  lines <- c(paste(names(frame), collapse = "\t"),
             apply(frame, 1, function(r) paste(r, collapse = "\t")))
  writeLines(lines, file.path(out, name), sep = "\n", useBytes = TRUE)
}

data <- data.frame(
  y = c(5.1, 6.0, 5.5, 7.9, 8.2, 7.4, 9.8, 10.9, 10.1, 6.3, 8.8, 11.6),
  cov = c(1.0, 2.0, 1.5, 3.0, 3.5, 2.5, 4.0, 5.0, 4.5, 2.2, 3.8, 5.5),
  group = factor(c("A", "A", "A", "B", "B", "B", "C", "C", "C", "A", "B", "C"))
)
write_tsv(data.frame(y = fmt(data$y), cov = fmt(data$cov), group = as.character(data$group)),
          "lm-smoke.data.tsv")

fit <- lm(y ~ cov + group, data = data)
coefficients <- summary(fit)$coefficients
table <- data.frame(
  term = rownames(coefficients),
  estimate = fmt(coefficients[, 1]),
  se = fmt(coefficients[, 2]),
  t = fmt(coefficients[, 3]),
  p = fmt(coefficients[, 4])
)
write_tsv(table, "lm-smoke.coefficients.tsv")

sequential <- anova(fit)
write_tsv(data.frame(
  term = rownames(sequential), df = fmt(sequential[["Df"]]), ss = fmt(sequential[["Sum Sq"]]),
  ms = fmt(sequential[["Mean Sq"]]), f = fmt(sequential[["F value"]]), p = fmt(sequential[["Pr(>F)"]])
), "lm-smoke.anova.tsv")

marginal <- drop1(fit, test = "F")
write_tsv(data.frame(
  term = rownames(marginal), df = fmt(marginal[["Df"]]), ss = fmt(marginal[["Sum of Sq"]]),
  rss = fmt(marginal[["RSS"]]), f = fmt(marginal[["F value"]]), p = fmt(marginal[["Pr(>F)"]])
), "lm-smoke.drop1.tsv")

# --- tassel-6#186: the model layer's designs -------------------------------------------------------
#
# One data file, one model description, and the lm()/anova()/drop1() output per design. Factor levels
# are sorted in the C locale and the reference level of treatment coding is the LAST present level,
# which is TASSEL 5's convention (the last level's effect is zero), so coefficients line up with
# TASSEL 6's `Contrast.treatment`. `keep = 0` rows are the ones a caller leaves out, so a level can
# be empty after row selection; R drops unused levels (lm's drop.unused.levels), as TASSEL 6 does.
# log10p is computed from R's own log-scale tail, in base 10. Test sizes and seeds are fixed here.
set.seed(186)

# TASSEL 6 stores numeric phenotype values as Float (ADR-0023), so the data R fits are the values
# rounded to a 32-bit float, which print exactly with %.17g and read back as the same Double.
to_float <- function(x) readBin(writeBin(as.numeric(x), raw(), size = 4), "numeric", size = 4, n = length(x))
radix_factor <- function(x) factor(x, levels = sort(unique(x), method = "radix"))
log10_upper_f <- function(f, df1, df2) pf(f, df1, df2, lower.tail = FALSE, log.p = TRUE) / log(10)
log10_two_sided_t <- function(t, df) (log(2) + pt(abs(t), df, lower.tail = FALSE, log.p = TRUE)) / log(10)

emit_model <- function(name, data, formula, keep = rep(TRUE, nrow(data)), contrast = "treatment",
                       weights = NULL, with_drop1 = TRUE) {
  frame <- data
  for (n in names(frame)) if (is.numeric(frame[[n]])) frame[[n]] <- to_float(frame[[n]])
  frame$keep <- as.numeric(keep)
  if (!is.null(weights)) frame$w <- to_float(weights)
  columns <- lapply(names(frame), function(n) {
    v <- frame[[n]]
    if (is.numeric(v)) fmt(v) else as.character(v)
  })
  names(columns) <- names(frame)
  write_tsv(as.data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE),
            paste0(name, ".data.tsv"))
  writeLines(c(paste0("formula = ", formula),
               paste0("contrast = ", contrast),
               paste0("weights = ", if (is.null(weights)) "none" else "w")),
             file.path(out, paste0(name, ".model.txt")), sep = "\n", useBytes = TRUE)

  sub <- frame[frame$keep == 1, , drop = FALSE]
  for (n in names(sub)) {
    if (is.character(sub[[n]])) {
      f <- droplevels(radix_factor(sub[[n]]))
      contrasts(f) <- if (contrast == "sum") contr.sum(levels(f)) else
        contr.treatment(levels(f), base = nlevels(f))
      sub[[n]] <- f
    }
  }
  fit <- if (is.null(weights)) lm(as.formula(formula), data = sub) else
    lm(as.formula(formula), data = sub, weights = w)
  dfres <- fit$df.residual

  estimate <- coef(fit)
  summary_table <- summary(fit)$coefficients
  se <- t <- p <- rep(NA_real_, length(estimate))
  names(se) <- names(t) <- names(p) <- names(estimate)
  se[rownames(summary_table)] <- summary_table[, 2]
  t[rownames(summary_table)] <- summary_table[, 3]
  p[rownames(summary_table)] <- summary_table[, 4]
  write_tsv(data.frame(term = names(estimate), estimate = fmt(estimate), se = fmt(se), t = fmt(t),
                       p = fmt(p), log10p = fmt(log10_two_sided_t(t, dfres))),
            paste0(name, ".coefficients.tsv"))

  sequential <- anova(fit)
  write_tsv(data.frame(
    term = rownames(sequential), df = fmt(sequential[["Df"]]), ss = fmt(sequential[["Sum Sq"]]),
    ms = fmt(sequential[["Mean Sq"]]), f = fmt(sequential[["F value"]]),
    p = fmt(sequential[["Pr(>F)"]]),
    log10p = fmt(log10_upper_f(sequential[["F value"]], sequential[["Df"]], dfres))
  ), paste0(name, ".anova.tsv"))

  if (with_drop1) {
    marginal <- drop1(fit, test = "F")
    write_tsv(data.frame(
      term = rownames(marginal), df = fmt(marginal[["Df"]]), ss = fmt(marginal[["Sum of Sq"]]),
      rss = fmt(marginal[["RSS"]]), f = fmt(marginal[["F value"]]), p = fmt(marginal[["Pr(>F)"]]),
      log10p = fmt(log10_upper_f(marginal[["F value"]], marginal[["Df"]], dfres))
    ), paste0(name, ".drop1.tsv"))
  }
}

levels_of <- function(n, labels, prob = NULL) {
  repeat {
    x <- sample(labels, n, replace = TRUE, prob = prob)
    if (length(unique(x)) == length(labels)) return(x)
  }
}

# Covariates only.
n <- 30
x1 <- round(rnorm(n, 10, 2), 3)
x2 <- round(runif(n, 0, 5), 3)
emit_model("lm-covariates",
           data.frame(y = round(3 + 1.5 * x1 - 2 * x2 + rnorm(n), 4), x1 = x1, x2 = x2),
           "y ~ x1 + x2")

# Factors only, two factors, unequal cell counts.
n <- 36
g1 <- levels_of(n, c("A", "B", "C"))
g2 <- levels_of(n, c("P", "Q", "R", "S"))
emit_model("lm-factors",
           data.frame(y = round(20 + as.numeric(factor(g1)) * 2 + as.numeric(factor(g2)) * 0.5 +
                                  rnorm(n, 0, 1.5), 4), g1 = g1, g2 = g2),
           "y ~ g1 + g2")

# Mixed covariates and factors, interleaved; once with treatment and once with sum coding.
n <- 40
x1 <- round(rnorm(n, 50, 8), 2)
x2 <- round(rnorm(n, 0, 1), 4)
g1 <- levels_of(n, c("L1", "L2", "L3"))
g2 <- levels_of(n, c("E1", "E2"))
mixed <- data.frame(
  y = round(5 + 0.2 * x1 + 3 * x2 + as.numeric(factor(g1)) - 2 * (g2 == "E2") + rnorm(n, 0, 2), 4),
  x1 = x1, g1 = g1, x2 = x2, g2 = g2)
emit_model("lm-mixed", mixed, "y ~ x1 + g1 + x2 + g2")
emit_model("lm-mixed-sum", mixed, "y ~ x1 + g1 + x2 + g2", contrast = "sum")

# A strong signal, so the p-values are tiny and the log-scale tail matters.
n <- 60
x1 <- round(rnorm(n, 0, 1), 4)
g1 <- levels_of(n, c("a", "b", "c", "d", "e"))
emit_model("lm-strong",
           data.frame(y = round(40 * x1 + 25 * as.numeric(factor(g1)) + rnorm(n, 0, 0.5), 5),
                      x1 = x1, g1 = g1),
           "y ~ x1 + g1")

# A covariate that is an exact combination of two others: x3 is aliased, so it has no coefficient.
n <- 24
x1 <- round(rnorm(n, 5, 1) * 8) / 8   # multiples of 1/8, so x3 is exact in 32-bit floats too
x2 <- round(rnorm(n, 2, 1) * 8) / 8
g1 <- levels_of(n, c("A", "B", "C"))
emit_model("lm-collinear",
           data.frame(y = round(1 + x1 - x2 + rnorm(n, 0, 0.7), 4), x1 = x1, x2 = x2,
                      x3 = 2 * x1 + x2, g1 = g1),
           "y ~ x1 + x2 + x3 + g1", with_drop1 = FALSE)

# A factor level left empty by row selection (all of "C" and "D" is excluded; "D" stays a level of
# the column), plus a second factor that loses a level's worth of rows.
n <- 48
x1 <- round(rnorm(n, 0, 1), 4)
g1 <- levels_of(n, c("A", "B", "C", "D"))
g2 <- levels_of(n, c("u", "v", "w"))
emit_model("lm-empty-level",
           data.frame(y = round(2 + x1 + as.numeric(factor(g1)) * 0.8 + rnorm(n, 0, 1), 4),
                      x1 = x1, g1 = g1, g2 = g2),
           "y ~ x1 + g1 + g2", keep = !(g1 %in% c("C", "D")))

# Unbalanced replicates: repeated taxa-like blocks with unequal cell sizes and factors that are
# correlated, so sequential and marginal sums of squares differ.
n <- 54
g1 <- levels_of(n, c("A", "B", "C"), prob = c(0.6, 0.3, 0.1))
g2 <- ifelse(runif(n) < ifelse(g1 == "A", 0.2, 0.75), "N", "M")
if (length(unique(g2)) < 2) g2[1] <- "N"
x1 <- round(rnorm(n, 3, 1), 3)
emit_model("lm-unbalanced",
           data.frame(y = round(10 + 2 * (g1 == "B") + 4 * (g1 == "C") + 3 * (g2 == "N") + x1 +
                                  rnorm(n, 0, 1.2), 4), g1 = g1, g2 = g2, x1 = x1),
           "y ~ g1 + g2 + x1")

# Weights, which the GLS kernel reproduces with V = diag(1 / w).
n <- 32
x1 <- round(rnorm(n, 4, 1.5), 3)
g1 <- levels_of(n, c("A", "B", "C"))
w <- round(runif(n, 0.5, 3), 3)
emit_model("lm-weighted",
           data.frame(y = round(1 + 2 * x1 + as.numeric(factor(g1)) + rnorm(n, 0, 1 / sqrt(w)), 4),
                      x1 = x1, g1 = g1),
           "y ~ x1 + g1", weights = w)
