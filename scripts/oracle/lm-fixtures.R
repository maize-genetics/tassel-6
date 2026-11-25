# R lm() / anova() / drop1() reference fixtures (tassel-6#183, ADR-0029).
#
# Usage: Rscript lm-fixtures.R OUT_DIR
#
# Seed fixture only: one small design (a covariate and a three-level factor) that proves the
# pinned R container, the file layout, and the test-side loader. tassel-6#186 adds the designs,
# contrasts, collinear and missing-row cases it needs by extending this script. Every number is
# written with %.17g so the files are byte-identical across runs of the same image.
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
