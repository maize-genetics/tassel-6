# R reference fixtures for the GLM kernel (tassel-6#191, ADR-0029).
#
# Usage: Rscript glm-fixtures.R OUT_DIR ASSOCIATION_DIR
#
# ASSOCIATION_DIR is testdata/association/ (the small hand-built inputs of the glm-synthetic goldens). Three
# groups of fixtures, all fitted with base R lm() and anova(), never with TASSEL code:
#
#   glm-replicates.*    GLM with replicated taxa, which TASSEL 5 refuses (decision Q1): the model is
#                       y ~ covariate + marker + taxon, the marker is tested against the taxon term added after
#                       it (taxa nested within the marker), as ADR-0026 describes.
#   glm-synthetic.*     the departures from TASSEL 5 in the glm-synthetic goldens: allele-copy coding of add_F and
#                       dom_F at a multiallelic site, and the addEffect/domEffect formula applied to the marker's
#                       own coefficients (TASSEL 5 indexes them by the number of effects, not columns).
#   glm-collinear.*     the collinear-covariate golden with c3 recognized as aliased (R's lm uses tol = 1e-7; the
#                       values are single precision, so c3 = c1 + c2 holds only to about 1e-7).
#
# Every number is written with %.17g so the files are byte-identical across runs of the same image. Phenotype
# values are rounded to 32-bit floats first, as TASSEL 6 stores them (ADR-0023).
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
input <- args[[2]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
write_tsv <- function(frame, name) {
  lines <- c(paste(names(frame), collapse = "\t"),
             apply(frame, 1, function(r) paste(r, collapse = "\t")))
  writeLines(lines, file.path(out, name), sep = "\n", useBytes = TRUE)
}
to_float <- function(x) readBin(writeBin(as.numeric(x), raw(), size = 4), "numeric", size = 4, n = length(x))
log10_upper_f <- function(f, df1, df2) pf(f, df1, df2, lower.tail = FALSE, log.p = TRUE) / log(10)
upper_f <- function(f, df1, df2) pf(f, df1, df2, lower.tail = FALSE)

# --- genotype classes and allele copies -----------------------------------------------------------------------

iupac <- c(AC = "M", AG = "R", AT = "W", CG = "S", CT = "Y", GT = "K")
class_label <- function(call) {
  a <- substr(call, 1, 1)
  b <- substr(call, 2, 2)
  if (a == b) return(a)
  pair <- paste(sort(c(a, b)), collapse = "")
  if (grepl("[+-]", pair)) "0" else iupac[[pair]]
}
# TASSEL 5 sorts classes with A, C, G, T first, then the rest by character; the last class is the reference.
class_order <- function(labels) {
  unique_labels <- unique(labels)
  first <- unique_labels[unique_labels %in% c("A", "C", "G", "T")]
  first <- first[order(match(first, c("A", "C", "G", "T")))]
  rest <- sort(setdiff(unique_labels, first), method = "radix")
  c(first, rest)
}
copies <- function(calls, allele) {
  as.numeric(substr(calls, 1, 1) == allele) + as.numeric(substr(calls, 2, 2) == allele)
}
# Alleles by decreasing count over the calls, ties by TASSEL's nucleotide code (A, C, G, T, then the rest).
allele_order <- function(calls) {
  chars <- unlist(strsplit(calls[!grepl("N", calls)], ""))
  counts <- table(chars)
  codes <- match(names(counts), c("A", "C", "G", "T", "+", "-"))
  names(counts)[order(-as.numeric(counts), codes)]
}

# The sequential F tests of one site for the model `base` (a formula string of covariates and factors), data `d`
# with columns y and calls. `taxon` names the taxon column when taxa are nested in the marker (replicates).
site_tests <- function(d, base, nested = FALSE) {
  d <- d[!grepl("N", d$call), , drop = FALSE]
  if (all(d$call == d$call[1])) return(NULL)
  labels <- vapply(d$call, class_label, character(1))
  levels <- class_order(labels)
  d$marker <- factor(labels, levels = levels)
  contrasts(d$marker) <- contr.treatment(levels, base = length(levels))
  d$taxon <- factor(d$taxon)
  rhs <- function(extra) paste(c(base, extra, if (nested) "taxon"), collapse = " + ")
  full <- lm(as.formula(paste("y ~", rhs("marker"))), data = d)
  a <- anova(full)
  marker_row <- a["marker", ]
  if (nested) {
    error_ss <- a["taxon", "Sum Sq"]
    error_df <- a["taxon", "Df"]
  } else {
    error_ss <- a["Residuals", "Sum Sq"]
    error_df <- a["Residuals", "Df"]
  }
  model_ss <- sum(a[setdiff(rownames(a), "Residuals"), "Sum Sq"])
  total_ss <- sum(a[, "Sum Sq"])
  marker_df <- marker_row[["Df"]]
  marker_f <- (marker_row[["Sum Sq"]] / marker_df) / (error_ss / error_df)

  alleles <- allele_order(d$call)
  add_f <- add_df <- dom_f <- dom_df <- add_error_df <- NA_real_
  if (length(alleles) >= 2) {
    columns <- character(0)
    for (i in seq_len(length(alleles) - 1)) {
      name <- paste0("a", i)
      d[[name]] <- copies(d$call, alleles[i])
      columns <- c(columns, name)
    }
    additive <- lm(as.formula(paste("y ~", rhs(columns))), data = d)
    b <- anova(additive)
    add_ss <- sum(b[columns, "Sum Sq"])
    add_df <- sum(b[columns, "Df"])
    add_error_ss <- if (nested) b["taxon", "Sum Sq"] else b["Residuals", "Sum Sq"]
    add_error_df <- if (nested) b["taxon", "Df"] else b["Residuals", "Df"]
    add_f <- (add_ss / add_df) / (add_error_ss / add_error_df)
    dom_df <- marker_df - add_df
    # No dominance degrees of freedom (no heterozygote class): the test does not exist, as TASSEL 5 reports NaN.
    dom_f <- if (dom_df > 0) ((marker_row[["Sum Sq"]] - add_ss) / dom_df) / (error_ss / error_df) else NA_real_
  }
  list(n = nrow(d), classes = length(levels), marker_df = marker_df, error_df = error_df, marker_f = marker_f,
       marker_ss = marker_row[["Sum Sq"]], error_ss = error_ss, total_ss = total_ss, model_ss = model_ss,
       model_df = sum(a[setdiff(rownames(a), "Residuals"), "Df"]), add_df = add_df, add_f = add_f,
       add_error_df = add_error_df,
       dom_df = dom_df, dom_f = dom_f, fit = full, alleles = alleles)
}

stats_row <- function(trait, site, t) {
  data.frame(
    trait = trait, marker = site, n = t$n, classes = t$classes, marker_df = fmt(t$marker_df),
    error_df = fmt(t$error_df), marker_F = fmt(t$marker_f),
    marker_p = fmt(upper_f(t$marker_f, t$marker_df, t$error_df)),
    marker_log10p = fmt(log10_upper_f(t$marker_f, t$marker_df, t$error_df)),
    marker_Rsq = fmt(t$marker_ss / t$total_ss), marker_MS = fmt(t$marker_ss / t$marker_df),
    error_MS = fmt(t$error_ss / t$error_df), model_df = fmt(t$model_df), model_MS = fmt(t$model_ss / t$model_df),
    add_df = fmt(t$add_df), add_F = fmt(t$add_f),
    add_p = fmt(ifelse(is.na(t$add_f), NA, upper_f(t$add_f, t$add_df, t$add_error_df))),
    dom_df = fmt(t$dom_df), dom_F = fmt(t$dom_f),
    dom_p = fmt(ifelse(is.na(t$dom_f) || t$dom_df <= 0, NA, upper_f(t$dom_f, t$dom_df, t$error_df))),
    stringsAsFactors = FALSE)
}

# --- replicated taxa -------------------------------------------------------------------------------------------

set.seed(191)
n_taxa <- 24
taxa <- sprintf("T%02d", seq_len(n_taxa))
reps <- sample(c(2, 3), n_taxa, replace = TRUE)
r1 <- sample(c("AA", "AC", "CC"), n_taxa, replace = TRUE, prob = c(0.4, 0.2, 0.4))
r2 <- sample(c("AA", "CC"), n_taxa, replace = TRUE)
r3 <- sample(c("AA", "AC", "CC"), n_taxa, replace = TRUE)
r3[c(4, 11, 19)] <- "NN"
r4 <- sample(c("AA", "AG", "GG"), n_taxa, replace = TRUE, prob = c(0.6, 0.3, 0.1))
genotypes <- data.frame(taxon = taxa, r1 = r1, r2 = r2, r3 = r3, r4 = r4, stringsAsFactors = FALSE)
write_tsv(genotypes, "glm-replicates.genotypes.tsv")

taxon_effect <- rnorm(n_taxa, 0, 1.5)
rows <- data.frame(taxon = rep(taxa, reps), stringsAsFactors = FALSE)
index <- match(rows$taxon, taxa)
rows$c1 <- to_float(round(rnorm(nrow(rows), 5, 1.5), 3))
shift <- c(AA = 0, AC = 0.8, CC = 1.6)[r1[index]]
rows$y <- to_float(round(10 + 0.6 * rows$c1 + shift + taxon_effect[index] + rnorm(nrow(rows), 0, 1), 4))
write_tsv(data.frame(taxon = rows$taxon, y = fmt(rows$y), c1 = fmt(rows$c1)), "glm-replicates.phenotype.tsv")

out_rows <- list()
for (site in c("r1", "r2", "r3", "r4")) {
  d <- data.frame(y = rows$y, c1 = rows$c1, taxon = rows$taxon, call = genotypes[[site]][index],
                  stringsAsFactors = FALSE)
  out_rows[[site]] <- stats_row("y", site, site_tests(d, "c1", nested = TRUE))
}
write_tsv(do.call(rbind, out_rows), "glm-replicates.tests.tsv")

# --- the glm-synthetic departures ------------------------------------------------------------------------------

hmp <- read.delim(file.path(input, "glm_synthetic.hmp.txt"), header = TRUE, colClasses = "character",
                  check.names = FALSE)
syn_taxa <- colnames(hmp)[12:ncol(hmp)]
read_phenotype <- function(file) {
  lines <- readLines(file.path(input, file))
  header <- strsplit(lines[3], "\t")[[1]]
  body <- do.call(rbind, strsplit(lines[4:length(lines)], "\t"))
  frame <- as.data.frame(body, stringsAsFactors = FALSE)
  names(frame) <- header
  frame
}
traits <- read_phenotype("glm_synthetic_traits.txt")
numeric_column <- function(x) to_float(suppressWarnings(as.numeric(ifelse(x == "NA", NA, x))))

synthetic_rows <- list()
effect_rows <- list()
for (trait in c("y1", "y2")) {
  d <- data.frame(taxon = traits$Taxon, y = numeric_column(traits[[trait]]), c1 = numeric_column(traits$c1),
                  c2 = numeric_column(traits$c2), grp = traits$grp, stringsAsFactors = FALSE)
  keep <- !is.na(d$y) & !is.na(d$c1) & !is.na(d$c2)
  d <- d[keep, , drop = FALSE]
  d$grp <- factor(d$grp, levels = sort(unique(d$grp), method = "radix"))
  contrasts(d$grp) <- contr.treatment(levels(d$grp), base = nlevels(d$grp))
  for (i in seq_len(nrow(hmp))) {
    site <- hmp[i, 1]
    calls_all <- as.character(unlist(hmp[i, syn_taxa]))
    names(calls_all) <- syn_taxa
    d$call <- calls_all[d$taxon]
    d$taxon <- d$taxon
    t <- site_tests(d, "grp + c1 + c2")
    if (is.null(t)) next
    synthetic_rows[[paste(trait, site)]] <- stats_row(trait, site, t)
    # addEffect and domEffect from the marker's own coefficients, for a biallelic site (TASSEL 5's formulas).
    if (length(t$alleles) == 2 && t$classes %in% c(2, 3)) {
      b <- coef(t$fit)
      marker_coefficients <- b[grep("^marker", names(b))]
      b0 <- unname(marker_coefficients[1])
      b1 <- if (length(marker_coefficients) >= 2) unname(marker_coefficients[2]) else NA_real_
      add_effect <- if (t$classes == 2) abs(b0 / 2) else abs((b0 - b1) / 2)
      dom_effect <- if (t$classes == 3) -(b0 + b1) / 2 else NA_real_
      effect_rows[[paste(trait, site)]] <- data.frame(
        trait = trait, marker = site, classes = t$classes, b0 = fmt(b0), b1 = fmt(b1),
        addEffect = fmt(add_effect), domEffect = fmt(dom_effect), stringsAsFactors = FALSE)
    }
  }
}
write_tsv(do.call(rbind, synthetic_rows), "glm-synthetic.tests.tsv")
write_tsv(do.call(rbind, effect_rows), "glm-synthetic.effects.tsv")

# --- the collinear-covariate golden -----------------------------------------------------------------------------

collinear <- read_phenotype("glm_synthetic_collinear.txt")
d <- data.frame(taxon = collinear$Taxon, y = numeric_column(collinear$y1), c1 = numeric_column(collinear$c1),
                c2 = numeric_column(collinear$c2), c3 = numeric_column(collinear$c3), stringsAsFactors = FALSE)
collinear_rows <- list()
for (i in seq_len(nrow(hmp))) {
  site <- hmp[i, 1]
  calls_all <- as.character(unlist(hmp[i, syn_taxa]))
  names(calls_all) <- syn_taxa
  d$call <- calls_all[d$taxon]
  t <- site_tests(d, "c1 + c2 + c3")
  if (is.null(t)) next
  collinear_rows[[site]] <- stats_row("y1", site, t)
}
write_tsv(do.call(rbind, collinear_rows), "glm-collinear.tests.tsv")
