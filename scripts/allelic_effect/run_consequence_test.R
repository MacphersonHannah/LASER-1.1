#!/usr/bin/env Rscript
## LASER consequence module: sample x locus functional characterisation of abTU-positive loci only. 
## So only test abTU positive ones. For each abTU-significant locus, project its K-transcript (number of transcripts)
## composition onto each feature eg ORF, NMD and tests whether the two haplotypes differ
## on that feature more than overdispersion allows.
##
## Null: for each axis, draw B Dirichlet-Multinomial tables under the shared-
## proportions null at the SAME theta used for detection, each haplotype drawn at
## its OWN observed total (so abTE is held fixed and cannot generate a hit), then
## project each null draw onto the axis. Two-sided p = fraction of |null stat| >= |observed stat|.
## This uses the same alpha = props/theta idiom as
## compute_dm_lr in run_allelic_effect_test.R.  consequence and detection share
## one overdispersion model.

## features currently available are:
## categorical (within-allele class fraction diff): coding, nmd
## continuous  (within-allele distribution-shape test): tss, tes, orf_length
## An axis is tested only if it varies across the locus's transcripts; otherwise
## it is silently skipped (ie a locus whose transcripts are all coding on both alleles is not tested).
## coding: is_coding (ORFannotate's raw coding_class call, untouched) vs everything else: "is there a shift from coding to non-coding."
## nmd: is_nmd vs not, tested ONLY among the transcripts ORFannotate called coding (is_coding): "of the coding ones, what proportion switches to NMD-sensitive":  NOT a fraction over every transcript
## at the locus, which would dilute the answer whenever genuinely non-coding transcripts are also present.
##
## Continuous axes (tss/tes/orf_length) are tested using a two-sample Anderson-Darling or Kolmogorov-Smirnov statistic (--conseq_dist_test), user gets to choose which but AD recommended.
## Both statistics are computed with an EXACT tie correction operating directly on the per-transcript (feature_value, hap1_reads, hap2_reads)
## table , never expanded to a per-read vector. This matters here because almost all of the "sample size" at a locus is read-depth-driven
## duplication at a handful (K = number of transcripts) of distinct
## feature values, not independent continuous draws; the untied textbook AD/KS formulas assume no repeated values and are undefined (or depend on arbitrary tie-breaking) the moment two reads share a transcript, ie
## almost always here. The tie-corrected formulas used below (Scholz & Stephens 1987 for AD; the direct empirical-CDF definition for KS)
## collapse cleanly onto the untied formulas as a special case, verified numerically against textbook references before use.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(jsonlite)
  library(readr)
})

opt_list <- list(
  make_option("--sample_rds",        type = "character"),
  make_option("--allelic_results",   type = "character"),  # allelic_effect_results.tsv
  make_option("--transcript_features", type = "character"),# transcript_features.tsv
  make_option("--theta_json",        type = "character"),
  make_option("--sample",            type = "character"),
  make_option("--conseq_perm_b",     type = "integer", default = 2000),
  make_option("--conseq_perm_seed",  type = "integer", default = 123),
  make_option("--conseq_fdr_threshold", type = "double", default = 0.05),
  make_option("--conseq_dist_test",  type = "character", default = "ad",
              help = "Distribution test for continuous axes (tss/tes/orf_length): 'ad' (Anderson-Darling, default) or 'ks' (Kolmogorov-Smirnov)"),
  make_option("--out_results",  type = "character", default = "consequence_results.tsv"),
  make_option("--out_manifest", type = "character", default = "consequence_manifest.tsv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

if (!(opt$conseq_dist_test %in% c("ad", "ks")))
  stop(sprintf("--conseq_dist_test must be 'ad' or 'ks', got '%s'", opt$conseq_dist_test))

d        <- readRDS(opt$sample_rds)
tx_input <- d$tx_input
results  <- readr::read_tsv(opt$allelic_results, show_col_types = FALSE)
features <- readr::read_tsv(opt$transcript_features, show_col_types = FALSE) %>%
  mutate(bulk_pbid = as.character(bulk_pbid))

## gene_name lookup. Falls back to the Ensembl ID whenever no GENCODE symbol is available (mirrors the fallback applied in
## prepare_loci_tables.R, in case an older allelic_results file predates it).
gene_name_lookup <- results %>%
  distinct(associated_gene, gene_name) %>%
  mutate(gene_name = if_else(is.na(gene_name) | gene_name %in% c("", "NA"),
                              associated_gene, gene_name))

theta_info   <- jsonlite::read_json(opt$theta_json, simplifyVector = TRUE)
sample_theta <- theta_info$theta
if (!is.null(theta_info$per_sample_theta) &&
    opt$sample %in% names(theta_info$per_sample_theta)) {
  sample_theta <- theta_info$per_sample_theta[[opt$sample]]
}
theta_provenance <- theta_info$theta_mode

## abTU-positive loci for this sample
abtu_pos_genes <- results %>%
  filter(!is.na(abtu_sig), abtu_sig) %>%
  distinct(associated_gene) %>%
  pull(associated_gene)

## transcript counts + features, restricted to abTU-positive loci
conseq_input <- tx_input %>%
  filter(associated_gene %in% abtu_pos_genes) %>%
  left_join(features, by = "bulk_pbid")

## DM null draw: allocate N reads across K transcripts at (props, theta)
## mirrors alpha = props/theta from compute_dm_lr; Dirichlet via gamma normalisation
draw_dm_counts <- function(N, props, theta) {
  props[is.na(props)] <- 0
  props[props <= 0]   <- 1e-10
  props <- props / sum(props)
  alpha <- props / theta
  g <- stats::rgamma(length(props), shape = alpha, rate = 1)
  if (sum(g) <= 0) g <- rep(1, length(g))
  p <- g / sum(g)
  as.integer(stats::rmultinom(1, size = N, prob = p))
}

## projections (return NA when the axis doesn't vary or is unusable)
proj_class_diff <- function(n1, n2, cls) {           # cls: logical, TRUE = class-of-interest
  ok <- !is.na(cls)
  if (sum(ok) < 1 || length(unique(cls[ok])) < 2) return(NA_real_)
  n1 <- n1[ok]; n2 <- n2[ok]; cls <- cls[ok]
  s1 <- sum(n1); s2 <- sum(n2)
  if (s1 <= 0 || s2 <= 0) return(NA_real_)
  sum(n1[cls]) / s1 - sum(n2[cls]) / s2
}
## Same test, but restricted to a subset of transcripts
## first. Used for the nmd axis: The DM-bootstrap null in run_axis() still draws over the
## FULL locus (same theta, same pooled proportions across all K transcripts,
## matching every other axis); restricting happens only in the projection,
## applied identically to the observed data and every null draw - exactly
## the same pattern the continuous axes already use internally to drop NA
## feature values.
proj_class_diff_within <- function(n1, n2, cls, universe) {
  keep <- !is.na(universe) & universe
  proj_class_diff(n1[keep], n2[keep], cls[keep])
}
proj_ad_diff <- function(n1, n2, feat) {          # feat: numeric per-transcript
  ok <- !is.na(feat)
  if (sum(ok) < 1) return(NA_real_)
  n1 <- n1[ok]; n2 <- n2[ok]; feat <- feat[ok]
  if (length(unique(feat)) < 2) return(NA_real_)
  s1 <- sum(n1); s2 <- sum(n2)
  if (s1 <= 0 || s2 <= 0) return(NA_real_)
  ## Collapse to one row per distinct feature value, summing read counts
  ## within each group first. required whenever two different transcripts
  ## happen to share an exact TSS/TES/orf_length value, and required
  ## in effect since most "observations" here are repeated reads on the same
  ## transcript. rowsum() returns rows sorted ascending by group value.
  grp <- rowsum(cbind(n1, n2), group = feat)
  f1 <- grp[, 1]; f2 <- grp[, 2]
  fj <- f1 + f2
  N  <- s1 + s2
  Bj <- cumsum(fj)          # cumulative total reads through group j
  Mj <- cumsum(f1)          # cumulative hap1 reads through group j
  L  <- length(fj)
  if (L < 2) return(NA_real_)
  j <- seq_len(L - 1)       # the final group (j=L) always contributes 0, drop it
  ## Scholz & Stephens (1987) two-sample Anderson-Darling statistic with the
  ## tie correction (weight each distinct-value group by its size fj, using
  ## group-level cumulative counts Bj/Mj in place of per-observation ranks).
  ## Verified numerically against the untied textbook formula (agrees exactly
  ## whenever there are no ties) and against base R's ks.test for the KS
  ## counterpart below.
  sum(fj[j] * (N * Mj[j] - s1 * Bj[j])^2 / (Bj[j] * (N - Bj[j]))) / (s1 * s2)
}
proj_ks_diff <- function(n1, n2, feat) {          # feat: numeric per-transcript
  ok <- !is.na(feat)
  if (sum(ok) < 1) return(NA_real_)
  n1 <- n1[ok]; n2 <- n2[ok]; feat <- feat[ok]
  if (length(unique(feat)) < 2) return(NA_real_)
  s1 <- sum(n1); s2 <- sum(n2)
  if (s1 <= 0 || s2 <= 0) return(NA_real_)
  grp <- rowsum(cbind(n1, n2), group = feat)
  f1 <- grp[, 1]; f2 <- grp[, 2]
  F1 <- cumsum(f1) / s1     # read-weighted empirical CDF, haplotype 1
  F2 <- cumsum(f2) / s2     # read-weighted empirical CDF, haplotype 2
  max(abs(F1 - F2))         # classic two-sample KS D statistic
}
dist_fn <- switch(opt$conseq_dist_test, ad = proj_ad_diff, ks = proj_ks_diff)
## median reported as descriptive effect size only (NOT the tested statistic)
wmedian <- function(n, feat) {
  ok <- !is.na(feat) & n > 0
  if (!any(ok)) return(NA_real_)
  o <- order(feat[ok]); f <- feat[ok][o]; w <- n[ok][o]
  cw <- cumsum(w) / sum(w)
  f[which(cw >= 0.5)[1]]
}

## one axis test at one locus
run_axis <- function(n1, n2, projection, theta, B) {
  obs <- projection(n1, n2)
  if (is.na(obs)) return(NULL)
  N1 <- sum(n1); N2 <- sum(n2)
  pooled <- (n1 + n2) / (N1 + N2)
  null <- vapply(seq_len(B), function(b) {
    projection(draw_dm_counts(N1, pooled, theta),
               draw_dm_counts(N2, pooled, theta))
  }, numeric(1))
  null <- null[!is.na(null)]
  if (length(null) == 0) return(NULL)
  list(
    delta     = obs,
    direction = sign(obs),
    p_value   = (sum(abs(null) >= abs(obs)) + 1) / (length(null) + 1),
    null_sd   = stats::sd(null),
    n_tx      = length(n1)
  )
}

##  iterate loci x axes
set.seed(opt$conseq_perm_seed)

genes <- split(conseq_input, conseq_input$associated_gene)

rows <- list()
for (g in names(genes)) {
  dd <- genes[[g]]
  n1 <- dd$tx_hap1_reads; n2 <- dd$tx_hap2_reads

  axis_specs <- list(
    list(axis = "coding",    type = "categorical",
         proj = function(a, b) proj_class_diff(a, b, dd$is_coding), med = NA),
    list(axis = "nmd",       type = "categorical",
         proj = function(a, b) proj_class_diff_within(a, b, dd$is_nmd, dd$is_coding), med = NA),
    list(axis = "tss",       type = "continuous",
         proj = function(a, b) dist_fn(a, b, dd$tss),       med = dd$tss),
    list(axis = "tes",       type = "continuous",
         proj = function(a, b) dist_fn(a, b, dd$tes),       med = dd$tes),
    list(axis = "orf_length", type = "continuous",
         proj = function(a, b) dist_fn(a, b, dd$orf_length), med = dd$orf_length)
  )

  for (spec in axis_specs) {
    r <- run_axis(n1, n2, spec$proj, sample_theta, opt$conseq_perm_b)
    if (is.null(r)) next
    med_delta <- if (spec$type == "continuous")
      wmedian(n1, spec$med) - wmedian(n2, spec$med) else NA_real_
    ## AD/KS are unsigned magnitudes (always >= 0), so sign(r$delta) is
    ## meaningless for continuous axes - direction there comes from the
    ## signed median shift instead. Categorical axes are unaffected: their
    ## proj_class_diff() statistic is already a signed fraction difference.
    direction <- if (spec$type == "continuous") sign(med_delta) else r$direction
    rows[[length(rows) + 1]] <- tibble(
      sample                = opt$sample,
      associated_gene       = g,
      axis                  = spec$axis,
      axis_type             = spec$type,
      conseq_delta          = r$delta,          # tested statistic: AD/KS stat (continuous) or frac diff (categorical)
      conseq_delta_median   = med_delta,        # descriptive only
      conseq_direction      = direction,
      conseq_p_value        = r$p_value,
      conseq_null_sd        = r$null_sd,
      n_transcripts_tested  = r$n_tx,
      conseq_dist_test      = if (spec$type == "continuous") opt$conseq_dist_test else NA_character_,
      conseq_theta          = sample_theta,
      conseq_theta_provenance = theta_provenance,
      conseq_perm_b         = opt$conseq_perm_b
    )
  }
}

conseq_results <- if (length(rows)) bind_rows(rows) else tibble(
  sample = character(), associated_gene = character(), axis = character(),
  axis_type = character(), conseq_delta = double(), conseq_delta_median = double(),
  conseq_direction = double(), conseq_p_value = double(), conseq_null_sd = double(),
  n_transcripts_tested = integer(), conseq_dist_test = character(),
  conseq_theta = double(), conseq_theta_provenance = character(), conseq_perm_b = integer(),
  conseq_fdr = double(), conseq_sig = logical()
)

## per-sample FDR, computed WITHIN axis
if (nrow(conseq_results)) {
  conseq_results <- conseq_results %>%
    group_by(axis) %>%
    mutate(conseq_fdr = stats::p.adjust(conseq_p_value, method = "BH"),
           conseq_sig = !is.na(conseq_fdr) & conseq_fdr < opt$conseq_fdr_threshold) %>%
    ungroup() %>%
    arrange(associated_gene, axis)
}

## attach gene_name (fallback to associated_gene if no GENCODE symbol)
conseq_results <- conseq_results %>%
  left_join(gene_name_lookup, by = "associated_gene") %>%
  mutate(gene_name = if_else(is.na(gene_name), associated_gene, gene_name)) %>%
  relocate(gene_name, .after = associated_gene)

readr::write_tsv(conseq_results, opt$out_results)

## manifest
axis_tally <- conseq_results %>%
  group_by(axis) %>%
  summarise(n_tested = n(), n_sig = sum(conseq_sig, na.rm = TRUE), .groups = "drop")

manifest_row <- tibble(
  sample              = opt$sample,
  theta_used          = sample_theta,
  theta_provenance    = theta_provenance,
  conseq_dist_test    = opt$conseq_dist_test,
  n_abtu_pos_loci     = length(abtu_pos_genes),
  n_axis_tests        = nrow(conseq_results),
  n_sig_total         = sum(conseq_results$conseq_sig, na.rm = TRUE),
  n_tested_coding     = axis_tally$n_tested[match("coding",    axis_tally$axis)] %||% 0L,
  n_tested_nmd        = axis_tally$n_tested[match("nmd",       axis_tally$axis)] %||% 0L,
  n_tested_tss        = axis_tally$n_tested[match("tss",       axis_tally$axis)] %||% 0L,
  n_tested_tes        = axis_tally$n_tested[match("tes",       axis_tally$axis)] %||% 0L,
  n_tested_orf_length = axis_tally$n_tested[match("orf_length", axis_tally$axis)] %||% 0L,
  n_sig_coding        = axis_tally$n_sig[match("coding",    axis_tally$axis)] %||% 0L,
  n_sig_nmd           = axis_tally$n_sig[match("nmd",       axis_tally$axis)] %||% 0L,
  n_sig_tss           = axis_tally$n_sig[match("tss",       axis_tally$axis)] %||% 0L,
  n_sig_tes           = axis_tally$n_sig[match("tes",       axis_tally$axis)] %||% 0L,
  n_sig_orf_length    = axis_tally$n_sig[match("orf_length", axis_tally$axis)] %||% 0L
) %>%
  mutate(across(starts_with("n_tested_") | starts_with("n_sig_"), ~replace_na(., 0L))) %>%
  mutate(
    # Same denominator (n_abtu_pos_loci) for every axis, on purpose, so the
    # five prop_sig_* numbers are directly comparable to each other at a
    # glance: "X% of this sample's abTU+ loci show a significant coding /
    # nmd / tss / tes / orf_length change". NA (not 0) when n_abtu_pos_loci
    # == 0, so an empty sample reads as "no data" rather than "0% everything".
    #
    # axes differ in what fraction of the abTU+ loci can even be tested on them,
    # and nmd in particular is tested on a much smaller slice than the
    # others (it only runs where a locus has >=2 transcripts ORFannotate
    # called coding AND those coding transcripts vary in NMD status. A low prop_sig_nmd next to a low
    # prop_tested_nmd means "we mostly couldn't test this", not "NMD
    # enrichment is rare".
    prop_tested_coding     = if_else(n_abtu_pos_loci > 0, n_tested_coding     / n_abtu_pos_loci, NA_real_),
    prop_tested_nmd        = if_else(n_abtu_pos_loci > 0, n_tested_nmd        / n_abtu_pos_loci, NA_real_),
    prop_tested_tss        = if_else(n_abtu_pos_loci > 0, n_tested_tss        / n_abtu_pos_loci, NA_real_),
    prop_tested_tes        = if_else(n_abtu_pos_loci > 0, n_tested_tes        / n_abtu_pos_loci, NA_real_),
    prop_tested_orf_length = if_else(n_abtu_pos_loci > 0, n_tested_orf_length / n_abtu_pos_loci, NA_real_),
    prop_sig_coding        = if_else(n_abtu_pos_loci > 0, n_sig_coding        / n_abtu_pos_loci, NA_real_),
    prop_sig_nmd           = if_else(n_abtu_pos_loci > 0, n_sig_nmd           / n_abtu_pos_loci, NA_real_),
    prop_sig_tss           = if_else(n_abtu_pos_loci > 0, n_sig_tss           / n_abtu_pos_loci, NA_real_),
    prop_sig_tes           = if_else(n_abtu_pos_loci > 0, n_sig_tes           / n_abtu_pos_loci, NA_real_),
    prop_sig_orf_length    = if_else(n_abtu_pos_loci > 0, n_sig_orf_length    / n_abtu_pos_loci, NA_real_)
  )
readr::write_tsv(manifest_row, opt$out_manifest)

cat(sprintf(
  "[%s] consequence (dist_test=%s): %d abTU+ loci | %d axis tests | %d sig (coding %d, nmd %d, tss %d, tes %d, orf_len %d) | theta=%.6f (%s)\n",
  opt$sample, opt$conseq_dist_test, length(abtu_pos_genes), nrow(conseq_results),
  manifest_row$n_sig_total, manifest_row$n_sig_coding, manifest_row$n_sig_nmd,
  manifest_row$n_sig_tss, manifest_row$n_sig_tes, manifest_row$n_sig_orf_length,
  sample_theta, theta_provenance
))