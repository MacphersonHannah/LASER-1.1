#!/usr/bin/env Rscript
## Calibrate a single theta pooled across every sample/locus in this run.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(jsonlite)
})

opt_list <- list(
  make_option("--combined_rds", type = "character"),
  make_option("--max_loci",      type = "integer", default = 2000),
  make_option("--n_null_ref",    type = "integer", default = 500),
  make_option("--n_test",        type = "integer", default = 50),
  make_option("--nominal_fpr",   type = "double",  default = 0.05),
  make_option("--seed",          type = "integer", default = 456),
  make_option("--out",           type = "character", default = "theta_pooled.json")
)
opt <- parse_args(OptionParser(option_list = opt_list))

d <- readRDS(opt$combined_rds)
locus_by_sample <- d$locus_by_sample
tx_input        <- d$tx_input

compute_dm_lr <- function(h1, h2, theta) {
  dm_loglik <- function(counts, props, theta) {
    props[is.na(props)] <- 0
    props[props <= 0]   <- 1e-10
    props <- props / sum(props)
    n <- sum(counts, na.rm = TRUE)
    if (n <= 0) return(0)
    alpha <- props / theta
    lgamma(sum(alpha)) - lgamma(n + sum(alpha)) +
      sum(lgamma(counts + alpha) - lgamma(alpha))
  }
  n1 <- sum(h1, na.rm = TRUE); n2 <- sum(h2, na.rm = TRUE)
  if (n1 <= 0 || n2 <= 0) return(0)
  n_tot  <- n1 + n2
  pooled <- (h1 + h2) / n_tot
  p1 <- h1 / n1; p2 <- h2 / n2
  ll_h0 <- dm_loglik(h1, pooled, theta) + dm_loglik(h2, pooled, theta)
  ll_h1 <- dm_loglik(h1, p1, theta)     + dm_loglik(h2, p2, theta)
  max(2 * (ll_h1 - ll_h0), 0)
}

## Candidate loci: same "joint testable" definition as the testing step
## (tested_status == "tested"), pooled across every sample.
testable_keys <- locus_by_sample %>%
  filter(tested_status == "tested") %>%
  select(sample, sample_short, associated_gene)

tx_for_calib <- tx_input %>%
  semi_join(testable_keys, by = c("sample", "sample_short", "associated_gene")) %>%
  group_by(sample, sample_short, associated_gene) %>%
  summarise(
    hap1_counts = list(tx_hap1_reads),
    hap2_counts = list(tx_hap2_reads),
    .groups = "drop"
  )

theta_candidates <- exp(seq(log(1e-3), log(0.2), length.out = 20))

set.seed(opt$seed)

valid <- vapply(seq_len(nrow(tx_for_calib)), function(i) {
  h1 <- tx_for_calib$hap1_counts[[i]]; h2 <- tx_for_calib$hap2_counts[[i]]
  keep <- (h1 + h2) > 0
  sum(keep) >= 2 && sum(h1[keep]) > 0 && sum(h2[keep]) > 0
}, logical(1))

idx <- which(valid)
if (length(idx) > opt$max_loci) idx <- sample(idx, opt$max_loci)
n_loci  <- length(idx)
n_theta <- length(theta_candidates)

if (n_loci < 50) {
  warning(sprintf(
    "Only %d loci available for pooled theta calibration (recommend >=50 for a stable estimate).",
    n_loci
  ))
}

fpr_matrix <- matrix(NA_real_, nrow = n_loci, ncol = n_theta)

for (li in seq_along(idx)) {
  i  <- idx[li]
  h1 <- round(as.numeric(tx_for_calib$hap1_counts[[i]]))
  h2 <- round(as.numeric(tx_for_calib$hap2_counts[[i]]))
  keep <- (h1 + h2) > 0
  h1 <- h1[keep]; h2 <- h2[keep]
  n1 <- sum(h1); n2 <- sum(h2); total <- h1 + h2

  all_tables  <- stats::r2dtable(opt$n_null_ref + opt$n_test, total, c(n1, n2))
  ref_tables  <- all_tables[seq_len(opt$n_null_ref)]
  test_tables <- all_tables[seq(opt$n_null_ref + 1L, opt$n_null_ref + opt$n_test)]

  for (ti in seq_along(theta_candidates)) {
    th <- theta_candidates[ti]
    ref_lrs  <- vapply(ref_tables,  function(pt) compute_dm_lr(pt[, 1L], pt[, 2L], th), numeric(1))
    test_lrs <- vapply(test_tables, function(pt) compute_dm_lr(pt[, 1L], pt[, 2L], th), numeric(1))
    test_pvals <- vapply(test_lrs, function(lr) {
      (sum(ref_lrs >= lr, na.rm = TRUE) + 1L) / (opt$n_null_ref + 1L)
    }, numeric(1))
    fpr_matrix[li, ti] <- mean(test_pvals < 0.05, na.rm = TRUE)
  }
  if (li %% 100 == 0) cat(sprintf("  %d / %d loci done\n", li, n_loci))
}

mean_fpr <- colMeans(fpr_matrix, na.rm = TRUE)
best     <- which.min(abs(mean_fpr - opt$nominal_fpr))

result <- list(
  theta_mode        = "calibrate_pooled",
  theta             = theta_candidates[best],
  achieved_fpr      = mean_fpr[best],
  nominal_fpr       = opt$nominal_fpr,
  n_loci_used       = n_loci,
  n_samples_pooled  = n_distinct(tx_for_calib$sample_short[idx]),
  seed              = opt$seed,
  timestamp         = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  per_sample_theta  = NULL   ## always NULL in this mode - one theta for everyone
)

jsonlite::write_json(result, opt$out, auto_unbox = TRUE, pretty = TRUE)

cat(sprintf(
  "\nPooled theta = %.6f (achieved FPR = %.4f, target = %.2f) from %d loci across %d samples.\n",
  result$theta, result$achieved_fpr, opt$nominal_fpr, n_loci, result$n_samples_pooled
))