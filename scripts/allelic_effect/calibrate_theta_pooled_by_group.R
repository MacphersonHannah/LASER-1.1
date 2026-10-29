#!/usr/bin/env Rscript
## Calibrate theta separately for each named sample group (eg brain, heart, other_peripheral). run once per theta_group instead of collapsed into
## one pooled group. per_sample_theta is populated from the result, so RUN_ALLELIC_EFFECT_TEST / RUN_CONSEQUENCE_TEST need no changes, they
## already check per_sample_theta before falling back to a single global theta.
## should not be run on a single sample in a group, but that's explained in the README

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(purrr)
  library(readr)
  library(jsonlite)
})

opt_list <- list(
  make_option("--combined_rds", type = "character"),
  make_option("--group_table",  type = "character"),  # CSV: sample,group
  make_option("--max_loci",      type = "integer", default = 2000),
  make_option("--n_null_ref",    type = "integer", default = 500),
  make_option("--n_test",        type = "integer", default = 50),
  make_option("--nominal_fpr",   type = "double",  default = 0.05),
  make_option("--seed",          type = "integer", default = 456),
  make_option("--out",           type = "character", default = "theta_pooled_by_group.json")
)
opt <- parse_args(OptionParser(option_list = opt_list))

d <- readRDS(opt$combined_rds)
locus_by_sample <- d$locus_by_sample
tx_input        <- d$tx_input

group_table <- readr::read_csv(opt$group_table, show_col_types = FALSE) %>%
  mutate(across(everything(), as.character))
if (!all(c("sample", "group") %in% names(group_table))) {
  stop("--group_table must be a CSV with columns: sample,group")
}

run_samples <- unique(locus_by_sample$sample_short)
missing <- setdiff(run_samples, group_table$sample)
if (length(missing) > 0) {
  stop(sprintf(
    "--theta_group_table is missing group assignment(s) for: %s. Every sample in this run needs a row.",
    paste(missing, collapse = ", ")
  ))
}

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

theta_candidates <- exp(seq(log(1e-3), log(0.2), length.out = 20))

## calibration function
## unchanged apart from taking hap1/hap2 lists directly instead of pulling
## them off a pre-built abtu_base_with_tissue object
calibrate_theta_fpr <- function(hap1_list, hap2_list,
                                 theta_grid  = theta_candidates,
                                 n_test      = opt$n_test,
                                 n_null_ref  = opt$n_null_ref,
                                 max_loci    = opt$max_loci,
                                 seed        = opt$seed) {
  set.seed(seed)
  valid <- vapply(seq_along(hap1_list), function(i) {
    h1   <- as.numeric(hap1_list[[i]])
    h2   <- as.numeric(hap2_list[[i]])
    keep <- (h1 + h2) > 0
    sum(keep) >= 2 && sum(h1[keep]) > 0 && sum(h2[keep]) > 0
  }, logical(1))
  idx <- which(valid)
  if (length(idx) > max_loci) idx <- sample(idx, max_loci)
  n_loci  <- length(idx)
  n_theta <- length(theta_grid)
  fpr_matrix <- matrix(NA_real_, nrow = n_loci, ncol = n_theta)

  for (li in seq_along(idx)) {
    i  <- idx[li]
    h1 <- round(as.numeric(hap1_list[[i]]))
    h2 <- round(as.numeric(hap2_list[[i]]))
    keep <- (h1 + h2) > 0
    h1 <- h1[keep]; h2 <- h2[keep]
    n1 <- sum(h1); n2 <- sum(h2); total <- h1 + h2

    all_tables  <- stats::r2dtable(n_null_ref + n_test, total, c(n1, n2))
    ref_tables  <- all_tables[seq_len(n_null_ref)]
    test_tables <- all_tables[seq(n_null_ref + 1L, n_null_ref + n_test)]

    for (ti in seq_along(theta_grid)) {
      th <- theta_grid[ti]
      ref_lrs  <- vapply(ref_tables,  function(pt) compute_dm_lr(pt[, 1L], pt[, 2L], th), numeric(1))
      test_lrs <- vapply(test_tables, function(pt) compute_dm_lr(pt[, 1L], pt[, 2L], th), numeric(1))
      test_pvals <- vapply(test_lrs, function(lr) {
        (sum(ref_lrs >= lr, na.rm = TRUE) + 1L) / (n_null_ref + 1L)
      }, numeric(1))
      fpr_matrix[li, ti] <- mean(test_pvals < 0.05, na.rm = TRUE)
    }
    if (li %% 100 == 0) cat(sprintf("  %d / %d loci done\n", li, n_loci))
  }

  tibble::tibble(
    theta    = theta_grid,
    mean_fpr = colMeans(fpr_matrix, na.rm = TRUE),
    sd_fpr   = apply(fpr_matrix, 2, sd, na.rm = TRUE),
    n_loci   = n_loci
  )
}

## Build hap1/hap2 count lists per group, then run per group
groups_to_calibrate <- sort(unique(group_table$group[group_table$sample %in% run_samples]))

cat(sprintf("\nRunning per-group permutation FPR calibration...\n"))
cat(sprintf("  %d loci max per group | %d null observations | %d theta candidates\n",
            opt$max_loci, opt$n_test, length(theta_candidates)))

calib_results <- purrr::map_dfr(groups_to_calibrate, function(grp) {
  cat(sprintf("[%s]\n", grp))
  grp_samples <- group_table$sample[group_table$group == grp & group_table$sample %in% run_samples]

  testable_keys <- locus_by_sample %>%
    filter(tested_status == "tested", sample_short %in% grp_samples) %>%
    select(sample, sample_short, associated_gene)

  tx_for_calib <- tx_input %>%
    semi_join(testable_keys, by = c("sample", "sample_short", "associated_gene")) %>%
    group_by(sample, sample_short, associated_gene) %>%
    summarise(
      hap1_counts = list(tx_hap1_reads),
      hap2_counts = list(tx_hap2_reads),
      .groups = "drop"
    )

  cat(sprintf("  %d loci available (%d samples)\n", nrow(tx_for_calib), length(grp_samples)))
  if (nrow(tx_for_calib) == 0) {
    stop(sprintf("Group '%s' has zero testable loci - check its samples passed the tested/read-count filters.", grp))
  }

  calibrate_theta_fpr(tx_for_calib$hap1_counts, tx_for_calib$hap2_counts) %>%
    mutate(theta_group = grp)
})

## Select calibrated theta per group (closest achieved FPR to nominal - in params but likely 0.05)
calibrated_thetas <- calib_results %>%
  group_by(theta_group) %>%
  slice_min(abs(mean_fpr - opt$nominal_fpr), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(theta_group, calibrated_theta = theta, calibrated_fpr = mean_fpr, n_loci) %>%
  mutate(
    fpr_deviation      = abs(calibrated_fpr - opt$nominal_fpr),
    calibration_status = case_when(
      fpr_deviation <= 0.01 ~ "well-calibrated (+/-1%)",
      fpr_deviation <= 0.02 ~ "acceptable (+/-2%)",
      TRUE                  ~ "review recommended"
    )
  )

cat("\n============================================================\n")
cat(" CALIBRATED THETA PER GROUP\n")
cat("============================================================\n")
print(calibrated_thetas)

## per_sample_theta for the pipeline's theta.json contract
per_sample_theta <- setNames(
  as.list(calibrated_thetas$calibrated_theta[match(group_table$group[group_table$sample %in% run_samples], calibrated_thetas$theta_group)]),
  group_table$sample[group_table$sample %in% run_samples]
)

result <- list(
  theta_mode        = "calibrate_pooled_by_group",
  theta             = NA,   # no single global value - always per-sample/per-group
  groups            = calibrated_thetas,
  per_sample_theta  = per_sample_theta,
  group_table       = opt$group_table,
  seed              = opt$seed,
  timestamp         = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)

jsonlite::write_json(result, opt$out, auto_unbox = TRUE, pretty = TRUE, na = "null")