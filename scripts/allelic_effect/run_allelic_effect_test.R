#!/usr/bin/env Rscript
## Run abTE (binomial) + abTU (permutation Dirichlet Multinomial) testing for one sample, using the theta assigned to it (from whichever
## params.theta_mode produced theta.json), and classify each locus
## into the four-way allelic-effect scheme. Loci that failed the tested/not-tested gate in prepare_loci_tables.R are kept in the output
## with tested_status set and all test columns NA,  not dropped.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(purrr)
  library(jsonlite)
  library(readr)
})

opt_list <- list(
  make_option("--sample_rds",   type = "character"),
  make_option("--theta_json",   type = "character"),
  make_option("--sample",       type = "character"),
  make_option("--abte_fdr_threshold",       type = "double",  default = 0.05),
  make_option("--abte_effect_guard",        type = "logical", default = TRUE),
  make_option("--abte_major_frac_threshold", type = "double", default = 0.55),
  make_option("--abte_fdr_scope",           type = "character", default = "sample"),
  make_option("--abtu_fdr_threshold",       type = "double",  default = 0.05),
  make_option("--abtu_effect_guard",        type = "logical", default = FALSE),
  make_option("--abtu_tvd_threshold",       type = "double",  default = 0.20),
  make_option("--abtu_perm_b",              type = "integer", default = 1000),
  make_option("--dm_perm_seed",             type = "integer", default = 123),
  make_option("--out_results",  type = "character", default = "allelic_effect_results.tsv"),
  make_option("--out_manifest", type = "character", default = "sample_manifest.tsv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

d <- readRDS(opt$sample_rds)
locus_by_sample <- d$locus_by_sample
tx_input        <- d$tx_input

theta_info <- jsonlite::read_json(opt$theta_json, simplifyVector = TRUE)

## per_sample_theta (if the eventual fixed/lift modes populate it) takes priority over the single pooled/global value
sample_theta <- theta_info$theta
if (!is.null(theta_info$per_sample_theta) &&
    opt$sample %in% names(theta_info$per_sample_theta)) {
  sample_theta <- theta_info$per_sample_theta[[opt$sample]]
}
theta_provenance <- theta_info$theta_mode

## abTE: binomial test on tested loci only

tested_loci <- filter(locus_by_sample, tested_status == "tested")
untested_loci <- filter(locus_by_sample, tested_status != "tested")

run_abte_binom_test <- function(hap1_reads, hap2_reads) {
  n <- hap1_reads + hap2_reads
  if (is.na(hap1_reads) || is.na(hap2_reads) || is.na(n) || n <= 0) return(NA_real_)
  stats::binom.test(round(hap1_reads), round(n), p = 0.5, alternative = "two.sided")$p.value
}

abte_results <- tested_loci %>%
  mutate(
    hap1_fraction               = locus_hap1_reads / locus_phased_reads,
    hap2_fraction               = locus_hap2_reads / locus_phased_reads,
    major_allele_fraction       = pmax(hap1_fraction, hap2_fraction, na.rm = TRUE),
    allelic_fraction_difference = abs(hap1_fraction - hap2_fraction),
    abte_direction = case_when(
      hap1_fraction > hap2_fraction ~ "haplotype_1_higher",
      hap2_fraction > hap1_fraction ~ "haplotype_2_higher",
      TRUE                          ~ "balanced"
    ),
    abte_p_value = map2_dbl(locus_hap1_reads, locus_hap2_reads, run_abte_binom_test)
  ) %>%
  mutate(
    abte_fdr = if (opt$abte_fdr_scope == "sample") {
      stats::p.adjust(abte_p_value, method = "BH")
    } else {
      stats::p.adjust(abte_p_value, method = "BH")  ## single sample per process => same as 'global' here
    },
    abte_sig_stat_only         = !is.na(abte_fdr) & abte_fdr < opt$abte_fdr_threshold,
    abte_sig_with_effect_guard = abte_sig_stat_only & major_allele_fraction >= opt$abte_major_frac_threshold,
    abte_sig = if (isTRUE(opt$abte_effect_guard)) abte_sig_with_effect_guard else abte_sig_stat_only
  )

## abTU: permutation Dirichlet Multinomial test on tested loci only

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

calculate_tvd <- function(h1, h2) {
  n1 <- sum(h1, na.rm = TRUE); n2 <- sum(h2, na.rm = TRUE)
  if (n1 <= 0 || n2 <= 0) return(NA_real_)
  0.5 * sum(abs(h1 / n1 - h2 / n2), na.rm = TRUE)
}
calculate_js_distance <- function(h1, h2) {
  n1 <- sum(h1, na.rm = TRUE); n2 <- sum(h2, na.rm = TRUE)
  if (n1 <= 0 || n2 <= 0) return(NA_real_)
  p <- h1 / n1; q <- h2 / n2; m <- 0.5 * (p + q)
  kl_pm <- sum(ifelse(p > 0, p * log(p / m), 0), na.rm = TRUE)
  kl_qm <- sum(ifelse(q > 0, q * log(q / m), 0), na.rm = TRUE)
  sqrt(0.5 * kl_pm + 0.5 * kl_qm)
}
get_top_driver_delta <- function(h1, h2) {
  n1 <- sum(h1, na.rm = TRUE); n2 <- sum(h2, na.rm = TRUE)
  if (n1 <= 0 || n2 <= 0) return(NA_real_)
  max(abs(h1 / n1 - h2 / n2), na.rm = TRUE)
}
get_top_driver_transcript <- function(tx, h1, h2) {
  n1 <- sum(h1, na.rm = TRUE); n2 <- sum(h2, na.rm = TRUE)
  if (length(tx) == 0 || n1 <= 0 || n2 <= 0) return(NA_character_)
  tx[which.max(abs(h1 / n1 - h2 / n2))]
}

run_abtu_dm_test <- function(h1, h2, theta, B) {
  h1 <- round(as.numeric(h1)); h2 <- round(as.numeric(h2))
  keep <- (h1 + h2) > 0
  h1 <- h1[keep]; h2 <- h2[keep]
  if (length(h1) < 2 || sum(h1) <= 0 || sum(h2) <= 0) return(NA_real_)
  n1 <- sum(h1); n2 <- sum(h2); total <- h1 + h2
  observed_lr <- compute_dm_lr(h1, h2, theta)
  null_lr <- replicate(B, {
    pt <- stats::r2dtable(1, total, c(n1, n2))[[1]]
    compute_dm_lr(pt[, 1], pt[, 2], theta)
  })
  null_lr <- null_lr[!is.na(null_lr)]
  if (length(null_lr) == 0) return(NA_real_)
  (sum(null_lr >= observed_lr) + 1) / (length(null_lr) + 1)
}

abtu_input <- tx_input %>%
  semi_join(select(tested_loci, sample, sample_short, associated_gene),
            by = c("sample", "sample_short", "associated_gene")) %>%
  group_by(sample, sample_short, associated_gene) %>%
  summarise(
    transcript_ids = list(bulk_pbid),
    hap1_counts    = list(tx_hap1_reads),
    hap2_counts    = list(tx_hap2_reads),
    .groups = "drop"
  )

set.seed(opt$dm_perm_seed)

abtu_results <- abtu_input %>%
  mutate(
    abtu_p_value = pmap_dbl(
      list(hap1_counts, hap2_counts),
      function(h1, h2) run_abtu_dm_test(h1, h2, theta = sample_theta, B = opt$abtu_perm_b)
    ),
    abtu_tvd                   = map2_dbl(hap1_counts, hap2_counts, calculate_tvd),
    abtu_js_distance           = map2_dbl(hap1_counts, hap2_counts, calculate_js_distance),
    abtu_max_delta_prop        = map2_dbl(hap1_counts, hap2_counts, get_top_driver_delta),
    abtu_top_driver_transcript = pmap_chr(list(transcript_ids, hap1_counts, hap2_counts), get_top_driver_transcript),
    abtu_theta       = as.double(sample_theta),
    abtu_theta_provenance = theta_provenance,
    abtu_perm_b      = opt$abtu_perm_b
  ) %>%
  mutate(
    abtu_fdr = stats::p.adjust(abtu_p_value, method = "BH"),
    abtu_sig_stat_only         = !is.na(abtu_fdr) & abtu_fdr < opt$abtu_fdr_threshold,
    abtu_sig_with_effect_guard = abtu_sig_stat_only & abtu_tvd >= opt$abtu_tvd_threshold,
    abtu_sig = if (isTRUE(opt$abtu_effect_guard)) abtu_sig_with_effect_guard else abtu_sig_stat_only
  ) %>%
  select(-transcript_ids, -hap1_counts, -hap2_counts)

## Four way classification, joined back onto tested loci

tested_classified <- tested_loci %>%
  left_join(select(abte_results, sample, sample_short, associated_gene,
                    abte_p_value, abte_fdr, abte_sig,
                    major_allele_fraction, allelic_fraction_difference, abte_direction),
            by = c("sample", "sample_short", "associated_gene")) %>%
  left_join(select(abtu_results, sample, sample_short, associated_gene,
                    abtu_p_value, abtu_fdr, abtu_sig, abtu_tvd, abtu_js_distance,
                    abtu_max_delta_prop, abtu_top_driver_transcript,
                    abtu_theta, abtu_theta_provenance),
            by = c("sample", "sample_short", "associated_gene")) %>%
  mutate(
    allelic_effect_class = case_when(
      !abte_sig & !abtu_sig ~ "No detectable allelic effect",
      abte_sig  & !abtu_sig ~ "abTE only",
      !abte_sig & abtu_sig  ~ "abTU without detectable abTE",
      abte_sig  & abtu_sig  ~ "Both abTE and abTU",
      TRUE                  ~ NA_character_
    )
  )

## Untested loci: keep the row, NA the test columns

na_test_cols <- tibble(
  abte_p_value = NA_real_, abte_fdr = NA_real_, abte_sig = NA,
  major_allele_fraction = NA_real_, allelic_fraction_difference = NA_real_,
  abte_direction = NA_character_,
  abtu_p_value = NA_real_, abtu_fdr = NA_real_, abtu_sig = NA,
  abtu_tvd = NA_real_, abtu_js_distance = NA_real_, abtu_max_delta_prop = NA_real_,
  abtu_top_driver_transcript = NA_character_,
  abtu_theta = NA_real_, abtu_theta_provenance = NA_character_,
  allelic_effect_class = NA_character_
)

untested_classified <- bind_cols(untested_loci, na_test_cols[rep(1, nrow(untested_loci)), ])

final_results <- bind_rows(tested_classified, untested_classified) %>%
  arrange(chrom, gene_start, associated_gene)

readr::write_tsv(final_results, opt$out_results)

## Sample-level manifest row

manifest_row <- tibble(
  sample                = opt$sample,
  theta_used            = sample_theta,
  theta_provenance      = theta_provenance,
  n_loci_total          = nrow(final_results),
  n_tested              = sum(final_results$tested_status == "tested"),
  n_not_tested_gene_not_in_reference = sum(final_results$tested_status == "not_tested_gene_not_in_reference"),
  n_not_tested_too_few_tx    = sum(final_results$tested_status == "not_tested_too_few_transcripts"),
  n_not_tested_insufficient_reads = sum(final_results$tested_status == "not_tested_insufficient_reads"),
  # "sig" everywhere below means FDR-corrected significance (abte_fdr/abtu_fdr < the *_fdr_threshold param), NOT a nominal/raw p-value cutoff - and,
  # if abte_effect_guard/abtu_effect_guard was turned on for this run in config, it also requires passing the effect-size guard (major_allele_fraction/
  # abtu_tvd threshold). The exact threshold/guard settings actually used are recorded in the run-level params block appended to allelic_effect_summary.tsv by
  # MERGE_RUN_MANIFEST, so the definition is always recoverable from there. THis part should also be fairly fast to rerun
  n_no_effect           = sum(final_results$allelic_effect_class == "No detectable allelic effect", na.rm = TRUE),
  n_abte_sig_only       = sum(final_results$allelic_effect_class == "abTE only", na.rm = TRUE),
  n_abtu_sig_only       = sum(final_results$allelic_effect_class == "abTU without detectable abTE", na.rm = TRUE),
  n_both_sig            = sum(final_results$allelic_effect_class == "Both abTE and abTU", na.rm = TRUE),
  # Inclusive totals (each already includes the "Both" cases above) - use n_abte_sig_only/n_abtu_sig_only instead if you want mutually-exclusive
  # counts that sum to n_tested alongside n_no_effect and n_both_sig.
  n_abte_sig            = sum(final_results$abte_sig, na.rm = TRUE),
  n_abtu_sig            = sum(final_results$abtu_sig, na.rm = TRUE),
  # Proportions of TESTED loci in the four mutually exclusive categories above (n_no_effect + n_abte_sig_only + n_abtu_sig_only + n_both_sig ==
  # n_tested, so these four proportions always sum to 1 when n_tested > 0).
  # NA (not 0) when n_tested == 0, so an empty sample reads as "no data" rather than silently as "0% everything".
  prop_no_effect        = if_else(n_tested > 0, n_no_effect     / n_tested, NA_real_),
  prop_abte_sig_only    = if_else(n_tested > 0, n_abte_sig_only / n_tested, NA_real_),
  prop_abtu_sig_only    = if_else(n_tested > 0, n_abtu_sig_only / n_tested, NA_real_),
  prop_both_sig         = if_else(n_tested > 0, n_both_sig      / n_tested, NA_real_)
)
readr::write_tsv(manifest_row, opt$out_manifest)

cat(sprintf(
  "[%s] %d loci (%d tested, %d not tested) | theta=%.6f (%s) | no effect=%d abTE-sig only=%d abTU-sig only=%d both=%d\n",
  opt$sample, manifest_row$n_loci_total, manifest_row$n_tested,
  manifest_row$n_loci_total - manifest_row$n_tested,
  sample_theta, theta_provenance,
  manifest_row$n_no_effect, manifest_row$n_abte_sig_only, manifest_row$n_abtu_sig_only, manifest_row$n_both_sig
))