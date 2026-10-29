#!/usr/bin/env Rscript
## theta_mode = 'lift': assign each sample the theta from a previously calibrated reference source,
## defaulting to --default_source unless overridden per sample.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(jsonlite)
})

opt_list <- list(
  make_option("--reference_csv",   type = "character"),
  make_option("--default_source",  type = "character", default = "all_encode"),
  make_option("--source_table",    type = "character", default = "NA"),  # path or "NA": CSV sample,source
  make_option("--samples",         type = "character"),                  # comma-separated sample names in this run
  make_option("--out",             type = "character", default = "theta_lift.json")
)
opt <- parse_args(OptionParser(option_list = opt_list))

run_samples <- strsplit(opt$samples, ",")[[1]]

reference <- readr::read_csv(opt$reference_csv, show_col_types = FALSE)

dup_sources <- reference$source[duplicated(reference$source)]
if (length(dup_sources) > 0) {
  stop(sprintf(
    "Duplicate source name(s) in %s: %s. Each 'source' must be unique - rename one before continuing.",
    opt$reference_csv, paste(unique(dup_sources), collapse = ", ")
  ))
}

if (!opt$default_source %in% reference$source) {
  stop(sprintf(
    "--theta_lift_source '%s' not found in reference asset (%s). Available: %s",
    opt$default_source, opt$reference_csv, paste(reference$source, collapse = ", ")
  ))
}

## table override, otherwise the default for everyone
sample_source <- setNames(rep(opt$default_source, length(run_samples)), run_samples)

has_table <- !is.na(opt$source_table) && opt$source_table != "NA" && nzchar(opt$source_table)
if (has_table) {
  tbl <- readr::read_csv(opt$source_table, show_col_types = FALSE)
  if (!all(c("sample", "source") %in% names(tbl))) {
    stop("--theta_lift_table must be a CSV with columns: sample,source")
  }
  unknown_sources <- setdiff(unique(tbl$source), reference$source)
  if (length(unknown_sources) > 0) {
    stop(sprintf(
      "--theta_lift_table references unknown source(s): %s. Available: %s",
      paste(unknown_sources, collapse = ", "), paste(reference$source, collapse = ", ")
    ))
  }
  for (i in seq_len(nrow(tbl))) {
    if (tbl$sample[i] %in% run_samples) sample_source[[tbl$sample[i]]] <- tbl$source[i]
  }
}

## Resolve theta value per sample, and fail loudly on any unfilled ones
per_sample_theta <- list()
per_sample_provenance <- list()

for (s in run_samples) {
  src <- sample_source[[s]]
  row <- dplyr::filter(reference, source == src)
  theta_val <- suppressWarnings(as.numeric(row$theta[1]))
  if (is.na(theta_val)) {
    stop(sprintf(
      "Reference theta for source '%s' (needed by sample '%s') is not a valid number in %s - it looks like a placeholder that hasn't been filled in yet.",
      src, s, opt$reference_csv
    ))
  }
  per_sample_theta[[s]] <- theta_val
  per_sample_provenance[[s]] <- list(
    source = src, achieved_fpr = row$achieved_fpr[1], n_loci = row$n_loci[1],
    calibration_date = row$calibration_date[1], laser_version = row$laser_version[1]
  )
}

result <- list(
  theta_mode            = "lift",
  theta                 = NA,   # no single global value in lift mode - always by sample
  per_sample_theta      = per_sample_theta,
  per_sample_provenance = per_sample_provenance,
  default_source        = opt$default_source,
  reference_asset       = opt$reference_csv,
  seed                  = NA,
  timestamp             = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)

jsonlite::write_json(result, opt$out, auto_unbox = TRUE, pretty = TRUE, na = "null")

cat(sprintf(
  "Lifted theta for %d sample(s). Sources used: %s\n",
  length(run_samples),
  paste(unique(unlist(sample_source)), collapse = ", ")
))
