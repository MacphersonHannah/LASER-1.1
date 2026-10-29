#!/usr/bin/env Rscript
## theta_mode = 'fixed'. User gives a theta to every sample or a table. No computation needed just pass the values through like the other scripts do

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(jsonlite)
})

opt_list <- list(
  make_option("--theta_fixed", type = "character", default = "NA"),   # single numeric or "NA"
  make_option("--theta_table", type = "character", default = "NA"),   # path or "NA": CSV sample,theta
  make_option("--samples",     type = "character"),                  # comma-separated sample names in this run
  make_option("--out",         type = "character", default = "theta_fixed.json")
)
opt <- parse_args(OptionParser(option_list = opt_list))

run_samples <- strsplit(opt$samples, ",")[[1]]

has_global <- !is.na(opt$theta_fixed) && opt$theta_fixed != "NA" && nzchar(opt$theta_fixed)
has_table  <- !is.na(opt$theta_table) && opt$theta_table != "NA" && nzchar(opt$theta_table)

if (!has_global && !has_table) {
  stop(
    "theta_mode = 'fixed' requires at least one of --theta_fixed (single value) ",
    "or --theta_table (CSV: sample,theta). Neither was supplied."
  )
}

per_sample_theta <- NULL

if (has_table) {
  tbl <- readr::read_csv(opt$theta_table, show_col_types = FALSE)
  if (!all(c("sample", "theta") %in% names(tbl))) {
    stop("--theta_table must be a CSV with columns: sample,theta")
  }
  per_sample_theta <- setNames(as.list(tbl$theta), tbl$sample)

  missing_samples <- setdiff(run_samples, names(per_sample_theta))
  if (length(missing_samples) > 0 && !has_global) {
    stop(sprintf(
      "Samples missing from --theta_table and no --theta_fixed fallback set: %s",
      paste(missing_samples, collapse = ", ")
    ))
  }
}

global_theta <- if (has_global) as.numeric(opt$theta_fixed) else NA_real_

result <- list(
  theta_mode       = "fixed",
  theta            = if (has_global) global_theta else NA,
  per_sample_theta = per_sample_theta,   # NULL if no table given -> global applies to all
  seed             = NA,
  timestamp        = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)

jsonlite::write_json(result, opt$out, auto_unbox = TRUE, pretty = TRUE, na = "null")

cat(sprintf(
  "Fixed theta resolved: global=%s, per-sample overrides for %d sample(s).\n",
  if (has_global) as.character(global_theta) else "none",
  if (has_table) length(per_sample_theta) else 0
))
