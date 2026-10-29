#!/usr/bin/env Rscript
## merge individual sample consequence manifests into one run-level summary table

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(purrr)
})

option_list <- list(
  make_option("--manifests", type = "character",
              help = "Comma-separated list of per-sample consequence_manifest.tsv paths"),
  make_option("--run_label", type = "character", default = NA_character_,
              help = "Optional tag (eg 'mouse', 'brain') added as a column for cross-run merges"),
  make_option("--out", type = "character",
              help = "Output path for merged consequence manifest")
)
opt <- parse_args(OptionParser(option_list = option_list))

manifest_paths <- strsplit(opt$manifests, ",")[[1]]

merged <- purrr::map_dfr(manifest_paths, readr::read_tsv, show_col_types = FALSE) %>%
  dplyr::arrange(sample)

if (!is.na(opt$run_label)) {
  merged <- merged %>% dplyr::mutate(run = opt$run_label, .before = sample)
}

readr::write_tsv(merged, opt$out)

cat(sprintf("[merge_consequence_manifest] merged %d sample manifests -> %s\n",
            nrow(merged), opt$out))