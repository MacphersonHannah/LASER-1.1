#!/usr/bin/env Rscript
## Build per-transcript (bulk_pbid) feature table for the consequence module, from bulk ORFannotate_summary.tsv. Run once per run; features are
## independent of sample since ORFannotate is run on all samples together
##
## Emits transcript_features.tsv keyed by bulk_pbid:
##   bulk_pbid, tss, tes, orf_length, is_coding, is_nmd
##
## If column names in orfannotate change, these will need changing:
COL_PBID       <- "transcript_id"
COL_STRAND     <- "strand"
COL_TX_START   <- "transcript_start"
COL_TX_END     <- "transcript_end"
COL_ORF_LENGTH <- "orf_nt_len"
COL_CODING     <- "coding_class"
COL_NMD        <- "NMD_sensitive"

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
})

opt_list <- list(
  make_option("--orf_summary", type = "character"),
  make_option("--out", type = "character", default = "transcript_features.tsv"),
  make_option("--coding_true_values", type = "character",
              default = "coding,protein_coding,TRUE,True,true,1,yes,Y"),
  make_option("--nmd_true_values",    type = "character",
              default = "TRUE,True,true,1,yes,Y,NMD,nmd,nmd_sensitive")
)
opt <- parse_args(OptionParser(option_list = opt_list))

raw <- readr::read_tsv(opt$orf_summary, show_col_types = FALSE,
                       guess_max = 1e6)

## pull a column by name if present, else NA
pull_col <- function(df, nm) if (nm %in% names(df)) df[[nm]] else rep(NA, nrow(df))

to_logical <- function(x, true_vals) {
  tv <- str_split(true_vals, ",")[[1]] |> str_trim()
  out <- as.character(x) %in% tv
  out[is.na(x)] <- NA
  out
}

if (!(COL_PBID %in% names(raw)))
  stop(sprintf("pbid column '%s' not found in %s. Columns present: %s",
               COL_PBID, opt$orf_summary, paste(names(raw), collapse = ", ")))

tx_start <- suppressWarnings(as.numeric(pull_col(raw, COL_TX_START)))
tx_end   <- suppressWarnings(as.numeric(pull_col(raw, COL_TX_END)))
strand   <- as.character(pull_col(raw, COL_STRAND))

# strand-aware boundaries: TSS = 5' end, TES = 3' end in transcription direction
tss <- ifelse(strand == "-", tx_end,   tx_start)
tes <- ifelse(strand == "-", tx_start, tx_end)

features <- tibble(
  bulk_pbid = as.character(pull_col(raw, COL_PBID)),
  tss       = tss,
  tes       = tes,
  orf_length = suppressWarnings(as.numeric(pull_col(raw, COL_ORF_LENGTH))),
  is_coding = to_logical(pull_col(raw, COL_CODING), opt$coding_true_values),
  is_nmd    = to_logical(pull_col(raw, COL_NMD),    opt$nmd_true_values)
) %>%
  filter(!is.na(bulk_pbid), bulk_pbid != "") %>%
  distinct(bulk_pbid, .keep_all = TRUE)

readr::write_tsv(features, opt$out)

cat(sprintf(
  "transcript_features: %d transcripts | tss %d | tes %d | orf_length %d | coding %d | nmd %d non-NA\n",
  nrow(features),
  sum(!is.na(features$tss)), sum(!is.na(features$tes)), sum(!is.na(features$orf_length)),
  sum(!is.na(features$is_coding)), sum(!is.na(features$is_nmd))
))