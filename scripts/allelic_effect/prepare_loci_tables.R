#!/usr/bin/env Rscript
## Build per-locus, per-sample summary tables from LASER's combined read_counts_summary.csv + gene_loci_named.bed, and split into one
## table per sample. Every locus with at least one read is kept, tagged with a tested_status, rather than silently dropping untestable loci.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
})

opt_list <- list(
  make_option("--counts",  type = "character"),
  make_option("--bed",     type = "character"),
  make_option("--min_transcripts_per_locus",  type = "integer",   default = 2),
  make_option("--min_reads_per_transcript",   type = "integer",   default = 1),
  make_option("--min_total_reads_per_locus",  type = "integer",   default = 30),
  make_option("--min_phased_reads_per_locus", type = "integer",   default = 30),
  make_option("--min_reads_per_haplotype",    type = "integer",   default = 5),
  make_option("--max_unphased_prop",          type = "double",    default = 0.30),
  make_option("--out_combined", type = "character", default = "combined_annotated.rds"),
  make_option("--out_sample_dir", type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = opt_list))

gene_loci <- readr::read_tsv(
  opt$bed,
  col_names = c("chrom", "gene_start", "gene_end", "gene_name", "associated_gene", "strand"),
  show_col_types = FALSE
) %>%
  mutate(across(c(associated_gene, gene_name, chrom, strand), as.character)) %>%
  distinct(associated_gene, .keep_all = TRUE) %>%
  ## Fall back to the Ensembl ID whenever the GTF has no gene_name attribute
  ## (make_gene_bed.nf writes the literal string "NA" in that case, not a true NA).
  mutate(gene_name = if_else(is.na(gene_name) | gene_name %in% c("", "NA"), associated_gene, gene_name))

counts <- readr::read_csv(opt$counts, show_col_types = FALSE) %>%
  mutate(
    sample = as.character(sample),
    ## sample_short used to be derived by stripping everything after the
    ## first underscore (eg ENCFF422XLS_ovary -> ENCFF422XLS). But that assumed an ACCESSION_label naming convention from the test data that isn't guaranteed for other datasets. 
    ## so sample_short is now just the real sample name, unmodified. Kept as its own column (rather than
    ## removed) because every downstream join in calibrate_theta_pooled.R, calibrate_theta_pooled_by_group.R, and run_allelic_effect_test.R
    ## already joins on c("sample", "sample_short", ...); aliasing it here keeps every one of those joins valid with no other changes required.
    sample_short = sample,
    bulk_pbid = as.character(bulk_pbid),
    haplotype = as.character(haplotype),
    associated_gene = as.character(associated_gene),
    associated_transcript = as.character(associated_transcript),
    count = as.numeric(count)
  ) %>%
  filter(
    !is.na(bulk_pbid), bulk_pbid != "",
    !is.na(associated_gene), associated_gene != "", associated_gene != "NA",
    count > 0
  ) %>%
  left_join(gene_loci, by = "associated_gene")

## transcript then locus-level summary (kept for ALL loci)

tx_by_locus_sample <- counts %>%
  group_by(sample, sample_short,
           associated_gene, gene_name, chrom, gene_start, gene_end, strand, bulk_pbid) %>%
  summarise(
    tx_hap1_reads     = sum(count[haplotype == "1"],        na.rm = TRUE),
    tx_hap2_reads     = sum(count[haplotype == "2"],        na.rm = TRUE),
    tx_unphased_reads = sum(count[haplotype == "unphased"], na.rm = TRUE),
    tx_total_reads    = sum(count,                          na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tx_phased_reads = tx_hap1_reads + tx_hap2_reads,
    tx_detected     = tx_total_reads >= opt$min_reads_per_transcript
  )

locus_by_sample <- tx_by_locus_sample %>%
  group_by(sample, sample_short,
           associated_gene, gene_name, chrom, gene_start, gene_end, strand) %>%
  summarise(
    n_transcripts_any_reads = n_distinct(bulk_pbid[tx_total_reads > 0]),
    n_detected_transcripts  = n_distinct(bulk_pbid[tx_detected]),
    locus_hap1_reads        = sum(tx_hap1_reads,     na.rm = TRUE),
    locus_hap2_reads        = sum(tx_hap2_reads,     na.rm = TRUE),
    locus_unphased_reads    = sum(tx_unphased_reads, na.rm = TRUE),
    locus_phased_reads      = sum(tx_phased_reads,   na.rm = TRUE),
    locus_total_reads       = sum(tx_total_reads,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    locus_unphased_prop = if_else(locus_total_reads > 0,
                                   locus_unphased_reads / locus_total_reads, NA_real_),
    ## associated_gene came from SQANTI's sqanti_associated_gene call, which can assign novel/synthetic gene IDs (eg intergenic novelGene_* calls) that
    ## don't exist in the reference GTF used to build gene_loci_named.bed. Those rows survive the left_join with chrom/gene_start/gene_end left NA - catch
    ## that here rather than let them silently pass the read/transcript gates.
    pass_reference_gene = !is.na(chrom),
    pass_transcript_count = n_detected_transcripts >= opt$min_transcripts_per_locus,
    pass_read_count =
      locus_total_reads   >= opt$min_total_reads_per_locus  &
      locus_phased_reads  >= opt$min_phased_reads_per_locus &
      locus_hap1_reads    >= opt$min_reads_per_haplotype    &
      locus_hap2_reads    >= opt$min_reads_per_haplotype    &
      !is.na(locus_unphased_prop) &
      locus_unphased_prop <= opt$max_unphased_prop,
    tested_status = case_when(
      !pass_reference_gene                      ~ "not_tested_gene_not_in_reference",
      !pass_transcript_count                    ~ "not_tested_too_few_transcripts",
      pass_transcript_count & !pass_read_count   ~ "not_tested_insufficient_reads",
      TRUE                                       ~ "tested"
    )
  )

## Transcript-level detail needed by the abTU DM test, kept alongside
tx_input <- tx_by_locus_sample %>%
  select(sample, sample_short, associated_gene, bulk_pbid,
         tx_hap1_reads, tx_hap2_reads, tx_phased_reads)

saveRDS(
  list(locus_by_sample = locus_by_sample, tx_input = tx_input),
  opt$out_combined
)

##  split by sample for parallel per-sample testing

dir.create(opt$out_sample_dir, showWarnings = FALSE, recursive = TRUE)

for (s in unique(locus_by_sample$sample_short)) {
  saveRDS(
    list(
      locus_by_sample = filter(locus_by_sample, sample_short == s),
      tx_input        = filter(tx_input,        sample_short == s)
    ),
    file.path(opt$out_sample_dir, paste0(s, ".loci.rds"))
  )
}

cat(sprintf(
  "Prepared %d loci across %d samples (%d tested, %d not tested).\n",
  nrow(locus_by_sample),
  n_distinct(locus_by_sample$sample_short),
  sum(locus_by_sample$tested_status == "tested"),
  sum(locus_by_sample$tested_status != "tested")
))