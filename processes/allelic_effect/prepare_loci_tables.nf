process PREPARE_LOCI_TABLES {
  tag "allelic_effect_prep"
  publishDir "${params.outdir}/allelic_effect_testing/_prepared", mode: 'link'
  conda "${projectDir}/envs_prebuilt/stats_r_env"

  input:
    path read_counts_csv
    path gene_loci_bed

  output:
    path "combined_annotated.rds",          emit: combined_annotated_rds
    path "per_sample/*.loci.rds",           emit: per_sample_rds

  script:
  """

  
  set -euo pipefail
  mkdir -p per_sample
  Rscript ${projectDir}/scripts/allelic_effect/prepare_loci_tables.R \\
    --counts ${read_counts_csv} \\
    --bed ${gene_loci_bed} \\
    --min_transcripts_per_locus  ${params.min_transcripts_per_locus} \\
    --min_reads_per_transcript   ${params.min_reads_per_transcript} \\
    --min_total_reads_per_locus  ${params.min_total_reads_per_locus} \\
    --min_phased_reads_per_locus ${params.min_phased_reads_per_locus} \\
    --min_reads_per_haplotype    ${params.min_reads_per_haplotype} \\
    --max_unphased_prop          ${params.max_unphased_prop} \\
    --out_combined combined_annotated.rds \\
    --out_sample_dir per_sample
  """
}
