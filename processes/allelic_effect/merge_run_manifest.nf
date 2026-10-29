process MERGE_RUN_MANIFEST {
  tag "run_manifest"
  publishDir "${params.outdir}/allelic_effect_testing/overview", mode: 'link'
  conda "${projectDir}/envs_prebuilt/core_utils_env"

  input:
    path sample_manifests
    path theta_json

  output:
    path "allelic_effect_summary.tsv", emit: manifest_tsv

  script:
  """
  set -euo pipefail
  head -n 1 \$(ls sample_manifest*.tsv | head -n 1) > allelic_effect_summary.tsv
  for f in sample_manifest*.tsv; do
    tail -n +2 "\$f" >> allelic_effect_summary.tsv
  done

  {
    echo ""
    echo "# --- run-level params ---"
    echo "# theta_mode=${params.theta_mode}"
    echo "# theta_json=\$(cat ${theta_json})"
    echo "# abte_fdr_threshold=${params.abte_fdr_threshold}"
    echo "# abte_effect_guard=${params.abte_effect_guard}"
    echo "# abte_major_frac_threshold=${params.abte_major_frac_threshold}"
    echo "# abte_fdr_scope=${params.abte_fdr_scope}"
    echo "# abtu_fdr_threshold=${params.abtu_fdr_threshold}"
    echo "# abtu_effect_guard=${params.abtu_effect_guard}"
    echo "# abtu_tvd_threshold=${params.abtu_tvd_threshold}"
    echo "# abtu_perm_b=${params.abtu_perm_b}"
    echo "# min_transcripts_per_locus=${params.min_transcripts_per_locus}"
    echo "# min_reads_per_transcript=${params.min_reads_per_transcript}"
    echo "# min_total_reads_per_locus=${params.min_total_reads_per_locus}"
    echo "# min_phased_reads_per_locus=${params.min_phased_reads_per_locus}"
    echo "# min_reads_per_haplotype=${params.min_reads_per_haplotype}"
    echo "# max_unphased_prop=${params.max_unphased_prop}"
  } >> allelic_effect_summary.tsv
  """
}
