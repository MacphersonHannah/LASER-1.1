process RUN_ALLELIC_EFFECT_TEST {
  tag { sample }
  publishDir "${params.outdir}/allelic_effect_testing/${sample}", mode: 'link'
  conda "${projectDir}/envs_prebuilt/stats_r_env"

  input:
    tuple val(sample), path(sample_rds), path(theta_json)

  output:
    tuple val(sample), path("allelic_effect_results.tsv"), emit: results_tsv
    path "sample_manifest_${sample}.tsv", emit: sample_manifest

  script:
  """
  set -euo pipefail
  Rscript ${projectDir}/scripts/allelic_effect/run_allelic_effect_test.R \\
    --sample_rds ${sample_rds} \\
    --theta_json ${theta_json} \\
    --sample ${sample} \\
    --abte_fdr_threshold        ${params.abte_fdr_threshold} \\
    --abte_effect_guard         ${params.abte_effect_guard} \\
    --abte_major_frac_threshold ${params.abte_major_frac_threshold} \\
    --abte_fdr_scope            ${params.abte_fdr_scope} \\
    --abtu_fdr_threshold        ${params.abtu_fdr_threshold} \\
    --abtu_effect_guard         ${params.abtu_effect_guard} \\
    --abtu_tvd_threshold        ${params.abtu_tvd_threshold} \\
    --abtu_perm_b               ${params.abtu_perm_b} \\
    --dm_perm_seed              ${params.dm_perm_seed} \\
    --out_results allelic_effect_results.tsv \\
    --out_manifest sample_manifest_${sample}.tsv
  """
}
