process RUN_CONSEQUENCE_TEST {
  tag { sample }
  publishDir "${params.outdir}/allelic_effect_testing/${sample}", mode: 'link'
  conda "${projectDir}/envs_prebuilt/stats_r_env"


  input:
    tuple val(sample), path(sample_rds), path(allelic_results), path(theta_json), path(transcript_features)

  output:
    tuple val(sample), path("consequence_results.tsv"), emit: results_tsv
    path "consequence_manifest.tsv", emit: sample_manifest

  script:
  def opt_flags = [
    conseq_perm_b:         '--conseq_perm_b',
    conseq_perm_seed:      '--conseq_perm_seed',
    conseq_fdr_threshold:  '--conseq_fdr_threshold',
    conseq_dist_test:      '--conseq_dist_test',
  ].findResults { key, flag -> params[key] != null ? "${flag} ${params[key]}" : null }
   .join(' \\\n    ')

  """
  
  set -euo pipefail
  Rscript ${projectDir}/scripts/allelic_effect/run_consequence_test.R \\
    --sample_rds          ${sample_rds} \\
    --allelic_results     ${allelic_results} \\
    --transcript_features ${transcript_features} \\
    --theta_json          ${theta_json} \\
    --sample              ${sample} \\
    ${opt_flags} \\
    --out_results  consequence_results.tsv \\
    --out_manifest consequence_manifest.tsv
  """
}