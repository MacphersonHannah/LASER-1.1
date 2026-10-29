process RESOLVE_FIXED_THETA {
  tag "theta_fixed"
  publishDir "${params.outdir}/allelic_effect_testing/_theta", mode: 'link'
  conda "${projectDir}/envs_prebuilt/stats_r_env"

  input:
    val sample_list   // comma joined string of every sample name in this run

  output:
    path "theta_fixed.json", emit: theta_value

  script:
  def theta_fixed_arg = params.theta_fixed ?: 'NA'
  def theta_table_arg = params.theta_table ?: 'NA'
  """
  set -euo pipefail
  Rscript ${projectDir}/scripts/allelic_effect/resolve_fixed_theta.R \\
    --theta_fixed "${theta_fixed_arg}" \\
    --theta_table "${theta_table_arg}" \\
    --samples "${sample_list}" \\
    --out theta_fixed.json
  """
}
