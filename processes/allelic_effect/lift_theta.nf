process LIFT_THETA {
  tag "theta_lift"
  publishDir "${params.outdir}/allelic_effect_testing/_theta", mode: 'link'
  conda "${projectDir}/envs_prebuilt/stats_r_env"

  input:
    val sample_list       // comma joined string of every sample name in this run
    path reference_csv    // assets/calibrated_thetas_reference.csv - other speciified in params

  output:
    path "theta_lift.json", emit: theta_value

  script:
  def source_table_arg = params.theta_lift_table ?: 'NA'
  """
  set -euo pipefail
  Rscript ${projectDir}/scripts/allelic_effect/lift_theta.R \\
    --reference_csv "${reference_csv}" \\
    --default_source "${params.theta_lift_source}" \\
    --source_table "${source_table_arg}" \\
    --samples "${sample_list}" \\
    --out theta_lift.json
  """
}
