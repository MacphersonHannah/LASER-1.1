process CALIBRATE_THETA_POOLED {
  tag "theta_calibration_pooled"
  publishDir "${params.outdir}/allelic_effect_testing/_theta", mode: 'link'
  conda "${projectDir}/envs_prebuilt/stats_r_env"

  input:
    path combined_annotated_rds
    path validation_marker   // from VALIDATE_THETA_GROUP_SAMPLE_COUNT - unused in script, present only
                              // to force this process to wait until the per group sample count check passes

  output:
    path "theta_pooled.json", emit: theta_value

  script:
  """
  set -euo pipefail
  Rscript ${projectDir}/scripts/allelic_effect/calibrate_theta_pooled.R \\
    --combined_rds ${combined_annotated_rds} \\
    --max_loci    ${params.theta_calib_max_loci} \\
    --n_null_ref  ${params.theta_calib_n_null_ref} \\
    --n_test      ${params.theta_calib_n_test} \\
    --nominal_fpr ${params.theta_calib_nominal_fpr} \\
    --seed        ${params.theta_calib_seed} \\
    --out theta_pooled.json
  """
}
