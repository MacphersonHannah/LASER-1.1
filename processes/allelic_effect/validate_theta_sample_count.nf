process VALIDATE_THETA_SAMPLE_COUNT {
  tag "theta_sample_count_check"
  publishDir "${params.outdir}/allelic_effect_testing/_validation", mode: 'link'
  conda "${projectDir}/envs_prebuilt/core_utils_env"
  
  input:
    val sample_list   // comma joined string of every sample name in this run

  output:
    path "theta_sample_count_ok.txt", emit: ok

  script:
  """
  set -euo pipefail
  IFS=',' read -ra SAMPLES <<< "${sample_list}"
  n=\${#SAMPLES[@]}

  if [ "\$n" -lt 2 ]; then
    echo "theta_mode = 'calibrate_pooled' needs at least 2 samples to pool a theta across (found \$n: ${sample_list})." >&2
    echo "A single sample has nothing to pool against - use theta_mode = 'fixed' (supply your own value) or 'lift' (reuse a previously calibrated one) instead." >&2
    exit 1
  fi

  if [ "\$n" -eq 2 ]; then
    echo "WARNING: theta_mode = 'calibrate_pooled' is running on only 2 samples. The calibrated theta rests on a very small basis and may not generalise well - proceeding, but consider 'lift' or adding more samples if that matters for your analysis." >&2
  fi

  echo "\${n} samples: ${sample_list}" > theta_sample_count_ok.txt
  """
}
