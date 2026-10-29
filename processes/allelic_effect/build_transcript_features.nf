process BUILD_TRANSCRIPT_FEATURES {
  publishDir "${params.outdir}/allelic_effect_testing/overview", mode: 'link'
  conda "${projectDir}/envs_prebuilt/stats_r_env"

  input:
    path orf_summary

  output:
    path "transcript_features.tsv", emit: features_tsv

  script:
  """
  set -euo pipefail
  Rscript ${projectDir}/scripts/allelic_effect/build_transcript_features.R \\
    --orf_summary ${orf_summary} \\
    --out transcript_features.tsv
  """
}