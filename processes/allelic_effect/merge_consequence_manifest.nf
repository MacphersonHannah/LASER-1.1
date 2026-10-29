process MERGE_CONSEQUENCE_MANIFEST {
  publishDir "${params.outdir}/allelic_effect_testing/overview", mode: 'link'
  conda "${projectDir}/envs_prebuilt/stats_r_env"

  input:
    path manifests, stageAs: 'manifest_input*/*'

  output:
    path "consequence_test_summary.tsv", emit: manifest_tsv

  script:
  def manifest_list = manifests instanceof List ? manifests.join(',') : manifests
  def run_flag = params.conseq_run_label ? "--run_label ${params.conseq_run_label}" : ''
  """
  set -euo pipefail
  Rscript ${projectDir}/scripts/allelic_effect/merge_consequence_manifest.R \\
    --manifests ${manifest_list} \\
    ${run_flag} \\
    --out consequence_test_summary.tsv
  """
}