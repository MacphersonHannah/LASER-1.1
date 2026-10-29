process SUCCESSFUL_LOCI {
  tag { "aggregate" }
  publishDir "${params.outdir}/identify_loci", mode: 'copy', overwrite: true

  input:
    path loci_csvs

  output:
    path "successful_loci.all.csv", emit: all_loci
    path "successful_loci.pass.csv", emit: pass_loci

  conda "${projectDir}/envs_prebuilt/core_utils_env"

  script:
  """
  set -euo pipefail

  python3 ${projectDir}/scripts/merge_loci_status.py \\
    --out-all successful_loci.all.csv \\
    --out-pass successful_loci.pass.csv \\
    ${loci_csvs.join(' ')}
  """
}
