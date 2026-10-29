// This process computes ASE counts and performs statistical tests using a provided SQLite database of read-level annotations.

process ASE_COUNTS_AND_TESTS {
  publishDir "${params.outdir}/ase_by_counts", mode: 'copy', overwrite: true

  input:
    path sqlite
    path passed_success_csv
    path gene_bed

  output:
    path "counts.all_samples.tsv", emit: all_samples

  conda "${projectDir}/envs_prebuilt/ase_laser"

  script:
  """
  set -euo pipefail
  echo "DEBUG sqlite: ${sqlite}"
  echo "DEBUG passes: ${passed_success_csv}"
  echo "DEBUG gene_bed: ${gene_bed}"
  ls -l

  python3 -u ${projectDir}/scripts/ase_counts_sqlite.py \
    --sqlite "${sqlite}" \
    --passes "${passed_success_csv}" \
    --gene_bed "${gene_bed}" \
    --out "counts.all_samples.tsv"
  """
}
