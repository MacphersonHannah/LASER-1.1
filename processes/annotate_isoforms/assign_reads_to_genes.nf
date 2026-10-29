// using featurecounts to define a gene

process ASSIGN_READS_TO_GENES {
  tag { sample }
  publishDir "${params.outdir}/annotate_isoforms/featurecounts/${sample}", mode: 'link'

  input:
    tuple val(sample), path(bam)
    path gtf

  output:
    tuple val(sample), path("${sample}.read2gene.tsv"), emit: read2gene

  conda "${projectDir}/envs_prebuilt/featurecounts_env"

  script:
  """
  set -euo pipefail

  featureCounts \\
    -L \\
    --largestOverlap \\
    -T ${task.cpus} \\
    -a ${gtf} \\
    -t exon -g gene_id \\
    --primary \\
    -s 1 \\
    -R CORE \\
    ${bam} \\
    -o ${sample}.featureCounts.txt

  # CORE output is <input_bam_basename>.featureCounts
  core="\$(basename ${bam}).featureCounts"
  if [ ! -f "\$core" ]; then
    echo "ERROR: CORE file \$core not found" >&2
    ls -lh >&2
    exit 1
  fi

  cat > parse_core.awk << 'AWK'
  BEGIN { OFS="\\t" }
  /^#/ { next }
  NF >= 4 {
    # read_id, status, n_targets, genes
    print \$1, \$2, \$3, \$4
  }
  AWK

  awk -F '\\t' -f parse_core.awk "\$core" > ${sample}.read2gene.tsv
  """
}
