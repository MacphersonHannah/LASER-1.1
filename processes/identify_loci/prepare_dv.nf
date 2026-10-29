// Prepares FLNC BAM for DeepVariant: creates regions BED and sorts/indexes the BAM.
process PREPARE_DV_BAM {
  tag { sample }
  publishDir "${params.outdir}/preprocessing/prepare_dv_bam/${sample}", mode: 'link'
  conda "${projectDir}/envs_prebuilt/isoform_collapse_env"

  input:
    tuple val(sample), path(flnc_mapped_bam)
    path gene_bed

  output:
    tuple val(sample),
          path("${sample}.hq.filtered.sorted.bam"),
          path("${sample}.hq.filtered.sorted.bam.bai"),
          path("${sample}.filtered_genes.bed"),
          emit: bam

  script:
  """
  set -euo pipefail

  # Reduce gene_bed to 4 columns for DeepVariant regions
  awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4}' "${gene_bed}" > "${sample}.filtered_genes.bed"

  # Check BAM sort order - pbmm2 --sort already sorts, so skip if coordinate-sorted
  sort_order=\$(samtools view -H "${flnc_mapped_bam}" | awk -F'\t' '/^@HD/ {for(i=1;i<=NF;i++) if(\$i ~ /^SO:/) {print substr(\$i,4); exit}}')
  echo "[DV-PREP] ${sample}: input BAM sort order: \${sort_order:-unknown}"
  if [ "\$sort_order" = "coordinate" ]; then
    ln -sf "\$(readlink -f "${flnc_mapped_bam}")" "${sample}.hq.filtered.sorted.bam"
  else
    echo "[DV-PREP] ${sample}: re-sorting (not coordinate-sorted)"
    samtools sort -@ ${task.cpus} -o "${sample}.hq.filtered.sorted.bam" "${flnc_mapped_bam}"
  fi
  samtools index -@ ${task.cpus} "${sample}.hq.filtered.sorted.bam"

  echo "[DV-PREP] ${sample}: \$(samtools view -c ${sample}.hq.filtered.sorted.bam) reads"
  """
}