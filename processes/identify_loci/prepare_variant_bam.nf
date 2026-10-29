process PREPARE_VARIANT_BAM {
  tag { sample }
  publishDir "${params.outdir}/preprocessing/prepare_dv_bam/${sample}", mode: 'link'
  conda "${projectDir}/envs_prebuilt/isoform_collapse_env"

  input:
    tuple val(sample), path(bam), path(bai)
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
  set -x

  regions="${sample}.filtered_genes.bed"

  echo "[DV-PREP-ONT] sample=${sample}"
  ls -lh "${bam}" "${gene_bed}" || true
  echo "[DV-PREP-ONT] gene_bed head:"; head -n 3 "${gene_bed}" || true
  echo "[DV-PREP-ONT] mapped reads: \$(samtools view -c \"${bam}\")"

  # Reduce gene_bed to 4 columns - matches PacBio PREPARE_DV_BAM output format
  awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4}' "${gene_bed}" > "\$regions"

  # BAM is already sorted (minimap2 pipes through samtools sort, pbmm2 has --sort).
  # Skip re-sorting - for 50GB+ BAMs that saves 30-60 minutes per sample.
  # Verify sort order from header; only re-sort if header says otherwise.
  sort_order=\$(samtools view -H "${bam}" | awk -F'\\t' '/^@HD/ {for(i=1;i<=NF;i++) if(\$i ~ /^SO:/) {print substr(\$i,4); exit}}')
  echo "[DV-PREP-ONT] input BAM sort order: \${sort_order:-unknown}"

  if [ "\$sort_order" = "coordinate" ]; then
    # Use a symlink so the BAM is staged into the work dir without copying
    ln -sf "\$(readlink -f "${bam}")" "${sample}.hq.filtered.sorted.bam"
  else
    echo "[DV-PREP-ONT] input not coordinate-sorted, sorting now"
    samtools sort -@ ${task.cpus} -o "${sample}.hq.filtered.sorted.bam" "${bam}"
  fi
  samtools index -@ ${task.cpus} "${sample}.hq.filtered.sorted.bam"
  """
}