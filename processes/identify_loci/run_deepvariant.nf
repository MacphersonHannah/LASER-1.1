// Runs DeepVariant on the prepared BAM files to call variants. Therefore is filtered by prepare_dv step.
// This makes it quicker but means that the VCFs are only for regions with sufficient coverage.

process RUN_DEEPVARIANT {
  tag { sample }
  scratch false
  publishDir "${params.outdir}/preprocessing/run_deepvariant/${sample}", mode: 'link'
  input:
    tuple val(sample), path(bam), path(bai), path(regions)
    tuple path(genome), path(genome_fai)
  output:
    tuple val(sample), path("DeepVariant/${sample}/${sample}.deepvariant.vcf.gz"), emit: vcf
    path  "DeepVariant/${sample}/${sample}.deepvariant.vcf.gz.tbi"
    path  "DeepVariant/${sample}/${sample}.deepvariant.g.vcf.gz"
    path  "DeepVariant/${sample}/${sample}.deepvariant.g.vcf.gz.tbi"
  container params.dv_container
  shell:
  '''
  set -euo pipefail
  export TMPDIR=$PWD
  outdir="DeepVariant/!{sample}"
  mkdir -p "$outdir"
  region_flag=""
  if [ -s "!{regions}" ]; then
    region_flag="--regions !{regions}"
  fi
  run_deepvariant \
    --model_type=!{ params.dv_model ?: 'MASSEQ' } \
    --ref !{genome} \
    --reads !{bam} \
    --output_vcf "$outdir/!{sample}.deepvariant.vcf.gz" \
    --output_gvcf "$outdir/!{sample}.deepvariant.g.vcf.gz" \
    --num_shards !{task.cpus} \
    $region_flag \
    --intermediate_results_dir "$outdir/intermediate_results_dir"
  '''
}
