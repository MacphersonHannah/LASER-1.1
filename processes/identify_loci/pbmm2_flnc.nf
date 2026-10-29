// Maps FLNC reads to the reference genome using pbmm2. Previously minimap2.
// This is slower than HQ but necessary to retain all FLNC reads for ASE.
// Runs on each individual sample but not bulk.

process PBMM2_FLNC {
  tag { sample }
  publishDir "${params.outdir}/preprocessing/pbmm2_flnc/${sample}", mode: 'link'

  input:
    tuple val(sample), path(flnc_bam)
    path genome

  output:
    tuple val(sample), path("${sample}.flnc.mapped.bam"), emit: flnc_mapped

  conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

  shell:
  '''
  set -euo pipefail

  bam="!{sample}.flnc.mapped.bam"

  pbmm2 align \
      --preset ISOSEQ \
      --sort \
      --num-threads !{task.cpus} \
      !{genome} \
      !{flnc_bam} \
      "$bam"
  '''
}



