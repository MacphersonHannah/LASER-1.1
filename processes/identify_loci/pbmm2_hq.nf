// just for HQ reads
process PBMM2_HQ {
    tag { sample }
    publishDir "${params.outdir}/preprocessing/pbmm2_hq/${sample}", mode: 'link'

    input:
        tuple val(sample), path(hq_bam), val(mode)
        path genome

    output:
        tuple val(sample), path("${sample}_pbmm2.sorted.bam")

    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    script:
    """
    pbmm2 align ${genome} ${hq_bam} ${sample}_pbmm2.sorted.bam --preset ISOSEQ --sort -j ${task.cpus}
    """
}