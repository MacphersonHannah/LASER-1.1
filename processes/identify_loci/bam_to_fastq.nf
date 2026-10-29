process BAM_TO_FASTQ {
    tag "${sample}"
    
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.fastq.gz"), emit: fastq

    script:
    """
    samtools fastq ${bam} | bgzip > ${sample}.fastq.gz
    """
}
