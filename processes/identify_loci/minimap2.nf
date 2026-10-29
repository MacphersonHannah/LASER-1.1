process MINIMAP2 {

    tag { sample }
    publishDir "${params.outdir}/preprocessing/minimap2/${sample}", mode: 'link'

    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    input:
        tuple val(sample), path(fastq)
        val(minimap_args)

    output:
        tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam

    script:
    """
    minimap2 \\
        ${minimap_args} \\
        --secondary=no \\
        --MD \\
        -G 500k \\
        -t ${task.cpus} \\
        ${params.genome} \\
        ${fastq} \\
    | samtools sort -@ ${task.cpus} -m 2G -o ${sample}.sorted.bam

    samtools index ${sample}.sorted.bam
    """
}
