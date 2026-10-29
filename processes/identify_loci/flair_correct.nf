// Corrects misaligned splice sites using GTF annotation.
// Novel splice sites will not be detected without orthogonal short-read support.

process FLAIR_CORRECT {

    tag { sample }
    publishDir "${params.outdir}/preprocessing/flair/correct/${sample}", mode: 'link'

    conda "${projectDir}/envs_prebuilt/flair_env"

    input:
        tuple val(sample), path(bed)
        path(genome)
        path(gtf)

    output:
        tuple val(sample), path("${sample}_all_corrected.bed"), emit: corrected_bed

    script:
    """
    flair correct \\
        --query ${bed} \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        --output ${sample} \\
        --threads ${task.cpus}
    """
}
