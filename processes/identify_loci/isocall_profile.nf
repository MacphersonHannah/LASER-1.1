process ISOCALL_PROFILE {
    tag { sample }
    publishDir "${params.outdir}/preprocessing/isocall/profiles", mode: 'link'

    input:
        tuple val(sample), path(flnc_mapped_bam)

    output:
        path "${sample}.profile.gz", emit: profile

    conda "${projectDir}/envs_prebuilt/isocall_laser"

    script:
    """
    ${params.isocall_bin} profile \\
        --reads ${flnc_mapped_bam} \\
        --sample ${sample} \\
        --output ${sample}.profile.gz
    """
}
