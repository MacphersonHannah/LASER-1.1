// Prepares known isoforms from a reference GTF annotation for isocall.
// Runs once per pipeline execution (not per sample).
// Requires params.isocall_bin

process ISOCALL_PREP_ISOFORMS {
    tag { "reference" }
    publishDir "${params.outdir}/preprocessing/isocall/reference", mode: 'link'

    input:
        path gtf

    output:
        path "ref_isoforms.gz", emit: ref_isoforms

    conda "${projectDir}/envs_prebuilt/isocall_laser"

    script:
    """
    ${params.isocall_bin} prep-isoforms \
        --gtf ${gtf} \
        --output ref_isoforms.gz
    """
}
