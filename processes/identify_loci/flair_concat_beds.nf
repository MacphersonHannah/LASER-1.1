// Concatenates corrected BED12 files before FLAIR collapse.
// FLAIR collapse requires all samples in a single BED file to build
// one consensus transcriptome across all samples.
nextflow.enable.dsl = 2

process FLAIR_CONCAT_BEDS {

    tag { 'all_samples' }
    publishDir "${params.outdir}/preprocessing/flair/concat", mode: 'link'

    conda "${projectDir}/envs_prebuilt/flair_env"

    input:
        path(beds)   

    output:
        path("all_samples_corrected.bed"), emit: combined_bed

    script:
    """
    cat ${beds} > all_samples_corrected.bed
    """
}
