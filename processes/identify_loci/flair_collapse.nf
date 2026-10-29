// High confidence isoform calling from all corrected reads.
// Run after FLAIR_CONCAT_BEDS with one combined BED and comma FASTQs.
// Direct RNA: remove --check_splice, --stringent, --trust_ends
nextflow.enable.dsl = 2

process FLAIR_COLLAPSE {

    tag { 'all_samples' }
    publishDir "${params.outdir}/flair/collapse", mode: 'link'

    conda "${projectDir}/envs_prebuilt/flair_env"

    input:
        path(combined_corrected_bed)   // from FLAIR_CONCAT_BEDS
        val(reads_csv)                 // comma FASTQ paths string
        path(genome)
        path(gtf)

    output:
        path("flair.isoforms.bed"),              emit: isoforms_bed
        path("flair.isoforms.gtf"),              emit: isoforms_gtf
        path("flair.isoforms.fa"),               emit: isoforms_fa
        path("flair.isoforms.read.map.txt"),     emit: read_map

    script:
    def max_ends = params.flair_max_ends    ?: 100
    def support  = params.flair_min_support ?: 2
    """
    flair collapse \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        --query ${combined_corrected_bed} \\
        --reads ${reads_csv} \\
        --output flair \\
        --check_splice \\
        --trust_ends \\
        --no_gtf_end_adjustment \\
        --annotation_reliant generate \\
        --generate_map \\
        --max_ends ${max_ends} \\
        --support ${support} \\
        --threads ${task.cpus}
    """
}
