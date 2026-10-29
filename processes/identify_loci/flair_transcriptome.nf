// Combined transcriptome from merged BAM of all samples.
// 3x faster and 20x less memory than flair correct + collapse.
// --stringent removed: asymmetric read loss between alleles creates false ASE signal.
// Direct RNA: see flair_transcriptome_drna.nf
nextflow.enable.dsl = 2

process FLAIR_TRANSCRIPTOME {

    tag { 'all_samples' }
    publishDir "${params.outdir}/preprocessing/flair/transcriptome", mode: 'link'

    conda "${projectDir}/envs_prebuilt/flair_env"

    input:
        path(merged_bam)
        path(merged_bai)
        path(genome)
        path(gtf)

    output:
        path("flair.isoforms.bed"), emit: isoforms_bed
        path("flair.isoforms.gtf"), emit: isoforms_gtf
        path("flair.isoforms.fa"),  emit: isoforms_fa

    script:
    def max_ends = params.flair_max_ends    ?: 100
    def support  = params.flair_min_support ?: 2
    """
    flair transcriptome \\
        --genomealignedbam ${merged_bam} \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        --output flair \\
        --check_splice \\
        --end_window 100 \\
        --max_ends ${max_ends} \\
        --support ${support} \\
        --threads ${task.cpus}
    """
}
