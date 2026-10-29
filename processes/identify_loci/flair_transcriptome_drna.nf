// FLAIR transcriptome for direct RNA sequencing.
// Key differences from cDNA version:
//   - No --stringent: direct RNA reads are fragmented/degraded at 5' end
//   - No --check_splice: splice site coverage unreliable with higher error rate
//   - No --trust_ends: 5' ends are unreliable (degradation); 3' ends are more reliable
//     but FLAIR uses both for TSS/TES - safer to let it use the data cautiously
//   - --annotation_reliant generate: anchors to known transcripts first to
//     reduce degradation products being called as novel isoforms
//   - flair_min_support default 5 (vs 3 for cDNA) to reduce noise
nextflow.enable.dsl = 2

process FLAIR_TRANSCRIPTOME_DRNA {

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
    def max_ends = params.flair_max_ends    ?: 1000
    def support  = params.flair_min_support ?: 2
    """
    flair transcriptome \\
        --genomealignedbam ${merged_bam} \\
        --genome ${genome} \\
        --gtf ${gtf} \\
        --output flair \\
        --annotation_reliant generate \\
        --end_window 100 \\
        --max_ends ${max_ends} \\
        --support ${support} \\
        --threads ${task.cpus}
    """
}
