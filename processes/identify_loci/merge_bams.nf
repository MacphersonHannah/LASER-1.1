// Merges all per-sample sorted BAMs into a single BAM for FLAIR transcriptome.
// FLAIR transcriptome is run once on the combined BAM to build a consensus
// transcriptome across all samples, then flair quantify runs per-sample.
nextflow.enable.dsl = 2

process MERGE_BAMS {

    tag { 'all_samples' }
    publishDir "${params.outdir}/preprocessing/flair/merged_bam", mode: 'link'

    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    input:
        path(bams)    // collected list of all sample BAMs

    output:
        path("all_samples_merged.bam"),     emit: bam
        path("all_samples_merged.bam.bai"), emit: bai

    script:
    """
    samtools merge -f -@ ${task.cpus} all_samples_merged.bam ${bams}
    samtools index all_samples_merged.bam
    """
}
