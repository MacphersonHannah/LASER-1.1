// Takes all genes from the bed created from the input gtf as 'total_genes'. Identifies which of these have at least one FLNC read mapped to them
// in the sample's FLNC bam file, outputs 'expressed_genes' bed and a tsv with counts per gene.
// Useful for baseline comparisons of what it is possible to do ASE on with the given data.

process COUNT_EXPRESSED_GENES {
    tag "${sample}"
    publishDir "${params.outdir}/identify_loci/count_genes/${sample}", mode: 'link'
    conda "${projectDir}/envs_prebuilt/isoform_collapse_env"
    input:
    tuple val(sample), path(flnc_bam), path(gene_bed)

    output:
    tuple val(sample), path("${sample}.gene_expression.tsv")

    script:
    """
    bash ${projectDir}/scripts/count_expressed_genes.sh \
        ${gene_bed} \
        ${flnc_bam} \
        ${sample}.gene_expression.tsv
    """
}
