// This is based on an old script from cDNA cupcake

process GET_ABUNDANCE {
    tag { sample }
    publishDir "${params.outdir}/preprocessing/Collapsed_isoforms/${sample}", mode: 'link'
    input:
        tuple val(sample), val(prefix), path(cluster_report)
    output:
        path "Collapsed_isoforms/${sample}/${sample}.collapsed.abundance.txt"
        tuple val(sample), path("Collapsed_isoforms/${sample}/${sample}_mapped_fl_count.txt")
        path "Collapsed_isoforms/${sample}/${sample}_get_abundance.log"
    conda "${projectDir}/envs_prebuilt/isoform_collapse_env"
    script:
    """
    mkdir -p Collapsed_isoforms/${sample}

    python3 ${projectDir}/scripts/get_abundance_post_collapse.py $prefix $cluster_report \
        > Collapsed_isoforms/${sample}/${sample}_get_abundance.log
    """

}
