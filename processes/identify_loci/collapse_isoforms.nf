process COLLAPSE_ISOFORMS {
    tag { sample }
    errorStrategy 'finish'
    maxRetries 1
    publishDir "${params.outdir}/preprocessing/Collapsed_isoforms/${sample}", mode: 'link'
    input:
        tuple val(sample), path(hq_bam)

    output:
        path("Collapsed_isoforms/${sample}/${sample}.collapsed.gff"),        emit: gff
        path("Collapsed_isoforms/${sample}/${sample}.collapsed.group.txt"), emit: group_txt
        path("Collapsed_isoforms/${sample}/${sample}_collapse_isoforms.log"), emit: log

    conda "${projectDir}/envs_prebuilt/isoform_collapse_env"

    script:
    """
    mkdir -p Collapsed_isoforms/${sample}
    # Run IsoSeq collapse
    isoseq collapse \
        --do-not-collapse-extra-5exons \
        --max-5p-diff 5 \
        --max-3p-diff 5 \
        ${hq_bam} \
        Collapsed_isoforms/${sample}/${sample}.collapsed.gff \
        > Collapsed_isoforms/${sample}/${sample}_collapse_isoforms.log 2>&1

    """
}

