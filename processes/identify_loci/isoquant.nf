nextflow.enable.dsl = 2

process ISOQUANT {
    tag "IsoQuant (${labels.size()} samples)"

    publishDir "${params.outdir}/preprocessing/isoquant", mode: 'link'

    conda "${projectDir}/envs_prebuilt/isoquant"

    input:
    tuple val(labels),
          path(bams),
          path(bais)
    path  reference
    path  annotation
    val   stranded
    val   data_type

    output:
    path "out/bulk/bulk.transcript_models.gtf",         emit: transcript_models_gtf
    path "out/bulk/bulk.transcript_counts.tsv",         emit: combined_counts
    path "out/bulk/bulk.gene_counts.tsv",               emit: combined_gene_counts
    path "out/bulk/bulk.transcript_model_reads.tsv.gz", emit: model_reads
    path "out/bulk/bulk.sqanti_classification.tsv",     emit: sqanti_classification
    path "out/bulk/bulk.read_assignments.tsv.gz",       emit: read_assignments, optional: true

    script:
    def bam_str   = (bams   instanceof List ? bams   : [bams]).collect { it.toString() }.join(' ')
    def label_str = (labels instanceof List ? labels : [labels]).collect { it.toString() }.join(' ')

    """
    for bam in ${bam_str}; do
        [ -f "\${bam}.bai" ] || samtools index -@ ${task.cpus} "\${bam}"
    done

    isoquant \\
        --data_type                  ${data_type}        \\
        --reference                  ${reference}        \\
        --genedb                     ${annotation}       \\
        --complete_genedb                                \\
        --bam                        ${bam_str}          \\
        --labels                     ${label_str}        \\
        --stranded                   ${stranded}         \\
        --output                     out                 \\
        --prefix                     bulk                \\
        --threads                    ${task.cpus}        \\
        --transcript_quantification  all                 \\
        --gene_quantification        all                 \\
        --matching_strategy          loose               \\
        --normalization_method       usable_reads        \\
        --count_exons                                    \\
        --sqanti_output                                  \\
        --min_mapq                   0                   \\
        --inconsistent_mapq_cutoff   0

    # IsoQuant names this file transcript_models_SQANTI.tsv - rename to match the expected output name used downstream
    if [ -f out/bulk/bulk.transcript_models_SQANTI.tsv ]; then
        mv out/bulk/bulk.transcript_models_SQANTI.tsv out/bulk/bulk.sqanti_classification.tsv
    else
        touch out/bulk/bulk.sqanti_classification.tsv
    fi

    # Guard for model reads in case of sparse/targeted libraries
    if [ ! -f out/bulk/bulk.transcript_model_reads.tsv.gz ]; then
        printf '' | gzip > out/bulk/bulk.transcript_model_reads.tsv.gz
    fi
    """

    stub:
    """
    mkdir -p out/bulk
    touch out/bulk/bulk.transcript_models.gtf
    touch out/bulk/bulk.transcript_counts.tsv
    touch out/bulk/bulk.gene_counts.tsv
    printf '' | gzip > out/bulk/bulk.transcript_model_reads.tsv.gz
    touch out/bulk/bulk.sqanti_classification.tsv
    touch out/bulk/bulk.read_assignments.tsv.gz
    """
}
