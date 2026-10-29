process ORFANNOTATE {
    tag { sample }
    publishDir "${params.outdir}/annotate_isoforms/ORFANNOTATE", mode: 'link'
    input:
        tuple val(sample), path(corrected_gtf)
        path genome
    output:
        tuple val(sample), path("ORFannotate/${sample}/ORFannotate_summary.tsv"), emit: summary
        tuple val(sample), path("ORFannotate/${sample}/ORFannotate_annotated.gtf"), emit: annotated_gtf
        tuple val(sample), path("ORFannotate/${sample}/uORFannotate_annotated.gtf"), emit: uorf_annotated_gtf, optional: true
        tuple val(sample), path("ORFannotate/${sample}/utr5.fa"), emit: utr5_fa, optional: true
        tuple val(sample), path("ORFannotate/${sample}/utr3.fa"), emit: utr3_fa, optional: true
        tuple val(sample), path("ORFannotate/${sample}/cds.fa"), emit: cds_fa, optional: true
        tuple val(sample), path("ORFannotate/${sample}/protein.fa"), emit: protein_fa, optional: true
        path "ORFannotate/${sample}/*", emit: per_sample
    conda "${projectDir}/envs_prebuilt/ORFannotate"
    script:
    def species = params.species ?: 'human'
    def cutoff  = params.coding_cutoff ? "--coding-cutoff ${params.coding_cutoff}" : ''
    """
    set -euo pipefail
    outdir="ORFannotate/${sample}"
    mkdir -p \$outdir

    python ${params.orf_script} \\
        --gtf ${corrected_gtf} \\
        --fa  ${genome} \\
        --outdir \$outdir \\
        --species ${species} \\
        ${cutoff} \\
        > \$outdir/${sample}.orfannotate.log 2>&1
    """
}