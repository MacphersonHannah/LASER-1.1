
// Adds five_prime_UTR and three_prime_UTR to the ORFannotate GTF.

process AGAT_ADD_UTR {
    tag { sample }
    publishDir "${params.outdir}/annotate_isoforms/utr_annotated", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(orf_gtf)    // ORFannotate_annotated.gtf

    output:
        tuple val(sample), path("${sample}_utr_annotated.gtf"), emit: utr_gtf

    conda "${projectDir}/envs_prebuilt/agat_env"

    script:
    """
    set -euo pipefail

    # Step 1: standardise GTF -> GFF3.
    # agat_convert_sp_gxf2gxf.pl infers UTR = exon_region − CDS_region per transcript
    # and writes five_prime_UTR, three_prime_UTR, start_codon, stop_codon into the output.
    # agat_sp_complement_annotations.pl is a different tool (merges two annotation files)
    # and does NOT do this -- do not confuse the two.
    agat_convert_sp_gxf2gxf.pl \\
        -g ${orf_gtf} \\
        -o intermediate.gff3

    # Step 2: convert GFF3 back to GTF format.
    agat_convert_sp_gff2gtf.pl \\
        --gff intermediate.gff3 \\
        -o ${sample}_utr_annotated.gtf
    """
}
