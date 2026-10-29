process FLNC_FASTQ_TO_BAM {
    tag { sample }
    publishDir "${params.outdir}/preprocessing/flnc_input/${sample}", mode: 'link'

    input:
        tuple val(sample), path(flnc_file)

    output:
        tuple val(sample), path("${sample}_flnc.bam")

    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    script:
    """
    set -euo pipefail
    case "${flnc_file}" in
        *.bam)
            cp "${flnc_file}" "${sample}_flnc.bam"
            ;;
        *)
            cp ${params.pacbio_template_bam} ./local_template.bam
            zcat "${flnc_file}" > tmp_reads.fastq
            python3 ${projectDir}/scripts/convert_to_pacbio_bam.py \
                tmp_reads.fastq \
                ./local_template.bam \
                ${sample}_converted
            mv ${sample}_converted.converted.bam ${sample}_flnc.bam
            rm tmp_reads.fastq local_template.bam
            ;;
    esac
    """
}