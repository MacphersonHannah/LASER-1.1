process PYCHOPPER {
    tag { sample }
    publishDir "${params.outdir}/preprocessing/pychopper/${sample}", mode: 'link'
    conda "${projectDir}/envs_prebuilt/isoform_collapse_env"

    input:
        tuple val(sample), path(fastq)

    output:
        tuple val(sample), path("${sample}.full_length.fastq.gz"),  emit: full_length
        tuple val(sample), path("${sample}.rescued.fastq.gz"),      emit: rescued

    script:
    def kit     = params.pychopper_kit
    def threads = task.cpus
    def chunk   = params.pychopper_reads_per_chunk ?: 200000
    """
    set -euo pipefail

    mkdir -p chunks out_full out_rescued out_qc

    # Split FASTQ into chunks (4 lines per read)
    zcat ${fastq} | split -l ${chunk * 4} --filter='gzip > chunks/\$FILE.fastq.gz' - chunk_

    # Run PyChopper and compress each chunk after successful completion
    ls chunks/*.fastq.gz | parallel -j ${threads} --halt soon,fail=1 \
        'pychopper -k ${kit} -m edlib -t 1 \
            -u /dev/null \
            -r out_qc/{/.}.pdf \
            -S out_qc/{/.}.tsv \
            -w out_rescued/{/.} \
            {} \
            out_full/{/.} && \
         gzip -f out_full/{/.} out_rescued/{/.}'

    # Concatenate compressed chunks
    cat out_full/*.fastq.gz > ${sample}.full_length.fastq.gz
    cat out_rescued/*.fastq.gz > ${sample}.rescued.fastq.gz

    # Validate gzip integrity
    gzip -t ${sample}.full_length.fastq.gz
    gzip -t ${sample}.rescued.fastq.gz
    """
}