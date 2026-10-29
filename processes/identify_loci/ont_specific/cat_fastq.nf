process CAT_FASTQ {
    tag { sample }

    input:
        tuple val(sample), path(fl), path(rescued)

    output:
        tuple val(sample), path("${sample}.combined.fastq.gz")

    script:
    """
    cat ${fl} ${rescued} > ${sample}.combined.fastq.gz
    """
}