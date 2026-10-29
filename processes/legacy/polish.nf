// Since cluster2 produces polished reads, this process just converts clustered BAM to FASTA format for downstream use. Used to be a separate step

process POLISH {
    tag { sample }
    publishDir "${params.outdir}/polish/${sample}", mode: 'copy', overwrite: true

    input:
        tuple val(sample), path(cluster_bam)           // {sample}_clustered.bam
    output:
        path("${sample}.hq.fasta"),          emit: hq_fasta
        path("${sample}.hq.fasta.gz"),       emit: hq_fasta_gz

    conda "${projectDir}/envs_prebuilt/refine_laser"


    script:
    """
    # cluster2 already polished; just dump FASTA
    samtools fasta ${cluster_bam} > ${sample}.hq.fasta
    gzip -c    ${sample}.hq.fasta > ${sample}.hq.fasta.gz
    """
}
