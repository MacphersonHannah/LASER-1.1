process CLUSTER {
    tag { sample }
    cpus   { mode == 'bulk' ? 128 : 16 }
    memory { mode == 'bulk' ? '256 GB' : '64 GB' }
    // maxForks controlled by base.config
    publishDir "${params.outdir}/preprocessing/cluster/${sample}", mode: 'link'
    input:
        tuple val(sample), path(flnc_bam), val(mode)
    output:
        path("${sample}_clustered.bam"),                    emit: clustered_bam
        path("${sample}_clustered.cluster_report.csv"),    emit: cluster_report
        path("${sample}_cluster.log"),                     emit: cluster_log
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"
    script:
    """
    isoseq cluster2 -j \$(( ${task.cpus} )) --log-file ${sample}_cluster.log $flnc_bam ${sample}_clustered.bam
    # IsoSeq3 will output ${sample}_clustered.cluster_report.csv by default
    """
}
