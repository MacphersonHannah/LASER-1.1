// Standard isoseq refine. Takes in BAM of full-length non-chimeric reads and primer file to produce refined FLNC BAM.
// Refine is now optional as input files can be refined already

process REFINE {
    tag { sample }
    publishDir "${params.outdir}/preprocessing/refine/${sample}", mode: 'link'
    input:
        tuple val(sample), path(input_bam), val(mode)
        path primers
    output:
        path "${sample}_flnc.bam"
        path "${sample}_flnc.report.csv"
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"
    script:
    """
    isoseq refine --require-polya -j \$(( ${task.cpus} )) $input_bam $primers ${sample}_flnc.bam > ${sample}_flnc.report.csv
    """
}
