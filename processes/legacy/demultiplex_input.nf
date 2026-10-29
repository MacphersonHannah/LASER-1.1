// October 2025 have added an option to not use this and input demuxed files directly.

process DEMULTIPLEX {
    tag { analysis_name }
    publishDir "${params.outdir}/demultiplex", mode: 'copy', overwrite: true
    input:
        path ccs_bam
        path primers
        path parse_py
    output:
        path "*.lima.report"
        path "*.bam"
        path "*lima.log"
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"
    script:
    """
    set -euo pipefail
    echo "Parsing primers..."
    primer_info=\$(python3 $parse_py $primers)
    echo "Primer info: \$primer_info"
    echo "Running lima..."
    lima $ccs_bam $primers ${params.analysis_name}.bam -j ${task.cpus} --isoseq --peek-guess --log-level INFO --log-file ${params.analysis_name}.lima.log

    # makes sure report file is named consistently
    if [ -f ${params.analysis_name}.bam.lima.report ]; then
        mv ${params.analysis_name}.bam.lima.report ${params.analysis_name}.lima.report
    fi
    """
}

