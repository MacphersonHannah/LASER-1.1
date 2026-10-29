process SQANTI_QC {
    tag { sample }
    publishDir "${params.outdir}/annotate_isoforms/SQANTI_QC", mode: 'link'
    input:
        tuple val(sample), path(isoforms), path(fl_count)
        path gtf
        path genome
        path cage
        path poly_peak
        path poly_motifs
    output:
        tuple val(sample), path("Sqanti/${sample}/${sample}_corrected.gtf"), emit: corrected
        tuple val(sample), path("Sqanti/${sample}/${sample}_classification.txt"), emit: classification
        path "Sqanti/${sample}/*_junctions.txt"

    conda "${projectDir}/envs_prebuilt/SQANTI3"

    script:
    def coverage_flag = (params.species?.toLowerCase() == 'human') ? "--coverage ${params.intropolis}" : ""
    """
    set -euo pipefail

    res_dir="Sqanti/${sample}"
    prefix="${sample}"
    mkdir -p \$res_dir

    if [ ! -s "${isoforms}" ] || [ ! -s "${fl_count}" ]; then
        echo "Error: Missing or empty input files for ${sample}." >&2
        exit 1
    fi

    export PYTHONPATH="${params.cupcake_path}/sequence:\${PYTHONPATH:-}"

    python ${params.sqanti_dir}/sqanti3_qc.py \\
        "${isoforms}" \\
        "${gtf}" \\
        "${genome}" \\
        --CAGE_peak "${cage}" \\
        --polyA_peak "${poly_peak}" \\
        --polyA_motif_list "${poly_motifs}" \\
        ${coverage_flag} \\
        -t ${task.cpus} -n ${params.sqanti_chunks} \\
        -o "\$prefix" -d "\$res_dir"
    """
}
