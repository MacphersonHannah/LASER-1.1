// extract transcript summary information from a read to annotated isoform SQLite database.

process EXTRACT_TRANSCRIPT_SUMMARY {
    tag { sqlite_file.baseName.replace('.sqlite','') }
    publishDir "${params.outdir}/annotate_isoforms/transcript_summary", mode: 'link'
    conda "${projectDir}/envs_prebuilt/core_utils_env"
    input:
        path sqlite_file

    output:
        path "read_counts_summary.csv", emit: summary_csv
        path "transcript_annotations.gtf", emit: transcript_gtf
        path "transcript_counts_by_library.csv", emit: transcript_counts_csv

    script:
    """
    python3 ${projectDir}/scripts/extract_transcript_summary_from_sqlite.py \\
        --sqlite ${sqlite_file} \\
        --summary_out read_counts_summary.csv \\
        --gtf_out transcript_annotations.gtf
    """
}
