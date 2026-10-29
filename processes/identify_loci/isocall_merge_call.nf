// Merges all per-sample profiles and calls isoforms jointly across all samples.
// Runs once after all ISOCALL_PROFILE jobs complete.
// Outputs a single GTF for all samples - equivalent to the 'bulk' collapsed GFF
// in the standard pacbio_preprocess workflow.
//
//   params.isocall_min_reads_per_isoform - min reads to report an isoform (default: 1)
//   params.isocall_use_all_chroms - include non-core contigs (default: false)
//   params.isocall_config - optional path to TOML config file

process ISOCALL_MERGE_CALL {
    tag { "all_samples" }
    publishDir "${params.outdir}/preprocessing/isocall/calls", mode: 'link'

    input:
        path profiles         // all *.profile.gz files collected into one channel item
        path ref_isoforms     // from ISOCALL_PREP_ISOFORMS
        path genome           // indexed reference FASTA

    output:
        path "isocall_merged.isoforms.gtf",        emit: isoforms_gtf
        path "isocall_merged.count_matrix.txt",    emit: count_matrix
        path "isocall_merged.closest_known.txt",   emit: closest_known

    conda "${projectDir}/envs_prebuilt/isocall_laser"

    script:
    def profile_args = profiles instanceof List
        ? profiles.collect { "--profiles ${it}" }.join(' ')
        : "--profiles ${profiles}"

    def min_reads = params.isocall_min_reads_per_isoform ?: 1
    def use_all_chroms = params.isocall_use_all_chroms ? '--use-all-chroms' : ''
    def config_arg = params.isocall_config ? "--config ${params.isocall_config}" : ''

    """
    # Merge all per-sample profiles
    ${params.isocall_bin} merge \\
        ${profile_args} \\
        --output merged.gz

    # Call isoforms from the merged profile
    # --min-reads-per-isoform 1: include all isoforms with any read support. Editable in config so up to the user

    ${params.isocall_bin} call \\
        --merged-profile merged.gz \\
        --known-isoforms ${ref_isoforms} \\
        --reference ${genome} \\
        --output-prefix isocall_merged \\
        --threads ${task.cpus} \\
        --min-reads-per-isoform ${min_reads} \\
        ${use_all_chroms} \\
        ${config_arg} \\
        --internal-priming-filter off
    
    # Decompress the GTF so downstream tools (ie SQANTI3) receive a plain file
    bgzip -d -c isocall_merged.isoforms.gtf.gz > isocall_merged.isoforms.gtf
    """
}
