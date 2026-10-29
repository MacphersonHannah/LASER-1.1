// makes downstream processing easier by creating a master file of file names (fofn) of all input bam files.

process MAKE_MASTER_FOFN {
    publishDir "${params.outdir}/preprocessing/fofn", mode: 'link'
    stageInMode 'symlink'

    input:
    path bam_files

    output:
    path "bulk.fofn"

    script:
    // Create a single space-separated list of quoted file paths
    def bam_list = bam_files.collect { it.toString() }.join(' ')

    """
    # Write absolute BAM paths into bulk.fofn
    realpath ${bam_list} | sort -u > bulk.fofn
    """
}



