// Concatenates per-sample FLAIR isoform.read.map.txt files into one bulk file.
// The bulk file is tagged 'bulk' to match channel key expected by annotate_isoforms
// Format is: isoform_id <tab> read1,read2,...  (one row per isoform)
// Note: if the same isoform appears in multiple samples, rows will be duplicated --
// build_read_isoform_lookup.py handles this by iterating all rows.
nextflow.enable.dsl = 2

process MERGE_FLAIR_MAPS {

    tag { 'all_samples' }
    publishDir "${params.outdir}/preprocessing/flair/quantify", mode: 'link'

    conda "${projectDir}/envs_prebuilt/core_utils_env"

    input:
        path(maps)   // collected list of per-sample read map files

    output:
        path("bulk.isoform.read.map.txt"), emit: bulk_map

    script:
    """
    cat ${maps} > bulk.isoform.read.map.txt
    """
}
