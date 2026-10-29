// Creates an empty stub file per sample to satisfy collapse_groups channel
// join in annotate_isoforms. build_read_isoform_lookup.py --platform flair
// ignores the group_txt argument entirely, so this file can be empty.
nextflow.enable.dsl = 2

process CREATE_STUB {

    tag { sample }

    conda "${projectDir}/envs_prebuilt/core_utils_env"

    input:
        tuple val(sample), path(any_file)   // piggybacks on read_map channel

    output:
        tuple val(sample), path("${sample}.group_stub.txt"), emit: stub

    script:
    """
    touch ${sample}.group_stub.txt
    """
}
