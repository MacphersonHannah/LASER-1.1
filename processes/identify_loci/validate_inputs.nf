// Validates input CCS BAM and primer files before proceeding with analysis.
// Not used if you input demuxed bams at the start and not currently called in identify_loci (14/11/2025). Basically legacy

process VALIDATE_INPUTS {
    tag { analysis_name }
    publishDir "${params.outdir}/preprocessing/validate_inputs", mode: 'link'

    input:
      val analysis_name
      path ccs_bam
      path primers
      path validate_py

    output:
      stdout emit: msg
      path "used_config.config", optional: true, emit: cfg

    conda "${projectDir}/envs_prebuilt/core_utils_env"

    script:
    """
    python3 $validate_py $analysis_name $ccs_bam $primers

    # Also copy the config into the work dir; Nextflow will publish it
    if [ -f ${params.config_location} ]; then
      cp ${params.config_location} used_config.config
    fi
    """
}
