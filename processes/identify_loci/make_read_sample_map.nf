process MAKE_READ_SAMPLE_MAP {
  publishDir "${params.outdir}/identify_loci/read_sample_map", mode: 'link'
  conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

  input:
    path(bam_dir)
    val(bam_pattern)   // eg '*/*.mapped.bam' for PacBio, '*/*.sorted.bam' for ONT

  output:
    path("read_sample_map.tsv")

  script:
  """
  bash ${projectDir}/scripts/make_read_sample_map.sh ${bam_dir} read_sample_map.tsv ${bam_pattern}
  """
}