// builds a read to annotated isoform lookup database
// Resulting SQLite database contains read-level annotations including isoform assignments, ORF annotations, SQANTI classifications, haplotype tags,
// and featureCounts read->gene assignments (separate table read_gene_fc).

// Isoquant outputs differ from isoseq3 whereby in isoseq info needed is split across two files rather than one. Therefore has platform as a requirement here
// and the called script slightly differs by platform to get those values.
// Platform mapping:
//   ont_flair -> 'flair' (FLAIR quantify read map format)
//   ont* / pacbio_isoquant / pacbio_isocall -> 'ont'  (IsoQuant model_reads format / isocall stubs)
//   pacbio -> 'pacbio' (isoseq3 cluster_report + group.txt format)

process BUILD_READ_ISOFORM_LOOKUP {
  tag { sample }
  publishDir "${params.outdir}/annotate_isoforms/lookups/", mode: 'link'
  input:
    tuple val(sample),
          path(cluster_report, stageAs: 'cluster_report_staged'),
          path(group_file,     stageAs: 'group_file_staged'),
          path(orfannotate_summary),
          path(sqanti_classification),
          path(read_sample_map),
          path(haplotype_tsvs),         
          path(featurecounts_tsvs)
    val(platform)      

  output:
    path("read_to_annotated_isoform.${sample}.sqlite"), emit: sqlite

  conda "${projectDir}/envs_prebuilt/core_utils_env"

  script:
  def hap_args = haplotype_tsvs.collect { "--haplotypes ${it}" }.join(' ')
  def fc_args  = featurecounts_tsvs.collect { "--featurecounts ${it}" }.join(' ')
  def platform_arg = (platform == 'ont_cdna' || platform == 'ont_cdna_minimap_entry' || platform == 'ont_directrna') ? 'flair'
                 : (platform.startsWith('ont') || platform.startsWith('pacbio_isoquant')) ? 'ont'
                 : 'pacbio'

  """
  
  python3 ${projectDir}/scripts/build_read_isoform_lookup.py \\
    --platform ${platform_arg} \\
    --cluster_report ${cluster_report} \\
    --group_txt ${group_file} \\
    --orf_summary ${orfannotate_summary} \\
    --sqanti_classification ${sqanti_classification} \\
    --sample ${sample} \\
    --read_sample_map ${read_sample_map} \\
    ${hap_args} \\
    ${fc_args} \\
    --out read_to_annotated_isoform.${sample}.sqlite
      
  """
}
