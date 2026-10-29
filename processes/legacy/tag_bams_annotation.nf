// not currently used, keeps read-level annotations in the BAM by tagging reads with their transcript and gene assignments from the SQLite database
// Nice idea but files ended up unnecessarily large and it's easy to manually do this afterwards if you want to view in IGV or something

process TAG_BAMS_ANNOTATION {
  tag { sample }
  publishDir "${params.outdir}/annotate_isoforms/tagged_bams/${sample}", mode: 'copy', overwrite: true

  input:
    tuple val(sample), path(hp_bam), path(sqlite_db)

  output:
    tuple val(sample), path("${sample}.tagged.bam"), path("${sample}.tagged.bam.bai"), emit: tagged


  conda "${projectDir}/envs_prebuilt/isoform_collapse_env"
  script:
  """
  python3 ${projectDir}/scripts/tag_bam_with_annotation.py \\
    --bam ${hp_bam} \\
    --sqlite ${sqlite_db} \\
    --out_bam ${sample}.tagged.bam
  """
}
