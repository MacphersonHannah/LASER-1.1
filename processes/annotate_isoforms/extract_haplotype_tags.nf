// extracts haplotype tags from output bam of whatshap
// haplotypes will be 1 and 2 (or unphased) but these will not necessarily match between samples or alleles
// if the user wants to match haplotypes across they'll need to do it manually afterwards

process EXTRACT_HAPLOTYPE_TAGS {
  tag { sample }
  publishDir "${params.outdir}/annotate_isoforms/alleles/${sample}", mode: 'link'

  input:
    tuple val(sample), path(hp_bam)

  output:
    tuple val(sample), path("haplotypes.${sample}.tsv")

  conda "${projectDir}/envs_prebuilt/isoform_collapse_env"

  script:
  """
  python3 ${projectDir}/scripts/extract_haplotype_tags.py \\
    --bam ${hp_bam} \\
    --out haplotypes.${sample}.tsv
  """
}
