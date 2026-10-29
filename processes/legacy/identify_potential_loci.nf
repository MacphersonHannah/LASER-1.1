process IDENTIFY_POTENTIAL_LOCI {
  tag { sample }
  cache false
  publishDir "${params.outdir}/identify_loci/${sample}", mode: 'copy', overwrite: true

  input:
    tuple val(sample), path(dv_vcf), path(phased_vcf), path(hp_bam), path(bed_file), path(sqlite_db)

  output:
    path "${sample}.loci.csv", emit: loci
    path "${sample}.variants_by_locus.tsv", emit: variants

  conda "${projectDir}/envs_prebuilt/isoform_collapse_env"

  script:
  """
  set -euo pipefail
  python3 ${projectDir}/scripts/potential_ASE_loci.py \
    --dv_vcf ${dv_vcf} \
    --phased_vcf ${phased_vcf} \
    --bam ${hp_bam} \
    --bed ${bed_file} \
    --sqlite ${sqlite_db} \
    --out ${sample}.loci.csv \
    --only_pass \
    --variants_out ${sample}.variants_by_locus.tsv \
    ${params.min_variant_sites ? "--min_variant_sites ${params.min_variant_sites}" : ""} \
    ${params.min_total_reads ? "--min_total_reads ${params.min_total_reads}" : ""} \
    ${params.min_minor_count ? "--min_minor_count ${params.min_minor_count}" : ""} \
    ${params.min_minor_frac ? "--min_minor_frac ${params.min_minor_frac}" : ""} \
    ${params.min_reads_per_haplotype ? "--min_reads_per_haplotype ${params.min_reads_per_haplotype}" : ""} \
    ${params.maximum_unphased_frac ? "--max_unphased_frac ${params.maximum_unphased_frac}" : ""}
  """
}
