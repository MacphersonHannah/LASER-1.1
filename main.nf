// main.nf
// Platform routing:
// SUPPORTED: 
//   ont_directrna       > ont_directrna_preprocess (direct RNA, Clair3-RNA, flair)
//   ont_cdna           > ont_cdna_preprocess    (cDNA, FLAIR)
//   pacbio              > pacbio_preprocess
// pacbio_fastqentry_preprocess > pacbio_fastqentry_preprocess.nf
// ont_cdna_minimap_entry > ont_cdna_minimap_entry_preprocess (cDNA from unmapped BAMs, FLAIR)

nextflow.enable.dsl = 2

include { pacbio_preprocess }                  from './workflows/pacbio_preprocess.nf'
include { ont_directrna_preprocess }           from './workflows/ont_directrna_preprocess.nf'
include { identify_loci }                      from './workflows/identify_loci.nf'
include { annotate_isoforms }                  from './workflows/annotate_isoforms.nf'
include { ont_cdna_preprocess }               from './workflows/ont_cdna_preprocess.nf'
include { pacbio_fastqentry_preprocess } from './workflows/pacbio_fastqentry_preprocess.nf'
include { ALLELIC_EFFECT_TESTING }       from './workflows/allelic_effect_testing.nf'
include { ont_cdna_minimap_entry_preprocess } from './workflows/ont_cdna_minimap_entry_preprocess.nf'

workflow {

  def platform = params.platform.toLowerCase()

  if (platform == 'pacbio') {
    pacbio_preprocess()
    identify_loci(
      pacbio_preprocess.out.mapped_bams,
      pacbio_preprocess.out.bam_dir,
      pacbio_preprocess.out.bam_pattern,
      pacbio_preprocess.out.vcf,
      pacbio_preprocess.out.isoform_gffs,
      pacbio_preprocess.out.isoform_abundance,
      pacbio_preprocess.out.cluster_reports,
      pacbio_preprocess.out.collapse_groups,
      pacbio_preprocess.out.gene_bed,
      pacbio_preprocess.out.sqanti_classification
    )


  } else if (platform == 'pacbio_fastqentry') {
    pacbio_fastqentry_preprocess()
    identify_loci(
      pacbio_fastqentry_preprocess.out.mapped_bams,
      pacbio_fastqentry_preprocess.out.bam_dir,
      pacbio_fastqentry_preprocess.out.bam_pattern,
      pacbio_fastqentry_preprocess.out.vcf,
      pacbio_fastqentry_preprocess.out.isoform_gffs,
      pacbio_fastqentry_preprocess.out.isoform_abundance,
      pacbio_fastqentry_preprocess.out.cluster_reports,
      pacbio_fastqentry_preprocess.out.collapse_groups,
      pacbio_fastqentry_preprocess.out.gene_bed,
      pacbio_fastqentry_preprocess.out.sqanti_classification
    )


  } else if (platform == 'ont_cdna') {
    ont_cdna_preprocess()
    identify_loci(
      ont_cdna_preprocess.out.mapped_bams,
      ont_cdna_preprocess.out.bam_dir,
      ont_cdna_preprocess.out.bam_pattern,
      ont_cdna_preprocess.out.vcf,
      ont_cdna_preprocess.out.isoform_gffs,
      ont_cdna_preprocess.out.isoform_abundance,
      ont_cdna_preprocess.out.cluster_reports,
      ont_cdna_preprocess.out.collapse_groups,
      ont_cdna_preprocess.out.gene_bed,
      ont_cdna_preprocess.out.sqanti_classification
    )


  } else if (platform == 'ont_cdna_minimap_entry') {
    ont_cdna_minimap_entry_preprocess()
    identify_loci(
      ont_cdna_minimap_entry_preprocess.out.mapped_bams,
      ont_cdna_minimap_entry_preprocess.out.bam_dir,
      ont_cdna_minimap_entry_preprocess.out.bam_pattern,
      ont_cdna_minimap_entry_preprocess.out.vcf,
      ont_cdna_minimap_entry_preprocess.out.isoform_gffs,
      ont_cdna_minimap_entry_preprocess.out.isoform_abundance,
      ont_cdna_minimap_entry_preprocess.out.cluster_reports,
      ont_cdna_minimap_entry_preprocess.out.collapse_groups,
      ont_cdna_minimap_entry_preprocess.out.gene_bed,
      ont_cdna_minimap_entry_preprocess.out.sqanti_classification
    )

  } else if (platform == 'ont_directrna') {
    ont_directrna_preprocess()
    identify_loci(
      ont_directrna_preprocess.out.mapped_bams,
      ont_directrna_preprocess.out.bam_dir,
      ont_directrna_preprocess.out.bam_pattern,
      ont_directrna_preprocess.out.vcf,
      ont_directrna_preprocess.out.isoform_gffs,
      ont_directrna_preprocess.out.isoform_abundance,
      ont_directrna_preprocess.out.cluster_reports,
      ont_directrna_preprocess.out.collapse_groups,
      ont_directrna_preprocess.out.gene_bed,
      ont_directrna_preprocess.out.sqanti_classification
    )

  } else {
    error "Unknown platform '${params.platform}'. Expected: pacbio | pacbio_fastqentry | ont_cdna | ont_cdna_minimap_entry | ont_directrna"
  }

  annotate_isoforms(
    identify_loci.out.isoform_gffs,
    identify_loci.out.isoform_abundance,
    identify_loci.out.cluster_reports,
    identify_loci.out.collapse_groups,
    identify_loci.out.haplotagged_bams,
    identify_loci.out.read_sample_map,
    identify_loci.out.sqanti_classification
  )

  ALLELIC_EFFECT_TESTING(
    annotate_isoforms.out.read_counts_csv,
    identify_loci.out.gene_bed,
    annotate_isoforms.out.orf_summary
  )
}
