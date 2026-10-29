// RUN_DEEPVARIANT removed from this workflow - variant calling now happens in each preprocess so we can offer Clair3 too
// vcf is taken as input and passed through to WHATSHAP_PHASE_TAG.

nextflow.enable.dsl = 2

include { MAKE_READ_SAMPLE_MAP }  from '../processes/identify_loci/make_read_sample_map.nf'
include { WHATSHAP_PHASE_TAG }    from '../processes/identify_loci/whatshap_phase_tag.nf'

workflow identify_loci {

  take:
    mapped_bams           // (sample, bam) - primary aligned reads per barcode
    bam_dir               // path - directory containing mapped BAMs
    bam_pattern           // string - glob for MAKE_READ_SAMPLE_MAP
    vcf                   // (sample, vcf) - variant calls from preprocess
    isoform_gffs          // path - bulk transcript models GTF
    isoform_abundance     // path - bulk transcript counts TSV
    cluster_reports       // (sample, tsv)
    collapse_groups       // (sample, tsv)
    gene_bed              // path
    sqanti_classification // (sample, tsv)

  main:

    def genome_pair_ch = Channel.value(tuple(file(params.genome), file(params.genome + '.fai')))


    // READ-SAMPLE MAP 
    MAKE_READ_SAMPLE_MAP(bam_dir, bam_pattern)

    // WHATSHAP PHASE + HAPLOTAG 
    // vcf comes from preprocess (DV, Clair3-RNA, or other caller).
    vcf
      .join(mapped_bams)
      .combine(genome_pair_ch)
      .map { s, vcf_file, bam, ref_fa, ref_fai -> tuple(s, vcf_file, bam, ref_fa, ref_fai) }
      .set { whatshap_inputs }

    WHATSHAP_PHASE_TAG(whatshap_inputs)

  emit:
    isoform_gffs          = isoform_gffs
    isoform_abundance     = isoform_abundance
    cluster_reports       = cluster_reports
    collapse_groups       = collapse_groups
    sqanti_classification = sqanti_classification
    haplotagged_bams      = WHATSHAP_PHASE_TAG.out.hp_bam
    gene_bed              = gene_bed
    read_sample_map       = MAKE_READ_SAMPLE_MAP.out
    dv_vcfs               = vcf                  // pass-through (name kept for annotate_isoforms compat)
    phased_vcf_output     = WHATSHAP_PHASE_TAG.out.phased_vcf
}
