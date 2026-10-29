// workflows/ont_directrna_preprocess.nf
// Direct RNA: no pychopper, different minimap2/IsoQuant flags.
// IsoQuant runs parallelised by chromosome (IsoQuant 3.13+ compatible).
// Variant calling: Clair3-RNA (in SIF container) -> FILTER_VCF (bcftools, conda).
nextflow.enable.dsl = 2

include { VALIDATE_ONT_INPUTS }  from '../processes/identify_loci/ont_specific/validate_ont_inputs.nf'
include { MINIMAP2 }         from '../processes/identify_loci/minimap2.nf'
include { PREPARE_VARIANT_BAM }   from '../processes/identify_loci/prepare_variant_bam.nf'
include { RUN_CLAIR3_RNA }       from '../processes/identify_loci/run_clair3_rna.nf'
include { FILTER_VCF }           from '../processes/identify_loci/filter_vcf.nf'
include { SPLIT_BAM_BY_CHROM }   from '../processes/identify_loci/isoquant_parallel.nf'
include { ISOQUANT_CHROM }       from '../processes/identify_loci/isoquant_parallel.nf'
include { MERGE_ISOQUANT }       from '../processes/identify_loci/isoquant_parallel.nf'
include { MAKE_GENE_BED }        from '../processes/identify_loci/make_gene_bed.nf'

workflow ont_directrna_preprocess {

  main:

    def genome_ch = Channel.value(file(params.genome))
    def gtf_ch    = Channel.value(file(params.gencode_gtf))
    def db_ch     = Channel.value(file(params.gencode_db))

    def chromosomes = Channel.of(
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
        'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
        'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
        'chr19', 'chr20', 'chr21', 'chr22', 'chrX'
    )

    // INPUT 
    def ont_input_ch = Channel
      .from(params.ont_input_files)
      .map { fastq_path ->
        def f      = file(fastq_path)
        def sample = f.getName().replaceFirst(/\.fastq(\.gz)?$/, '')
        if (!f.exists()) error "ONT fastq not found: ${fastq_path}"
        tuple(sample, f)
      }

    VALIDATE_ONT_INPUTS(ont_input_ch)

    VALIDATE_ONT_INPUTS.out.view { sample, fq ->
      "Direct RNA input validated: ${sample} -> ${fq}"
    }

    MINIMAP2(VALIDATE_ONT_INPUTS.out, Channel.value('-ax splice -k14'))

    MINIMAP2.out.bam.view { sample, bam, bai ->
      "Direct RNA mapped: ${sample} -> ${bam}"
    }

    MINIMAP2.out.bam
      .map { s, b, bai -> tuple(s, b) }
      .set { mapped_bams_ch }

    // GENE BED + DV BAM PREP 
    MAKE_GENE_BED(file(params.gencode_gtf))
    PREPARE_VARIANT_BAM(MINIMAP2.out.bam, MAKE_GENE_BED.out.bed_named)

    // PARALLELISED ISOQUANT
    SPLIT_BAM_BY_CHROM(MINIMAP2.out.bam, chromosomes)
    
    SPLIT_BAM_BY_CHROM.out
      .groupTuple(by: 0)
      .map { chrom, samples, bams, bais -> tuple(chrom, samples, bams, bais) }
      .set { isoquant_chrom_input }
    
    ISOQUANT_CHROM(isoquant_chrom_input, genome_ch, db_ch, 'forward', 'nanopore')
    
    MERGE_ISOQUANT(
      ISOQUANT_CHROM.out.gtf.map { chrom, file -> file }.collect(),
      ISOQUANT_CHROM.out.counts.map { chrom, file -> file }.collect(),
      ISOQUANT_CHROM.out.gene_counts.map { chrom, file -> file }.collect(),
      ISOQUANT_CHROM.out.model_reads.map { chrom, file -> file }.collect(),
      ISOQUANT_CHROM.out.sqanti_classification.map { chrom, file -> file }.collect(),
      ISOQUANT_CHROM.out.read_assignments.map { chrom, file -> file }.collect(),
      ISOQUANT_CHROM.out.discovered_counts.map { chrom, file -> file }.collect()
)

    MERGE_ISOQUANT.out.model_reads
      .map { tsv -> tuple('bulk', tsv) }
      .set { cluster_reports_ch }
      
    MERGE_ISOQUANT.out.model_reads
      .map { tsv -> tuple('bulk', tsv) }
      .set { collapse_groups_ch }

    MERGE_ISOQUANT.out.sqanti_classification
      .map { tsv -> tuple('bulk', tsv) }
      .set { sqanti_classification_ch }

    MINIMAP2.out.bam
      .collect()
      .map { ignored -> file("${params.outdir}/preprocessing/minimap2") }
      .set { bam_dir_ch }

    // VARIANT CALLING 
    def genome_pair_ch = Channel.value(tuple(file(params.genome), file(params.genome + '.fai')))

    PREPARE_VARIANT_BAM.out.bam
    .map { sample, bam, bai, bed -> tuple(sample, bam, bai) }
    .set { clair3_input_ch }

    RUN_CLAIR3_RNA(clair3_input_ch, genome_pair_ch, params.clair3_rna_model)

    FILTER_VCF(RUN_CLAIR3_RNA.out.vcf)

    FILTER_VCF.out.vcf
      .map { sample, vcf, tbi -> tuple(sample, vcf) }
      .set { filtered_vcf_ch }

  emit:
    mapped_bams           = mapped_bams_ch
    bam_dir               = bam_dir_ch
    dv_bams               = PREPARE_VARIANT_BAM.out.bam
    vcf                   = filtered_vcf_ch
    isoform_gffs          = MERGE_ISOQUANT.out.transcript_models_gtf.map { gtf -> tuple('bulk', gtf) }
    isoform_abundance     = MERGE_ISOQUANT.out.combined_counts
    cluster_reports       = cluster_reports_ch
    collapse_groups       = collapse_groups_ch
    sqanti_classification = sqanti_classification_ch
    gene_bed              = MAKE_GENE_BED.out.bed_named
    bam_pattern           = Channel.value('*/*.sorted.bam')
}
