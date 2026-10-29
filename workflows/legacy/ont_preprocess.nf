// workflows/ont_preprocess.nf
// cDNA: includes pychopper. IsoQuant runs parallelised by chromosome.
// Variant calling: Clair3-RNA (default) or DeepVariant.
// Chemistry: set params.ont_chemistry = 'r10' (default) or 'r9'
nextflow.enable.dsl = 2

include { VALIDATE_ONT_INPUTS }  from '../processes/identify_loci/ont_specific/validate_ont_inputs.nf'
include { PYCHOPPER }            from '../processes/identify_loci/ont_specific/pychopper.nf'
include { MINIMAP2 }             from '../processes/identify_loci/minimap2.nf'
include { PREPARE_VARIANT_BAM }  from '../processes/identify_loci/prepare_variant_bam.nf'
include { RUN_CLAIR3_RNA }       from '../processes/identify_loci/run_clair3_rna.nf'
include { FILTER_VCF }           from '../processes/identify_loci/filter_vcf.nf'
include { RUN_DEEPVARIANT }      from '../processes/identify_loci/run_deepvariant.nf'
include { SPLIT_BAM_BY_CHROM }   from '../processes/identify_loci/isoquant_parallel.nf'
include { ISOQUANT_CHROM }       from '../processes/identify_loci/isoquant_parallel.nf'
include { MERGE_ISOQUANT }       from '../processes/identify_loci/isoquant_parallel.nf'
include { MAKE_GENE_BED }        from '../processes/identify_loci/make_gene_bed.nf'

workflow ont_preprocess {

  main:

    def genome_ch = Channel.value(file(params.genome))
    def db_ch     = Channel.value(file(params.gencode_db))

    def chromosomes = Channel.of(
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
        'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
        'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
        'chr19', 'chr20', 'chr21', 'chr22', 'chrX'
    )

    // RESOLVE SETTINGS 
    def effective_caller  = params.variant_caller ?: 'clair3'
    
    if (!params.ont_chemistry) error "params.ont_chemistry is not set. Please set it to 'r9' or 'r10' in your config."

    def effective_chemistry = params.ont_chemistry.toLowerCase()

    def minimap2_opts
    if (effective_chemistry == 'r9') {
        minimap2_opts = '-ax splice -uf'
        log.info "ONT chemistry: R9 -> minimap2 preset: splice"
    } else if (effective_chemistry == 'r10') {
        minimap2_opts = '-ax splice:hq -uf'
        log.info "ONT chemistry: R10 -> minimap2 preset: splice:hq"
    } else {
        error "Unknown ONT chemistry '${params.ont_chemistry}'. Must be 'r9' or 'r10'."
    }
        

    log.info "ONT chemistry: ${effective_chemistry} -> minimap2 opts: ${minimap2_opts}"
    log.info "Variant caller: ${effective_caller}"

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
    PYCHOPPER(VALIDATE_ONT_INPUTS.out)

    PYCHOPPER.out.full_length.view { sample, fq ->
        "ONT pychopper complete: ${sample} -> ${fq}"
    }

    def reads_for_mapping = params.pychopper_rescued
        ? PYCHOPPER.out.full_length.mix(PYCHOPPER.out.rescued)
        : PYCHOPPER.out.full_length

    MINIMAP2(reads_for_mapping, Channel.value(minimap2_opts))

    MINIMAP2.out.bam.view { sample, bam, bai ->
        "ONT mapped: ${sample} -> ${bam}"
    }

    MINIMAP2.out.bam
        .map { s, b, bai -> tuple(s, b) }
        .set { mapped_bams_ch }

    // GENE BED + VARIANT BAM PREP-
    MAKE_GENE_BED(file(params.gencode_gtf))
    PREPARE_VARIANT_BAM(MINIMAP2.out.bam, MAKE_GENE_BED.out.bed_named)

    // PARALLELISED ISOQUANT
    SPLIT_BAM_BY_CHROM(MINIMAP2.out.bam, chromosomes)

    SPLIT_BAM_BY_CHROM.out
        .groupTuple(by: 0)
        .map { chrom, samples, bams, bais -> tuple(chrom, samples, bams, bais) }
        .set { isoquant_chrom_input }

    ISOQUANT_CHROM(isoquant_chrom_input, genome_ch, db_ch, 'reverse', 'nanopore')

    MERGE_ISOQUANT(
        ISOQUANT_CHROM.out.gtf.map             { chrom, f -> f }.collect(),
        ISOQUANT_CHROM.out.counts.map          { chrom, f -> f }.collect(),
        ISOQUANT_CHROM.out.gene_counts.map     { chrom, f -> f }.collect(),
        ISOQUANT_CHROM.out.model_reads.map     { chrom, f -> f }.collect(),
        ISOQUANT_CHROM.out.sqanti_classification.map { chrom, f -> f }.collect(),
        ISOQUANT_CHROM.out.read_assignments.map { chrom, f -> f }.collect(),
        ISOQUANT_CHROM.out.discovered_counts.map { chrom, f -> f }.collect()
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

    if (effective_caller == 'deepvariant') {
        RUN_DEEPVARIANT(PREPARE_VARIANT_BAM.out.bam, genome_pair_ch)
        RUN_DEEPVARIANT.out.vcf
            .set { filtered_vcf_ch }
    } else {
        PREPARE_VARIANT_BAM.out.bam
            .map { sample, bam, bai, bed -> tuple(sample, bam, bai) }
            .set { clair3_input_ch }
        RUN_CLAIR3_RNA(clair3_input_ch, genome_pair_ch, params.clair3_rna_model)
        FILTER_VCF(RUN_CLAIR3_RNA.out.vcf)
        FILTER_VCF.out.vcf
            .map { sample, vcf, tbi -> tuple(sample, vcf) }
            .set { filtered_vcf_ch }
    }

  emit:
    mapped_bams           = mapped_bams_ch
    bam_dir               = bam_dir_ch
    dv_bams               = PREPARE_VARIANT_BAM.out.bam
    vcf                   = filtered_vcf_ch
    isoform_gffs = MERGE_ISOQUANT.out.transcript_models_gtf
            .map { gtf -> tuple('bulk', gtf) }
    isoform_abundance     = MERGE_ISOQUANT.out.combined_counts
    cluster_reports       = cluster_reports_ch
    collapse_groups       = collapse_groups_ch
    sqanti_classification = sqanti_classification_ch
    gene_bed              = MAKE_GENE_BED.out.bed_named
    bam_pattern           = Channel.value('*/*.sorted.bam')
}