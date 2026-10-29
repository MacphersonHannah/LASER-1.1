// workflows/ont_stringtie_preprocess.nf
// ONT cDNA preprocessing using StringTie2 for transcript discovery/quantification.
// Replaces IsoQuant with StringTie2 - much faster, suitable for whole transcriptome.
// Variant calling: Clair3-RNA (default) or DeepVariant - identical to ont_preprocess.nf.
// Chemistry: set params.ont_chemistry = 'r10' or 'r9'
//
// To revert to IsoQuant: change platform = 'ont' in config.
// This file is self-contained and does not modify ont_preprocess.nf.
nextflow.enable.dsl = 2

include { VALIDATE_ONT_INPUTS }  from '../processes/identify_loci/ont_specific/validate_ont_inputs.nf'
include { PYCHOPPER }            from '../processes/identify_loci/ont_specific/pychopper.nf'
include { MINIMAP2 }             from '../processes/identify_loci/minimap2.nf'
include { PREPARE_VARIANT_BAM }  from '../processes/identify_loci/prepare_variant_bam.nf'
include { RUN_CLAIR3_RNA }       from '../processes/identify_loci/run_clair3_rna.nf'
include { FILTER_VCF }           from '../processes/identify_loci/filter_vcf.nf'
include { RUN_DEEPVARIANT }      from '../processes/identify_loci/run_deepvariant.nf'
include { MAKE_GENE_BED }        from '../processes/identify_loci/make_gene_bed.nf'
include { STRINGTIE }            from '../processes/identify_loci/stringtie.nf'
include { STRINGTIE_MERGE }      from '../processes/identify_loci/stringtie.nf'
include { STRINGTIE_QUANT }      from '../processes/identify_loci/stringtie.nf'
include { MERGE_STRINGTIE_READS } from '../processes/identify_loci/stringtie.nf'
include { STRINGTIE_COUNT_MATRIX } from '../processes/identify_loci/stringtie.nf'

workflow ont_stringtie_preprocess {

  main:

    def genome_ch = Channel.value(file(params.genome))

    // RESOLVE SETTINGS
    def effective_caller = params.variant_caller ?: 'clair3'

    if (!params.ont_chemistry) error "params.ont_chemistry is not set. Please set it to 'r9' or 'r10' in your config."
    def effective_chemistry = params.ont_chemistry.toLowerCase()

    def minimap2_opts
    def strand_flag     // StringTie2 strand flag
    def fc_strand       // featureCounts strand: 2=reverse (PCB114 cDNA), 1=forward (direct RNA)

    if (effective_chemistry == 'r9') {
        minimap2_opts = '-ax splice -uf'
        strand_flag   = '--rf'
        fc_strand     = 2
        log.info "ONT chemistry: R9 -> minimap2: splice, StringTie2: --rf"
    } else if (effective_chemistry == 'r10') {
        minimap2_opts = '-ax splice:hq -uf'
        strand_flag   = '--rf'
        fc_strand     = 2
        log.info "ONT chemistry: R10 -> minimap2: splice:hq, StringTie2: --rf"
    } else {
        error "Unknown ONT chemistry '${params.ont_chemistry}'. Must be 'r9' or 'r10'."
    }

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

    // GENE BED + VARIANT BAM PREP
    MAKE_GENE_BED(file(params.gencode_gtf))
    PREPARE_VARIANT_BAM(MINIMAP2.out.bam, MAKE_GENE_BED.out.bed_named)

    // STRINGTIE TRANSCRIPT DISCOVERY
    def gtf_ch = Channel.value(file(params.gencode_gtf))

    // Per-sample discovery
    STRINGTIE(
        mapped_bams_ch,
        gtf_ch,
        Channel.value(strand_flag)
    )

    // Merge into unified catalogue
    STRINGTIE_MERGE(
        STRINGTIE.out.gtf.collect(),
        gtf_ch
    )

    // Per-sample quantification against merged catalogue
    STRINGTIE_QUANT(
        mapped_bams_ch,
        STRINGTIE_MERGE.out.merged_gtf,
        Channel.value(strand_flag),
        Channel.value(fc_strand)
    )

    // Merge read -> transcript assignments into bulk file
    MERGE_STRINGTIE_READS(
        STRINGTIE_QUANT.out.transcript_reads
            .map { sample, f -> f }
            .collect()
    )

    // Build count matrix from per-sample quant GTFs
    STRINGTIE_COUNT_MATRIX(
        STRINGTIE_QUANT.out.quant_gtf
            .map { sample, f -> f }
            .collect()
    )

    // Channel wiring for annotate_isoforms
    // isoform_gffs: (sample, gtf) use merged GTF tagged as 'bulk'
    STRINGTIE_MERGE.out.merged_gtf
        .map { gtf -> tuple('bulk', gtf) }
        .set { isoform_gffs_ch }

    // cluster_reports and collapse_groups both use the bulk read -> transcript file
    MERGE_STRINGTIE_READS.out.model_reads
        .map { tsv -> tuple('bulk', tsv) }
        .set { cluster_reports_ch }

    MERGE_STRINGTIE_READS.out.model_reads
        .map { tsv -> tuple('bulk', tsv) }
        .set { collapse_groups_ch }

    MINIMAP2.out.bam
        .collect()
        .map { ignored -> file("${params.outdir}/preprocessing/minimap2") }
        .set { bam_dir_ch }

    // VARIANT CALLING--
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
    isoform_gffs          = isoform_gffs_ch
    isoform_abundance     = STRINGTIE_COUNT_MATRIX.out.counts
    cluster_reports       = cluster_reports_ch
    collapse_groups       = collapse_groups_ch
    sqanti_classification = Channel.empty()
    gene_bed              = MAKE_GENE_BED.out.bed_named
    bam_pattern           = Channel.value('*/*.sorted.bam')
}
