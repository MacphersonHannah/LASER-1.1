// ONT cDNA preprocessing with FLAIR. Design to keep it as similar to pacbio as possible: 
// FL-only reads are used for isoform discovery (FLAIR_TRANSCRIPTOME), FL+rescued combined for
// variant calling and quantification (FLAIR_QUANTIFY).
// If params.pychopper_rescued = false, use FL only.
// See ont_directrna_preprocess.nf for direct RNA adaptation.

nextflow.enable.dsl = 2

include { VALIDATE_ONT_INPUTS }            from '../processes/identify_loci/ont_specific/validate_ont_inputs.nf'
include { PYCHOPPER }                      from '../processes/identify_loci/ont_specific/pychopper.nf'
include { CAT_FASTQ }                      from '../processes/identify_loci/ont_specific/cat_fastq.nf'
include { MINIMAP2 as MINIMAP2_FL }        from '../processes/identify_loci/minimap2.nf'
include { MINIMAP2 as MINIMAP2_ALL }       from '../processes/identify_loci/minimap2.nf'
include { MERGE_BAMS }                     from '../processes/identify_loci/merge_bams.nf'
include { FLAIR_TRANSCRIPTOME }            from '../processes/identify_loci/flair_transcriptome.nf'
include { FLAIR_QUANTIFY }                 from '../processes/identify_loci/flair_quantify.nf'
include { PREPARE_VARIANT_BAM }            from '../processes/identify_loci/prepare_variant_bam.nf'
include { RUN_CLAIR3_RNA }                 from '../processes/identify_loci/run_clair3_rna.nf'
include { FILTER_VCF }                     from '../processes/identify_loci/filter_vcf.nf'
include { RUN_DEEPVARIANT }                from '../processes/identify_loci/run_deepvariant.nf'
include { MAKE_GENE_BED }                  from '../processes/identify_loci/make_gene_bed.nf'
include { CREATE_STUB }                    from '../processes/identify_loci/create_stub.nf'
include { MERGE_FLAIR_MAPS }               from '../processes/identify_loci/merge_flair_maps.nf'



workflow ont_cdna_preprocess {

  main:

    def genome_ch = Channel.value(file(params.genome))
    def gtf_ch    = Channel.value(file(params.gencode_gtf))
    def effective_caller = params.variant_caller ?: 'clair3'

    if (!params.ont_chemistry) error "params.ont_chemistry is not set. Please set it - most likely to 'r9' or 'r10'."
    def effective_chemistry = params.ont_chemistry.toLowerCase()

    // The RUN_DEEPVARIANT process reads params.dv_model directly (default 'MASSEQ' for PacBio). For ONT data ONT_R104 is trained on R10 and works reasonably on R9 too.
    // Clair3 is preferable for ONT data hence setting it as the default above but leaving it up to the user
    if (effective_caller == 'deepvariant') {
        def m = params.dv_model ?: 'MASSEQ'
        if (m == 'MASSEQ') {
            error "params.dv_model='${m}' is the PacBio MAS-Seq model and is not recommended for ONT data. Set params.dv_model='ONT_R104' in your run config (recommended for both R9 and R10 chemistries)."
        }
    }

    def minimap2_opts
    if (effective_chemistry == 'r9') {
        minimap2_opts = '-ax splice -uf -k 14'
        log.info "ONT chemistry: R9 > minimap2 preset: splice -k 14"
    } else if (effective_chemistry == 'r10') {
        minimap2_opts = '-ax splice:hq -uf -k 14'
        log.info "ONT chemistry: R10 > minimap2 preset: splice:hq -k 14"
    } else {
        error "Unknown ONT chemistry '${params.ont_chemistry}'. Must be 'r9' or 'r10'."
    }

    log.info "Variant caller: ${effective_caller}"
    if (effective_caller == 'deepvariant') {
        log.info "DeepVariant model: ${params.dv_model}"
    }
    log.info "FLAIR max_ends: ${params.flair_max_ends ?: 100}"
    log.info "FLAIR min_support: ${params.flair_min_support ?: 2}"

    def use_rescued = params.pychopper_rescued ?: false
    log.info "Pychopper rescued reads: ${use_rescued ? 'included (combined stream)' : 'excluded'}"

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
        "ONT pychopper complete: ${sample} > ${fq}"
    }

    // Split design (pychopper output assuming you have rescued reads enabled)
    //
    // 1. FL only > isoform discovery (FLAIR_TRANSCRIPTOME)
    // Full-length reads have reliable 5'/3' ends > clean transcript models.
    // Similar to PacBio HQ (clustered consensus) reads.
    //
    // 2. FL + rescued reads combined > everything else
    // Maximum read depth for variant calling, WhatsHap phasing, and FLAIR_QUANTIFY. Rescued reads have good sequence quality, just
    // more unreliable boundaries. Similar to PacBio FLNC reads.
    //
    // When pychopper_rescued = false, both streams are FL only (no CAT_FASTQ).

    def fl_reads = PYCHOPPER.out.full_length // (sample, fastq) x 4

    def combined_reads
    if (use_rescued) {
        PYCHOPPER.out.full_length
            .join(PYCHOPPER.out.rescued)
            .set { joined_ch }
        CAT_FASTQ(joined_ch)
        combined_reads = CAT_FASTQ.out // (sample, fastq) x 4
    } else {
        combined_reads = fl_reads
    }

    // FL only > FLAIR isoform discovery
    MINIMAP2_FL(fl_reads, Channel.value(minimap2_opts))

    MINIMAP2_FL.out.bam
        .map { sample, bam, bai -> bam }
        .collect()
        .set { fl_bams_ch }

    MERGE_BAMS(fl_bams_ch)

    FLAIR_TRANSCRIPTOME(
        MERGE_BAMS.out.bam,
        MERGE_BAMS.out.bai,
        genome_ch,
        gtf_ch
    )

    // combined > mapping, variant calling, quantification 
    MINIMAP2_ALL(combined_reads, Channel.value(minimap2_opts))

    MINIMAP2_ALL.out.bam.view { sample, bam, bai ->
        "ONT mapped (combined): ${sample} > ${bam}"
    }

    MINIMAP2_ALL.out.bam
        .map { s, b, bai -> tuple(s, b) }
        .set { mapped_bams_ch }

    // GENE BED + VARIANT BAM PREP 
    MAKE_GENE_BED(gtf_ch)
    PREPARE_VARIANT_BAM(MINIMAP2_ALL.out.bam, MAKE_GENE_BED.out.bed_named)

    // PER-SAMPLE QUANTIFY
    FLAIR_QUANTIFY(
        combined_reads,
        FLAIR_TRANSCRIPTOME.out.isoforms_fa,
        FLAIR_TRANSCRIPTOME.out.isoforms_bed
    )

    FLAIR_QUANTIFY.out.read_map.view { sample, map ->
        "FLAIR quantify complete: ${sample} > ${map}"
    }

    // Build per-sample counts into a combined abundance channel (tagged 'bulk')
    FLAIR_QUANTIFY.out.counts
        .map { sample, tsv -> tsv }
        .collect()
        .map { tsvs -> tuple('bulk', tsvs[0]) }
        .set { isoform_abundance_ch }
        
    // Collect all per-sample read maps into one bulk file tagged 'bulk'
    FLAIR_QUANTIFY.out.read_map
        .map { sample, map -> map }
        .collect()
        .set { all_maps_ch }

    MERGE_FLAIR_MAPS(all_maps_ch)

    // cluster_reports: ('bulk', bulk_map) feeds build_read_isoform_lookup.py
    // as cluster_report with platform flair
    MERGE_FLAIR_MAPS.out.bulk_map
        .map { f -> tuple('bulk', f) }
        .set { cluster_reports_ch }

    // collapse_groups: ('bulk', stub) build_read_isoform_lookup.py platform flair
    // ignores group_txt entirely but the channel join requires it to exist.
    CREATE_STUB(
        MERGE_FLAIR_MAPS.out.bulk_map
            .map { f -> tuple('bulk', f) }
    )
    CREATE_STUB.out.stub
        .set { collapse_groups_ch }

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

    MINIMAP2_ALL.out.bam
        .collect()
        .map { ignored -> file("${params.outdir}/preprocessing/minimap2") }
        .set { bam_dir_ch }

  emit:
    mapped_bams           = mapped_bams_ch
    bam_dir               = bam_dir_ch
    dv_bams               = PREPARE_VARIANT_BAM.out.bam
    vcf                   = filtered_vcf_ch
    gene_bed              = MAKE_GENE_BED.out.bed_named
    isoform_gffs          = FLAIR_TRANSCRIPTOME.out.isoforms_gtf
                                .map { gtf -> tuple('bulk', gtf) }
    isoform_abundance     = isoform_abundance_ch
    cluster_reports       = cluster_reports_ch
    collapse_groups       = collapse_groups_ch
    sqanti_classification = Channel.empty()
    bam_pattern           = Channel.value('*/*.sorted.bam')
}
