// ONT cDNA preprocessing from unmapped BAMs which skips pychopper but pychopper must have been run beforehand.
// Identical to ont_cdna_preprocess except:
//   - Input: unmapped BAMs via params.ont_input_files (same param name)
//   - Validates BAM integrity instead of fastq format
//   - Converts BAM > FASTQ via samtools fastq
//   - Skips PYCHOPPER entirely (no adapter trimming, no FL/rescued split)
//   - All reads treated together (equivalent to pychopper_rescued=false)
//
// Everything from minimap2 onward is identical to ont_cdna_preprocess.
// Config is the same as ont_cdna, only params.platform differs.
// Use params.ont_input_files to point at unmapped BAMs.

nextflow.enable.dsl = 2

include { VALIDATE_ONT_BAM_INPUTS }        from '../processes/identify_loci/ont_specific/validate_ont_bam_inputs.nf'
include { BAM_TO_FASTQ }                   from '../processes/identify_loci/bam_to_fastq.nf'
include { MINIMAP2 }                       from '../processes/identify_loci/minimap2.nf'
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



workflow ont_cdna_minimap_entry_preprocess {

  main:

    def genome_ch = Channel.value(file(params.genome))
    def gtf_ch    = Channel.value(file(params.gencode_gtf))

    // RESOLVE SETTINGS
    def effective_caller = params.variant_caller ?: 'clair3'

    if (!params.ont_chemistry) error "params.ont_chemistry is not set. Please set it to 'r9' or 'r10'."
    def effective_chemistry = params.ont_chemistry.toLowerCase()

    // DeepVariant model validation (same as ont_cdna_preprocess)
    if (effective_caller == 'deepvariant') {
        def m = params.dv_model ?: 'MASSEQ'
        if (m == 'MASSEQ') {
            error "params.dv_model='${m}' is the PacBio MAS-Seq model and will produce nonsense on ONT data. Set params.dv_model='ONT_R104' in your run config (recommended for both R9 and R10 chemistries)."
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

    log.info "Mode: ont_cdna_minimap_entry (unmapped BAM input, no pychopper)"
    log.info "Variant caller: ${effective_caller}"
    if (effective_caller == 'deepvariant') {
        log.info "DeepVariant model: ${params.dv_model}"
    }
    log.info "FLAIR max_ends: ${params.flair_max_ends ?: 100}"
    log.info "FLAIR min_support: ${params.flair_min_support ?: 2}"

    // INPUT
    def ont_input_ch = Channel
        .from(params.ont_input_files)
        .map { input_path ->
            def f      = file(input_path)
            def sample = f.getName().replaceFirst(/\.(bam|fastq(\.gz)?)$/, '')
            if (!f.exists()) error "ONT input not found: ${input_path}"
            tuple(sample, f)
        }

    // VALIDATE + CONVERT
    // Validates BAM integrity, then converts to fastq for minimap2.
    VALIDATE_ONT_BAM_INPUTS(ont_input_ch)
    BAM_TO_FASTQ(VALIDATE_ONT_BAM_INPUTS.out)

    BAM_TO_FASTQ.out.fastq.view { sample, fq ->
        "ONT BAM>FASTQ complete: ${sample} > ${fq}"
    }

    //  NO PYCHOPPER
    // Pychopper has already been run upstream, so there is no FL/rescued split
    // to make here: every read is full-length as far as this workflow is
    // concerned. One alignment per sample therefore serves every downstream
    // consumer: isoform discovery (via MERGE_BAMS > FLAIR_TRANSCRIPTOME),
    // quantification (FLAIR_QUANTIFY), and variant calling + phasing
    // (PREPARE_VARIANT_BAM). ont_cdna_preprocess keeps two separate alignments
    // because there the FL and FL+rescued read sets genuinely differ.

    def all_reads = BAM_TO_FASTQ.out.fastq

    // SINGLE ALIGNMENT, shared by every downstream branch
    MINIMAP2(all_reads, Channel.value(minimap2_opts))

    MINIMAP2.out.bam.view { sample, bam, bai ->
        "ONT mapped (minimap_entry): ${sample} > ${bam}"
    }

    MINIMAP2.out.bam
        .map { s, b, bai -> tuple(s, b) }
        .set { mapped_bams_ch }

    // aligned reads > FLAIR isoform discovery
    MINIMAP2.out.bam
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

    // GENE BED + VARIANT BAM PREP
    MAKE_GENE_BED(gtf_ch)
    PREPARE_VARIANT_BAM(MINIMAP2.out.bam, MAKE_GENE_BED.out.bed_named)

    // PER-SAMPLE QUANTIFY
    // Each sample's reads are quantified against the transcriptome.
    // --generate_map produces sample.isoform.read.map.txt
    // (read_id > isoform_id) for build_read_isoform_lookup.py downstream.
    FLAIR_QUANTIFY(
        all_reads,
        FLAIR_TRANSCRIPTOME.out.isoforms_fa,
        FLAIR_TRANSCRIPTOME.out.isoforms_bed
    )

    FLAIR_QUANTIFY.out.read_map.view { sample, map ->
        "FLAIR quantify complete: ${sample} > ${map}"
    }

    // Build per-sample counts into a combined abundance channel (tagged 'bulk')
    // to match the interface expected by identify_loci/annotate_isoforms.
    // SQANTI_QC is the only consumer and only stat-checks that this file is
    // non-empty, so any one sample's counts file will do.
    FLAIR_QUANTIFY.out.counts
        .map { sample, tsv -> tsv }
        .collect()
        .map { tsvs -> tuple('bulk', tsvs[0]) }
        .set { isoform_abundance_ch }

    // Collect all sample read maps into one bulk file tagged 'bulk'.
    // annotate_isoforms filters cluster_reports and collapse_groups for s == 'bulk'
    // before joining, so sample keys would be silently dropped.
    FLAIR_QUANTIFY.out.read_map
        .map { sample, map -> map }
        .collect()
        .set { all_maps_ch }

    MERGE_FLAIR_MAPS(all_maps_ch)

    // cluster_reports: ('bulk', bulk_map) - feeds build_read_isoform_lookup.py
    // as --cluster_report with --platform flair
    MERGE_FLAIR_MAPS.out.bulk_map
        .map { f -> tuple('bulk', f) }
        .set { cluster_reports_ch }

    // collapse_groups: ('bulk', stub) - build_read_isoform_lookup.py --platform flair
    // ignores --group_txt entirely but the channel join requires it to exist.
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

    MINIMAP2.out.bam
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
