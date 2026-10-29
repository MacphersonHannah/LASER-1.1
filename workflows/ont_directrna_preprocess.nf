// ONT direct RNA preprocessing with FLAIR. Forked from ont_cdna_preprocess.nf:
//   - No pychopper (no cDNA adapters)
//   - minimap2 splice preset + --splice-flank=no
//   - FLAIR_TRANSCRIPTOME: no --stringent, --check_splice, --trust_ends
//   - FLAIR_QUANTIFY: no --stringent, --check_splice (keeps --trust_ends for 3')
//   - flair_min_support default 5 to reduce degradation artefacts but up to the user

nextflow.enable.dsl = 2

include { VALIDATE_ONT_INPUTS }  from '../processes/identify_loci/ont_specific/validate_ont_inputs.nf'
include { MINIMAP2 }             from '../processes/identify_loci/minimap2.nf'
include { MERGE_BAMS }           from '../processes/identify_loci/merge_bams.nf'
include { FLAIR_TRANSCRIPTOME_DRNA } from '../processes/identify_loci/flair_transcriptome_drna.nf'
include { FLAIR_QUANTIFY_DRNA }  from '../processes/identify_loci/flair_quantify_drna.nf'
include { CREATE_STUB }          from '../processes/identify_loci/create_stub.nf'
include { MERGE_FLAIR_MAPS }     from '../processes/identify_loci/merge_flair_maps.nf'
include { PREPARE_VARIANT_BAM }  from '../processes/identify_loci/prepare_variant_bam.nf'
include { RUN_CLAIR3_RNA }       from '../processes/identify_loci/run_clair3_rna.nf'
include { FILTER_VCF }           from '../processes/identify_loci/filter_vcf.nf'
include { RUN_DEEPVARIANT }      from '../processes/identify_loci/run_deepvariant.nf'
include { MAKE_GENE_BED }        from '../processes/identify_loci/make_gene_bed.nf'

workflow ont_directrna_preprocess {

  main:

    def genome_ch = Channel.value(file(params.genome))
    def gtf_ch    = Channel.value(file(params.gencode_gtf))
    def effective_caller = params.variant_caller ?: 'clair3'

    if (!params.ont_chemistry) error "params.ont_chemistry is not set. Please set it - most likely to 'r9' or 'r10'."
    def effective_chemistry = params.ont_chemistry.toLowerCase()

    // Direct RNA uses standard splice preset regardless of chemistry
    // splice:hq is designed for high-accuracy cDNA reads only
    def minimap2_opts
    if (effective_chemistry == 'r9') {
        minimap2_opts = '-ax splice -uf -k 14splice-flank=no'
        log.info "ONT direct RNA chemistry: R9 > minimap2 preset: splice"
    } else if (effective_chemistry == 'r10') {
        minimap2_opts = '-ax splice -uf -k 14splice-flank=no'
        log.info "ONT direct RNA chemistry: R10 > minimap2 preset: splice"
    } else {
        error "Unknown ONT chemistry '${params.ont_chemistry}'. Must be 'r9' or 'r10'."
    }

    log.info "Mode: direct RNA (no pychopper, relaxed FLAIR flags)"
    log.info "Variant caller: ${effective_caller}"
    log.info "FLAIR max_ends: ${params.flair_max_ends ?: 10}"
    log.info "FLAIR min_support: ${params.flair_min_support ?: 5}"

    // INPUT 
    // No PYCHOPPER for direct RNA so reads go straight to mapping
    def ont_input_ch = Channel
        .from(params.ont_input_files)
        .map { fastq_path ->
            def f      = file(fastq_path)
            def sample = f.getName().replaceFirst(/\.fastq(\.gz)?$/, '')
            if (!f.exists()) error "ONT fastq not found: ${fastq_path}"
            tuple(sample, f)
        }

    VALIDATE_ONT_INPUTS(ont_input_ch)

    def reads_for_mapping = VALIDATE_ONT_INPUTS.out

    // GENOME ALIGNMENT
    MINIMAP2(reads_for_mapping, Channel.value(minimap2_opts))

    MINIMAP2.out.bam.view { sample, bam, bai ->
        "ONT direct RNA mapped: ${sample} > ${bam}"
    }

    MINIMAP2.out.bam
        .map { s, b, bai -> tuple(s, b) }
        .set { mapped_bams_ch }

    // GENE BED + VARIANT BAM PREP 
    MAKE_GENE_BED(gtf_ch)
    PREPARE_VARIANT_BAM(MINIMAP2.out.bam, MAKE_GENE_BED.out.bed_named)

    // MERGE ALL BAMS > FLAIR TRANSCRIPTOME
    MINIMAP2.out.bam
        .map { sample, bam, bai -> bam }
        .collect()
        .set { all_bams_ch }

    MERGE_BAMS(all_bams_ch)

    FLAIR_TRANSCRIPTOME_DRNA(
        MERGE_BAMS.out.bam,
        MERGE_BAMS.out.bai,
        genome_ch,
        gtf_ch
    )

    // PER-SAMPLE QUANTIFY
    FLAIR_QUANTIFY_DRNA(
        reads_for_mapping,
        FLAIR_TRANSCRIPTOME_DRNA.out.isoforms_fa,
        FLAIR_TRANSCRIPTOME_DRNA.out.isoforms_bed
    )

    FLAIR_QUANTIFY_DRNA.out.read_map.view { sample, map ->
        "FLAIR quantify complete: ${sample} > ${map}"
    }

    FLAIR_QUANTIFY_DRNA.out.counts
        .map { it[1] }
        .collect()
        .map { tsvs -> tuple('bulk', tsvs[0]) }   // tag as ('bulk', path) for annotate_isoforms.nf
        .set { isoform_abundance_ch }            

    FLAIR_QUANTIFY_DRNA.out.read_map
        .map { sample, map -> map }
        .collect()
        .set { all_maps_ch }

    MERGE_FLAIR_MAPS(all_maps_ch)

    MERGE_FLAIR_MAPS.out.bulk_map
        .map { f -> tuple('bulk', f) }
        .set { cluster_reports_ch }

    CREATE_STUB(
        MERGE_FLAIR_MAPS.out.bulk_map
            .map { f -> tuple('bulk', f) }
    )
    CREATE_STUB.out.stub
        .set { collapse_groups_ch }

    // VARIANT CALLING
    def genome_pair_ch = Channel.value(tuple(file(params.genome), file(params.genome + '.fai')))

    // clair3 prefered for ONT but leaving it to user's descretion
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
    isoform_gffs          = FLAIR_TRANSCRIPTOME_DRNA.out.isoforms_gtf
                                .map { gtf -> tuple('bulk', gtf) }
    isoform_abundance     = isoform_abundance_ch
    cluster_reports       = cluster_reports_ch
    collapse_groups       = collapse_groups_ch
    sqanti_classification = Channel.empty()
    bam_pattern           = Channel.value('*/*.sorted.bam')
}
