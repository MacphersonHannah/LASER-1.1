// workflows/pacbio_isocall_preprocess.nf
// PacBio preprocessing using isocall for joint isoform calling.
// Replaces CLUSTER + PBMM2_HQ + COLLAPSE_ISOFORMS + GET_ABUNDANCE with:
// ISOCALL_PREP_ISOFORMS -> ISOCALL_PROFILE (per sample) -> ISOCALL_MERGE_CALL
// Everything downstream (WhatsHap, DeepVariant, SQANTI3, annotation) is unchanged.
// Set params.platform = 'pacbio_isocall' and params.isocall_bin = '/path/to/isocall'
// Variant calling: DeepVariant by default, Clair3-RNA if params.variant_caller = 'clair3'

nextflow.enable.dsl = 2
import java.nio.file.Path
import groovy.json.JsonSlurper

include { MAKE_MASTER_FOFN }       from '../processes/identify_loci/make_master_fofn.nf'
include { VALIDATE_INPUTS }        from '../processes/identify_loci/validate_inputs.nf'
include { REFINE }                 from '../processes/identify_loci/refine.nf'
include { PBMM2_FLNC }            from '../processes/identify_loci/pbmm2_flnc.nf'
include { ISOCALL_PREP_ISOFORMS }  from '../processes/identify_loci/isocall_prep_isoforms.nf'
include { ISOCALL_PROFILE }        from '../processes/identify_loci/isocall_profile.nf'
include { ISOCALL_MERGE_CALL }     from '../processes/identify_loci/isocall_merge_call.nf'
include { CREATE_STUB }            from '../processes/identify_loci/create_stub.nf'
include { PREPARE_DV_BAM }         from '../processes/identify_loci/prepare_dv.nf'
include { MAKE_GENE_BED }          from '../processes/identify_loci/make_gene_bed.nf'
include { RUN_DEEPVARIANT }        from '../processes/identify_loci/run_deepvariant.nf'
include { RUN_CLAIR3_RNA }         from '../processes/identify_loci/run_clair3_rna.nf'
include { FILTER_VCF }             from '../processes/identify_loci/filter_vcf.nf'
include { FEATURECOUNTS_ISOCALL_READ_MAP } from '../processes/identify_loci/featurecounts_isocall.nf'

workflow pacbio_isocall_preprocess {

  main:

    log.info "Project dir: ${projectDir}"
    log.info "Scanning library directories... [PacBio + isocall mode]"

    if (!params.isocall_bin) {
        error "params.isocall_bin is not set. Please specify the full path to the isocall binary."
    }

    // RESOLVE VARIANT CALLER
    def effective_caller = params.variant_caller ?: 'deepvariant'
    def effective_dv_model = params.dv_model ?: 'MASSEQ'
    def effective_clair3_model = params.clair3_model ?: 'hifi_sequel2_minimap2'
    log.info "Variant caller: ${effective_caller}"

    def library_dirs = (params.library_dirs instanceof List) ? params.library_dirs : [params.library_dirs]
    def all_libraries = []

    library_dirs.each { lib_dir ->
        def lib_path = file(lib_dir)
        if (!lib_path.exists()) error "Library directory not found: ${lib_path}"

        def primers = lib_path.listFiles().find { it.name ==~ /.*\.(fa|fasta)$/ }
        if (!primers) error "No primer (.fa or .fasta) file found in ${lib_path}"

        def bam_files = lib_path.listFiles().findAll { it.name.endsWith('.bam') }
        if (bam_files.isEmpty()) error "No BAM files found in ${lib_path}"

        all_libraries << [lib_path.getName(), primers, bam_files]
    }

    log.info "Found ${all_libraries.size()} library(ies)."

    def genome_ch = Channel.value(file(params.genome))

    // BAM CHANNELS
    def all_barcode_bams = Channel
        .from(all_libraries)
        .flatMap { lib_name, primers, bam_files ->
            bam_files.collect { file(it) }
        }

    def all_bams = Channel
        .from(all_libraries)
        .flatMap { lib_name, primers, bam_files ->
            bam_files.collect { file(it) }
        }
        .collect()

    MAKE_MASTER_FOFN(all_bams)

    // SAMPLE AND PRIMER CHANNELS
    all_barcode_bams
        .map { bam ->
            def sample = bam.getFileName().toString().replaceFirst(/\.bam$/, '')
            tuple(sample, bam, 'barcode')
        }
        .set { per_barcode_tuples }

    def first_primers = file(all_libraries[0][1])
    Channel.value(first_primers).set { primers_channel }

    // REFINE
    REFINE(per_barcode_tuples, primers_channel)

    REFINE.out[0]
        .map { flnc_bam ->
            def sample = flnc_bam.getFileName().toString().replaceFirst(/_flnc\.bam$/, '')
            tuple(sample, flnc_bam)
        }
        .filter { s, _ -> s != 'bulk' }   // isocall works directly on per-barcode FLNC
        .set { flnc_per_barcode }

    // MAP FLNC (per barcode)
    PBMM2_FLNC(flnc_per_barcode, genome_ch)
    def flnc_mapped_ch = PBMM2_FLNC.out.flnc_mapped   // (sample, bam)

    // ISOCALL: prep reference, profile per sample, merge + call
    ISOCALL_PREP_ISOFORMS(file(params.gencode_gtf))

    ISOCALL_PROFILE(flnc_mapped_ch)

    // Collect all profiles into a single list before the joint call
    def all_profiles = ISOCALL_PROFILE.out.profile.collect()

    ISOCALL_MERGE_CALL(
        all_profiles,
        ISOCALL_PREP_ISOFORMS.out.ref_isoforms,
        genome_ch
    )

    // Emit the joint GTF as the 'bulk' isoform set (matches annotate_isoforms expectation)
    def isoform_gtf_ch = ISOCALL_MERGE_CALL.out.isoforms_gtf
        .map { gtf -> tuple('bulk', gtf) }

    // READ -> TRANSCRIPT MAP via featureCounts--
    // Temporary substitute for a native isocall read assignment output.
    // Runs featureCounts at transcript level against the isocall GTF on all
    // per-sample FLNC mapped BAMs and converts the CORE output to FLAIR
    // read-map format (transcript_id<TAB>read1,read2,...).

    FEATURECOUNTS_ISOCALL_READ_MAP(
        flnc_mapped_ch.map { s, bam -> bam }.collect(),
        ISOCALL_MERGE_CALL.out.isoforms_gtf
    )

    def stub_cluster_ch = FEATURECOUNTS_ISOCALL_READ_MAP.out.read_map
        .map { f -> tuple('bulk', f) }

    // group_file is not parsed by read_flair_map() in build_read_isoform_lookup.py;
    // count_matrix is a valid existing file that satisfies the input path requirement.
    def stub_group_ch = ISOCALL_MERGE_CALL.out.count_matrix
        .map { f -> tuple('bulk', f) }

    // VARIANT CALLING
    MAKE_GENE_BED(file(params.gencode_gtf))

    PREPARE_DV_BAM(flnc_mapped_ch, MAKE_GENE_BED.out.bed_named)

    def genome_pair_ch = Channel.value(tuple(file(params.genome), file(params.genome + '.fai')))

    if (effective_caller == 'clair3') {
        def clair3_bam_ch = PREPARE_DV_BAM.out.bam
            .map { sample, bam, bai, bed -> tuple(sample, bam, bai) }
        RUN_CLAIR3_RNA(clair3_bam_ch, genome_pair_ch, effective_clair3_model)
        FILTER_VCF(RUN_CLAIR3_RNA.out.vcf)
        FILTER_VCF.out.vcf
            .map { sample, vcf, tbi -> tuple(sample, vcf) }
            .set { vcf_ch }
    } else {
        RUN_DEEPVARIANT(PREPARE_DV_BAM.out.bam, genome_pair_ch)
        RUN_DEEPVARIANT.out.vcf
            .set { vcf_ch }
    }

    flnc_mapped_ch
        .collect()
        .map { _ -> file("${params.outdir}/preprocessing/pbmm2_flnc") }
        .set { bam_dir_ch }



  emit:
    mapped_bams           = flnc_mapped_ch
    bam_dir               = bam_dir_ch
    bam_pattern           = Channel.value('*/*.flnc.mapped.bam')
    vcf                   = vcf_ch
    dv_bams               = PREPARE_DV_BAM.out.bam
    isoform_gffs          = isoform_gtf_ch
    isoform_abundance     = ISOCALL_MERGE_CALL.out.count_matrix
    cluster_reports       = stub_cluster_ch
    collapse_groups       = stub_group_ch
    gene_bed              = MAKE_GENE_BED.out.bed_named
    sqanti_classification = Channel.empty()
}
