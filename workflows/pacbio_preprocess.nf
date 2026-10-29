// PacBio preprocessing using isoseq3 tools (CLUSTER > POLISH > COLLAPSE_ISOFORMS).
// Variant calling: DeepVariant by default, Clair3-RNA if params.variant_caller = 'clair3'.
// DeepVariant is recommended for PacBio targeted runs.
nextflow.enable.dsl = 2
import java.nio.file.Path
import groovy.json.JsonSlurper

include { MAKE_MASTER_FOFN }   from '../processes/identify_loci/make_master_fofn.nf'
include { VALIDATE_INPUTS }    from '../processes/identify_loci/validate_inputs.nf'
include { REFINE }             from '../processes/identify_loci/refine.nf'
include { PBMM2_FLNC }        from '../processes/identify_loci/pbmm2_flnc.nf'
include { CLUSTER }            from '../processes/identify_loci/cluster.nf'
include { PBMM2_HQ }          from '../processes/identify_loci/pbmm2_hq.nf'
include { COLLAPSE_ISOFORMS }  from '../processes/identify_loci/collapse_isoforms.nf'
include { GET_ABUNDANCE }      from '../processes/identify_loci/get_abundance.nf'
include { PREPARE_DV_BAM }     from '../processes/identify_loci/prepare_dv.nf'
include { MAKE_GENE_BED }      from '../processes/identify_loci/make_gene_bed.nf'
include { RUN_DEEPVARIANT }    from '../processes/identify_loci/run_deepvariant.nf'
include { RUN_CLAIR3_RNA }     from '../processes/identify_loci/run_clair3_rna.nf'
include { FILTER_VCF }         from '../processes/identify_loci/filter_vcf.nf'

workflow pacbio_preprocess {

  main:

    log.info "Project dir: ${projectDir}"
    log.info "Scanning library directories... [PacBio mode]"

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

    // BAM AND FOFN CHANNELS
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

    MAKE_MASTER_FOFN.out
        .map { fofn -> tuple('bulk', fofn, 'bulk') }
        .set { bulk_tuple }

    def first_primers = file(all_libraries[0][1])
    Channel.value(first_primers).set { primers_channel }

    // REFINE 
    REFINE(per_barcode_tuples.mix(bulk_tuple), primers_channel)

    REFINE.out[0]
        .map { flnc_bam ->
            def sample = flnc_bam.getFileName().toString().replaceFirst(/_flnc\.bam$/, '')
            tuple(sample, flnc_bam)
        }
        .set { flnc_for_map }

    def flnc_for_map_barcode = flnc_for_map.filter { s, _ -> s != 'bulk' }

    // MAP FLNC 
    PBMM2_FLNC(flnc_for_map_barcode, genome_ch)
    def flnc_mapped_ch = PBMM2_FLNC.out.flnc_mapped   // (sample, bam)

    // CLUSTER
    flnc_for_map
        .map { sample, flnc_bam ->
            def mode = sample == 'bulk' ? 'bulk' : 'barcode'
            tuple(sample, flnc_bam, mode)
        }
        .set { cluster_inputs }

    CLUSTER(cluster_inputs)

    // MAP HQ (direct from cluster2 BAM no FASTA conversion needed) 
    CLUSTER.out.clustered_bam
        .map { bam ->
            def sample = bam.getBaseName().replaceFirst(/_clustered\.bam$/, '')
            def mode   = sample == 'bulk' ? 'bulk' : 'barcode'
            tuple(sample, bam, mode)
        }
        .set { map_inputs }

    PBMM2_HQ(map_inputs, genome_ch)

    PBMM2_HQ.out
        .set { bam_ch }

    // COLLAPSE AND ABUNDANCE 
    COLLAPSE_ISOFORMS(bam_ch)

    COLLAPSE_ISOFORMS.out[0]
        .map { gff ->
            def sample = gff.getParent().getName().replaceFirst(/_clustered$/, '')
            def prefix = "${gff.getParent().toString()}/${gff.getParent().getName()}.collapsed"
            tuple(sample, prefix)
        }
        .set { collapse_prefix_ch }

    CLUSTER.out.cluster_report
        .map { cluster_report ->
            def sample = cluster_report.getName()
                        .replaceFirst(/_clustered\.cluster_report\.csv$/, '')
            tuple(sample, cluster_report)
        }
        .filter { sample, report -> sample && report }
        .set { cluster_report_ch }

    COLLAPSE_ISOFORMS.out[1]
        .map { group_txt ->
            def sample = group_txt.getParent().getName().replaceFirst(/_clustered$/, '')
            tuple(sample, group_txt)
        }
        .set { collapse_groups_ch }

    collapse_prefix_ch
        .filter { sample, prefix -> sample && prefix }
        .join(cluster_report_ch)
        .map { sample, prefix, cluster_report -> tuple(sample, prefix, cluster_report) }
        .set { abundance_inputs }

    GET_ABUNDANCE(abundance_inputs)

    // GENE BED 
    MAKE_GENE_BED(file(params.gencode_gtf))

    flnc_mapped_ch
        .filter { s, _ -> s != 'bulk' }
        .set { prepare_dv_inputs }

    PREPARE_DV_BAM(prepare_dv_inputs, MAKE_GENE_BED.out.bed_named)

    // VARIANT CALLING 
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

    // Wait for all FLNC BAMs before emitting the directory path
    flnc_mapped_ch
        .collect()
        .map { _ -> file("${params.outdir}/preprocessing/pbmm2_flnc") }
        .set { bam_dir_ch }

  emit:
    mapped_bams           = flnc_mapped_ch
    bam_dir               = bam_dir_ch
    bam_pattern           = Channel.value('*/*.mapped.bam')
    vcf                   = vcf_ch
    dv_bams               = PREPARE_DV_BAM.out.bam
    isoform_gffs = COLLAPSE_ISOFORMS.out[0]
        .map { gff -> tuple(gff.getParent().getName().replaceFirst(/_clustered$/, ''), gff) }
    isoform_abundance     = GET_ABUNDANCE.out[1]
    cluster_reports       = cluster_report_ch
    collapse_groups       = collapse_groups_ch
    gene_bed              = MAKE_GENE_BED.out.bed_named
    sqanti_classification = Channel.empty()
    }
