// Identical to pacbio_preprocess.nf from CLUSTER onward. Only the first part differs: input is already-refined FLNC reads (bam or fastq), so DEMULTIPLEX, VALIDATE_INPUTS
// and REFINE are skipped, and no primer fasta is required. Created because a lot of ENCODE samples are only available already slighlty processed
// 'bulk' is built the same way MAKE_MASTER_FOFN always built it - a plain fofn of per-sample FLNC bams handed to CLUSTER as-is, since cluster2
// accepts an flnc.fofn directly. No merge step needed.
// Variant calling: DeepVariant by default as recommended for PacBio but Clair3-RNA possible if params.variant_caller = 'clair3'.

nextflow.enable.dsl = 2

include { FLNC_FASTQ_TO_BAM }  from '../processes/identify_loci/flnc_fastq_to_bam.nf'
include { MAKE_MASTER_FOFN }   from '../processes/identify_loci/make_master_fofn.nf'
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

workflow pacbio_fastqentry_preprocess {

  main:

    log.info "Project dir: ${projectDir}"
    log.info "Scanning library directories... [PacBio FLNC input, cluster2 path]"

    def effective_caller       = params.variant_caller ?: 'deepvariant'
    def effective_dv_model     = params.dv_model ?: 'MASSEQ'
    def effective_clair3_model = params.clair3_model ?: 'hifi_sequel2_minimap2'
    log.info "Variant caller: ${effective_caller}"

    def library_dirs = (params.library_dirs instanceof List) ? params.library_dirs : [params.library_dirs]
    def flnc_files = []

    library_dirs.each { lib_dir ->
        def lib_path = file(lib_dir)
        if (!lib_path.exists()) error "Library directory not found: ${lib_path}"

        def found = lib_path.listFiles().findAll {
            it.name ==~ /.*\.(bam|fastq\.gz|fq\.gz|fastq|fq)$/
        }
        if (found.isEmpty()) error "No FLNC bam/fastq files found in ${lib_path}"
        flnc_files.addAll(found)
    }

    log.info "Found ${flnc_files.size()} FLNC file(s)."

    def genome_ch = Channel.value(file(params.genome))

    // PER-SAMPLE FLNC INPUT (no demux/refine already done upstream) 
    Channel.from(flnc_files)
        .map { f ->
            def sample = f.getFileName().toString()
                .replaceFirst(/\.(bam|fastq\.gz|fq\.gz|fastq|fq)$/, '')
                .replaceFirst(/[._-]?flnc$/, '')
            tuple(sample, file(f))
        }
        .set { raw_flnc_ch }

    FLNC_FASTQ_TO_BAM(raw_flnc_ch)
    def flnc_bam_ch = FLNC_FASTQ_TO_BAM.out   // (sample, flnc_bam) same shape REFINE.out had

    // 'bulk' same MAKE_MASTER_FOFN mechanism as always, pointed at the
    // already-refined per-sample bams instead of raw per-barcode bams 
    MAKE_MASTER_FOFN(flnc_bam_ch.map { s, bam -> bam }.collect())

    MAKE_MASTER_FOFN.out
        .map { fofn -> tuple('bulk', fofn) }
        .set { bulk_tuple }

    def flnc_for_map = flnc_bam_ch.mix(bulk_tuple) 

    def flnc_for_map_barcode = flnc_for_map.filter { s, ignored -> s != 'bulk' }

    // EVERYTHING BELOW IS UNCHANGED FROM pacbio_preprocess.nf -

    // MAP FLNC 
    PBMM2_FLNC(flnc_for_map_barcode, genome_ch)
    def flnc_mapped_ch = PBMM2_FLNC.out.flnc_mapped

    // CLUSTER 
    flnc_for_map
        .map { sample, flnc_bam ->
            def mode = sample == 'bulk' ? 'bulk' : 'barcode'
            tuple(sample, flnc_bam, mode)
        }
        .set { cluster_inputs }

    CLUSTER(cluster_inputs)

    // MAP HQ 
    CLUSTER.out.clustered_bam
        .map { bam ->
            def sample = bam.getBaseName().replaceFirst(/_clustered\.bam$/, '')
            def mode   = sample == 'bulk' ? 'bulk' : 'barcode'
            tuple(sample, bam, mode)
        }
        .set { map_inputs }

    PBMM2_HQ(map_inputs, genome_ch)
    PBMM2_HQ.out.set { bam_ch }

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
        .filter { s, ignored -> s != 'bulk' }
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

    flnc_mapped_ch
        .collect()
        .map { ignored -> file("${params.outdir}/preprocessing/pbmm2_flnc") }
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