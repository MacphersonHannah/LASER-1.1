// workflows/pacbio_isoquant_preprocess.nf
nextflow.enable.dsl = 2

include { MAKE_MASTER_FOFN }    from '../processes/identify_loci/make_master_fofn.nf'
include { REFINE }              from '../processes/identify_loci/refine.nf'
include { MINIMAP2 }           from '../processes/identify_loci/minimap2.nf'
include { PBMM2_FLNC }           from '../processes/identify_loci/pbmm2_flnc.nf'
include { PREPARE_VARIANT_BAM }   from '../processes/identify_loci/prepare_variant_bam.nf'
include { RUN_DEEPVARIANT }    from '../processes/identify_loci/run_deepvariant.nf'
include { MAKE_GENE_BED }      from '../processes/identify_loci/make_gene_bed.nf'
include { SPLIT_BAM_BY_CHROM }   from '../processes/identify_loci/isoquant_parallel.nf'
include { ISOQUANT_CHROM }       from '../processes/identify_loci/isoquant_parallel.nf'
include { MERGE_ISOQUANT }       from '../processes/identify_loci/isoquant_parallel.nf'
include { SQANTI_QC }          from '../processes/annotate_isoforms/sqanti_qc.nf'

workflow pacbio_isoquant_preprocess {

  main:

    log.info "Project dir: ${projectDir}"
    log.info "Scanning library directories... [PacBio + IsoQuant mode]"

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
    def gtf_ch    = Channel.value(file(params.gencode_gtf))
    def db_ch     = Channel.value(file(params.gencode_db))

    def chromosomes = Channel.of(
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
        'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
        'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
        'chr19', 'chr20', 'chr21', 'chr22', 'chrX'
    )

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

    def first_primers = file(all_libraries[0][1])
    Channel.value(first_primers).set { primers_channel }

    // REFINE (barcode only) 
    REFINE(per_barcode_tuples, primers_channel)

    REFINE.out[0]
        .map { flnc_bam ->
            def sample = flnc_bam.getFileName().toString().replaceFirst(/_flnc\.bam$/, '')
            tuple(sample, flnc_bam)
        }
        .filter { s, ignored -> s != 'bulk' }
        .set { flnc_for_map }

    // MAP FLNC 
    PBMM2_FLNC(flnc_for_map, genome_ch)

    PBMM2_FLNC.out.flnc_mapped
        .set { flnc_mapped_ch }

    PBMM2_FLNC.out.flnc_mapped
        .map { sample, bam ->
            def bai = file("${bam}.bai")
            tuple(sample, bam, bai)
        }
        .set { flnc_mapped_with_bai }

    // GENE BED 
    MAKE_GENE_BED(file(params.gencode_gtf))

    // PREPARE DV BAM 
    PREPARE_VARIANT_BAM(flnc_mapped_with_bai, MAKE_GENE_BED.out.bed_named)

    // PARALLELISED ISOQUANT 
    SPLIT_BAM_BY_CHROM(flnc_mapped_with_bai, chromosomes)
    
    SPLIT_BAM_BY_CHROM.out
      .groupTuple(by: 0)
      .map { chrom, samples, bams, bais -> tuple(chrom, samples, bams, bais) }
      .set { isoquant_chrom_input }
    
    ISOQUANT_CHROM(isoquant_chrom_input, genome_ch,  db_ch, 'forward', 'pacbio')
    
    MERGE_ISOQUANT(
        ISOQUANT_CHROM.out.gtf.map { chrom, file -> file }.collect(),
        ISOQUANT_CHROM.out.counts.map { chrom, file -> file }.collect(),
        ISOQUANT_CHROM.out.gene_counts.map { chrom, file -> file }.collect(),
        ISOQUANT_CHROM.out.model_reads.map { chrom, file -> file }.collect(),
        ISOQUANT_CHROM.out.sqanti_classification.map { chrom, file -> file }.collect(),
        ISOQUANT_CHROM.out.read_assignments.map { chrom, file -> file }.collect(),
        ISOQUANT_CHROM.out.discovered_counts.map { chrom, file -> file }.collect()
)

    // SQANTI QC (real SQANTI3 on the IsoQuant GTF) 
    def sqanti_input_ch = MERGE_ISOQUANT.out.transcript_models_gtf
        .map { gtf -> tuple('bulk', gtf) }
        .join(
            MERGE_ISOQUANT.out.combined_counts
                .map { tsv -> tuple('bulk', tsv) }
        )
        .map { s, gtf, tsv -> tuple(s, gtf, tsv) }

    SQANTI_QC(
        sqanti_input_ch,
        file(params.gencode_gtf),
        file(params.genome),
        file(params.cage),
        file(params.intropolis),
        file(params.polya_atlas),
        file(params.poly_motifs)
    )

    SQANTI_QC.out.classification
        .set { sqanti_classification_ch }

    MERGE_ISOQUANT.out.model_reads
        .map { tsv -> tuple('bulk', tsv) }
        .set { cluster_reports_ch }

    MERGE_ISOQUANT.out.model_reads
        .map { tsv -> tuple('bulk', tsv) }
        .set { collapse_groups_ch }

    flnc_mapped_ch
        .collect()
        .map { ignored -> file("${params.outdir}/preprocessing/pbmm2_flnc") }
        .set { bam_dir_ch }

    // VARIANT CALLING
    def genome_pair_ch = Channel.value(tuple(file(params.genome), file(params.genome + '.fai')))
    RUN_DEEPVARIANT(PREPARE_VARIANT_BAM.out.bam, genome_pair_ch)

  emit:
    mapped_bams           = flnc_mapped_ch
    bam_dir               = bam_dir_ch
    dv_bams               = PREPARE_VARIANT_BAM.out.bam
    vcf                   = RUN_DEEPVARIANT.out.vcf
    isoform_gffs          = SQANTI_QC.out.corrected        // SQANTI-corrected GTF for ORFANNOTATE
    isoform_abundance     = MERGE_ISOQUANT.out.combined_counts
    cluster_reports       = cluster_reports_ch
    collapse_groups       = collapse_groups_ch
    sqanti_classification = sqanti_classification_ch        // real SQANTI3 classification
    gene_bed              = MAKE_GENE_BED.out.bed_named
    bam_pattern           = Channel.value('*/*.mapped.bam')
}
