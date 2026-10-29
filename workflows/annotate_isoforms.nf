nextflow.enable.dsl = 2
import java.nio.file.Path
import groovy.json.JsonSlurper

include { SQANTI_QC }                                          from '../processes/annotate_isoforms/sqanti_qc.nf'
include { ORFANNOTATE }                                        from '../processes/annotate_isoforms/orf_annotate.nf'
include { BUILD_READ_ISOFORM_LOOKUP as BUILD_READ_ISOFORM_LOOKUP_PROC } from '../processes/annotate_isoforms/build_read_isoform_lookup.nf'
include { EXTRACT_HAPLOTYPE_TAGS }                             from '../processes/annotate_isoforms/extract_haplotype_tags.nf'
include { EXTRACT_TRANSCRIPT_SUMMARY }                         from '../processes/annotate_isoforms/extract_transcript_summary.nf'
include { ASSIGN_READS_TO_GENES }                              from '../processes/annotate_isoforms/assign_reads_to_genes.nf'

workflow annotate_isoforms {
  take:
    gff_files
    fl_count_files
    cluster_reports
    collapse_groups
    haplotagged_bams
    read_sample_map
    sqanti_classification

  main:
    def genome_ch = Channel.value(file(params.genome))


    def gff_ch = gff_files
      .map { sample, gff -> tuple(sample, gff) }
      .filter { s, _ -> s == 'bulk' }

    def fl_count_ch = fl_count_files
      .filter { s, _ -> s == 'bulk' }
    def sqanti_inputs = gff_ch
      .join(fl_count_ch)
      .map { s, gff, fl -> tuple(s, gff, fl) }

    SQANTI_QC(
      sqanti_inputs,
      file(params.gencode_gtf),
      file(params.genome),
      file(params.cage),
      file(params.polya_atlas),
      file(params.poly_motifs)
    )

    def orfannotate_inputs = SQANTI_QC.out.corrected.filter { s, _ -> s == 'bulk' }

    ORFANNOTATE(orfannotate_inputs, Channel.value(file(params.genome)))

    def bulk_sqanti_gene_map = SQANTI_QC.out.classification
  .filter { s, _ -> s == 'bulk' }
  .map { s, tsv -> tsv }

    def bulk_cluster = cluster_reports
      .filter { s, _ -> s == 'bulk' }

    def bulk_group = collapse_groups
      .filter { s, _ -> s == 'bulk' }

    def bulk_orf_sum = ORFANNOTATE.out.summary
      .filter { s, _ -> s == 'bulk' }

    def bulk_sqanti = SQANTI_QC.out.classification
      .filter { s, _ -> s == 'bulk' }

    def hp_bams_barcode = haplotagged_bams.filter { s, _ -> s != 'bulk' }

    EXTRACT_HAPLOTYPE_TAGS(hp_bams_barcode)

    def haplo_tsv_list =
      EXTRACT_HAPLOTYPE_TAGS.out
        .map { s, tsv -> tsv }
        .collect()
        .map { tsvs -> tuple('bulk', tsvs) }

    ASSIGN_READS_TO_GENES(
      haplotagged_bams,
      file(params.gencode_gtf)
    )

    def read2gene_tsv_list =
      ASSIGN_READS_TO_GENES.out.read2gene
        .map { s, tsv -> tsv }
        .collect()
        .map { tsvs -> tuple('bulk', tsvs) }

    def lookup_in = bulk_cluster
      .join(bulk_group)
      .join(bulk_orf_sum)
      .join(bulk_sqanti)
      .combine(read_sample_map)
      .join(haplo_tsv_list)
      .join(read2gene_tsv_list)
      .map { s, cluster, group, orf, sqanti, mapfile, haplist, fc_list ->
          tuple(s, cluster, group, orf, sqanti, mapfile, haplist, fc_list)
      }
      .view { "LOOKUP_IN: ${it[0]} | ${it[1].name} | ${it[2].name} | ${it[3].name} | ${it[4].name} | ${it[5].name} | ${it[6].size()} hap files | ${it[7].size()} fc files" }

    BUILD_READ_ISOFORM_LOOKUP_PROC(lookup_in, Channel.value(params.platform.toLowerCase()))

    def bulk_sqlite_one = BUILD_READ_ISOFORM_LOOKUP_PROC.out.sqlite
    def sqlite_ch       = BUILD_READ_ISOFORM_LOOKUP_PROC.out.sqlite

    EXTRACT_TRANSCRIPT_SUMMARY(sqlite_ch)

  emit:
    lookup_db       = BUILD_READ_ISOFORM_LOOKUP_PROC.out.sqlite
    read_counts_csv = EXTRACT_TRANSCRIPT_SUMMARY.out.summary_csv
    orf_summary     = ORFANNOTATE.out.summary.filter { s, _ -> s == 'bulk' }.map { s, f -> f }  // NEW
}
