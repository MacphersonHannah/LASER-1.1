nextflow.enable.dsl=2

//include { MAP_COUNTS }     from '../processes/reporting/map_counts.nf'
//include { SQANTI_SUMMARY } from '../processes/reporting/sqanti_summary.nf'

//workflow reporting {
//   take:
//     bams                      // tuple val(sample), path(bam), path(bai)
//     gene_bed                  // path
//     sqanti_classification     // path

//   main:
//     // per-sample/library counts
//     MAP_COUNTS(bams, gene_bed)

//     // once-per-run SQANTI summary
//     SQANTI_SUMMARY(sqanti_classification)

//   emit:
//     map_counts       = MAP_COUNTS.out
//     sqanti_by_cat    = SQANTI_SUMMARY.out[0]
//     sqanti_gene_cat  = SQANTI_SUMMARY.out[1]
//     sqanti_summary   = SQANTI_SUMMARY.out[2]
// }
