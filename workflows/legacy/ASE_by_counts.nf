nextflow.enable.dsl = 2

include { ASE_COUNTS_AND_TESTS } from '../processes/ASE_by_counts/ASE_by_counts_process.nf'
include { ADD_BINOMIAL_TESTS }  from '../processes/ASE_by_counts/add_binomial_tests.nf'
include { UNIQUE_TRANSCRIPTS } from '../processes/ASE_by_counts/unique_transcripts.nf'

workflow ASE_by_counts {
  take:
    passed_success_csv  // successful loci from identify loci
    gene_bed            // gene BED (name in 4th column)
    lookup_db

  main:
 
    ASE_COUNTS_AND_TESTS(lookup_db, passed_success_csv, gene_bed)

  emit:
    per_sample   = ASE_COUNTS_AND_TESTS.out.all_samples
    //with_tests   = ADD_BINOMIAL_TESTS.out.with_tests
    //unique_tx    = UNIQUE_TRANSCRIPTS.out.uniq_tsv           // per-barcode
    //all_unique   = file('haplotype_unique_transcripts.tsv')  // merged
}


