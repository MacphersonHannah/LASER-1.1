include { PREPARE_LOCI_TABLES }               from '../processes/allelic_effect/prepare_loci_tables.nf'
include { VALIDATE_THETA_SAMPLE_COUNT }       from '../processes/allelic_effect/validate_theta_sample_count.nf'
include { VALIDATE_THETA_GROUP_SAMPLE_COUNT } from '../processes/allelic_effect/validate_theta_group_sample_count.nf'
include { CALIBRATE_THETA_POOLED }            from '../processes/allelic_effect/calibrate_theta_pooled.nf'
include { CALIBRATE_THETA_POOLED_BY_GROUP }   from '../processes/allelic_effect/calibrate_theta_pooled_by_group.nf'
include { RESOLVE_FIXED_THETA }               from '../processes/allelic_effect/resolve_fixed_theta.nf'
include { LIFT_THETA }                        from '../processes/allelic_effect/lift_theta.nf'
include { RUN_ALLELIC_EFFECT_TEST }           from '../processes/allelic_effect/run_allelic_effect_test.nf'
include { MERGE_RUN_MANIFEST }                from '../processes/allelic_effect/merge_run_manifest.nf'
include { BUILD_TRANSCRIPT_FEATURES }         from '../processes/allelic_effect/build_transcript_features.nf'
include { RUN_CONSEQUENCE_TEST }              from '../processes/allelic_effect/run_consequence_test.nf'
include { MERGE_CONSEQUENCE_MANIFEST }        from '../processes/allelic_effect/merge_consequence_manifest.nf'

workflow ALLELIC_EFFECT_TESTING {
  take:
    read_counts_csv   // path: combined read_counts_summary.csv (all samples)
    gene_loci_bed     // path: gene_loci_named.bed
    orf_summary       // path: bulk ORFannotate_summary.tsv

  main:
    PREPARE_LOCI_TABLES(read_counts_csv, gene_loci_bed)

    // Separate per-sample files into (sample, path) channel
    per_sample_ch = PREPARE_LOCI_TABLES.out.per_sample_rds
      .flatten()
      .map { f -> tuple(f.baseName.replaceFirst(/\.loci$/, ''), f) }

    // Comma sample list, needed by fixed/lift/calibrate_pooled_by_group and by the sample count guards below.
    sample_list_ch = per_sample_ch
      .map { sample, f -> sample }
      .collect()
      .map { samples -> samples.join(',') }

    // Choosing a theta: every mode converges on the same theta.json shape so RUN_ALLELIC_EFFECT_TEST never needs to know which mode produced it
    if (params.theta_mode == 'calibrate_pooled') {
      // Fail fast, before any permutations if there aren't enough samples to meaningfully pool a theta across.
      VALIDATE_THETA_SAMPLE_COUNT(sample_list_ch)
      CALIBRATE_THETA_POOLED(
        PREPARE_LOCI_TABLES.out.combined_annotated_rds,
        VALIDATE_THETA_SAMPLE_COUNT.out.ok
      )
      theta_ch = CALIBRATE_THETA_POOLED.out.theta_value
    } else if (params.theta_mode == 'calibrate_pooled_by_group') {
      // Same calibration as above but run once per named group (eg brain/heart/other_peripheral) instead of pooling every sample in the run into one theta. 
      VALIDATE_THETA_GROUP_SAMPLE_COUNT(sample_list_ch, file(params.theta_group_table))
      CALIBRATE_THETA_POOLED_BY_GROUP(
        PREPARE_LOCI_TABLES.out.combined_annotated_rds,
        file(params.theta_group_table),
        VALIDATE_THETA_GROUP_SAMPLE_COUNT.out.ok
      )
      theta_ch = CALIBRATE_THETA_POOLED_BY_GROUP.out.theta_value
      // fixed is giving a theta to LASER
    } else if (params.theta_mode == 'fixed') {
      RESOLVE_FIXED_THETA(sample_list_ch)
      theta_ch = RESOLVE_FIXED_THETA.out.theta_value
      // lift is the same as fixed except you lift the number from a reference CSV. Thought it might be useful although the user will need to fill it themselves
    } else if (params.theta_mode == 'lift') {
      LIFT_THETA(sample_list_ch, file(params.theta_reference_csv))
      theta_ch = LIFT_THETA.out.theta_value
    } else {
      error "Unknown params.theta_mode '${params.theta_mode}'. Must be one of: calibrate_pooled, calibrate_pooled_by_group, fixed, lift"
    }

    // theta_ch carries a single value channel (one theta.json) - combine() broadcasts it onto every sample without needing a per-sample join.
    RUN_ALLELIC_EFFECT_TEST(
      per_sample_ch.combine(theta_ch)
    )

    MERGE_RUN_MANIFEST(
      RUN_ALLELIC_EFFECT_TEST.out.sample_manifest.collect(),
      theta_ch
    )

    // Consequence module: functional characterisation of abTU-positive loci only
    BUILD_TRANSCRIPT_FEATURES(orf_summary)

    conseq_in = per_sample_ch
      .join(RUN_ALLELIC_EFFECT_TEST.out.results_tsv)      // (sample, rds, results_tsv)
      .combine(theta_ch)                                  // + theta.json
      .combine(BUILD_TRANSCRIPT_FEATURES.out.features_tsv) // + transcript_features.tsv

    RUN_CONSEQUENCE_TEST(conseq_in)

    MERGE_CONSEQUENCE_MANIFEST(
      RUN_CONSEQUENCE_TEST.out.sample_manifest.collect()
    )

emit:
  per_sample_results   = RUN_ALLELIC_EFFECT_TEST.out.results_tsv
  run_manifest         = MERGE_RUN_MANIFEST.out.manifest_tsv
  consequence_results  = RUN_CONSEQUENCE_TEST.out.results_tsv
  consequence_manifest = MERGE_CONSEQUENCE_MANIFEST.out.manifest_tsv
}
