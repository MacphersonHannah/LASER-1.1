.. _key_outputs:

Key Outputs
===========

This page describes the files produced by LASER and the external annotation files required to run the pipeline.

----

Pipeline Outputs
----------------

Results are written to ``run_{RUN_NAME}/{RUN_NAME}/``. The config file used is automatically copied into this directory for reproducibility.

.. list-table::
   :widths: 45 55
   :header-rows: 1

   * - File
     - Description
   * - ``annotate_isoforms/transcript_summary/transcript_counts_by_library.csv``
     - Per-transcript read counts across all samples and alleles
   * - ``annotate_isoforms/transcript_summary/read_counts_summary.csv``
     - Per-sample allele-level read counts
   * - ``annotate_isoforms/transcript_summary/transcript_annotations.gtf``
     - GTF with SQANTI and ORF annotations
   * - ``annotate_isoforms/SQANTI_QC/``
     - Full SQANTI3 QC output
   * - ``annotate_isoforms/ORFANNOTATE/``
     - ORF annotation output
   * - ``identify_loci/whatshap/<sample>/``
     - Phased VCFs and haplotagged BAMs per sample
   * - ``allelic_effect_testing/<sample>/allelic_effect_results.tsv``
     - Per-locus abTE/abTU p-values, FDR, four-way classification, and tested/not-tested status for that sample
   * - ``allelic_effect_testing/overview/allelic_effect_summary.tsv``
     - Run-level summary: theta used per sample, thresholds, significant-locus counts
   * - ``allelic_effect_testing/overview/consequence_test_summary.tsv``
     - Consequence-test summary: functional characterisation of abTU-positive loci across samples
   * - ``allelic_effect_testing/overview/transcript_features.tsv``
     - Per-transcript ORF and coding features used by the consequence test
   * - ``allelic_effect_testing/_theta/*.json``
     - Theta provenance: calibration diagnostics, lift source, or fixed value used

See :doc:`nextflow_workflow_details` for details on the allelic effect testing stage that produces the ``allelic_effect_testing/`` outputs.

----

Required Annotation Files
--------------------------

The following external files must be downloaded and provided via your config file before running LASER. See :doc:`nextflow_workflow_details` for the corresponding config parameters.

.. list-table::
   :widths: 20 40 40
   :header-rows: 1

   * - Parameter
     - Where to get it
     - Typical filename
   * - ``genome``
     - `GENCODE downloads <https://www.gencodegenes.org/human/>`_ — primary assembly FASTA
     - ``GRCh38.primary_assembly.genome.fa`` (or ``GRCh38_noY_noMT.fa`` if stripped)
   * - ``gencode_gtf``
     - `GENCODE downloads <https://www.gencodegenes.org/human/>`_, matching release
     - ``gencode.v47_noY_noMT.annotation.gtf``
   * - ``cage``
     - `refTSS database <http://reftss.clst.riken.jp/reftss/Main_Page>`_
     - ``refTSS_v3.3_human_coordinate.hg38.bed``
   * - ``intropolis``
     - `Intropolis <https://github.com/nellore/intropolis>`_
     - ``intropolis.v1.hg19_with_liftover_to_hg38.min_count_10.tsv``
   * - ``polya_atlas``
     - `PolyA_DB / PolyASite atlas <https://polyasite.unibas.ch/>`_
     - ``atlas.clusters.2.0.GRCh38.96.mod.bed``
   * - ``poly_motifs``
     - **Bundled with LASER** at ``assets/human.polyA.list.txt``; also available in SQANTI3 at ``data/polyA_motifs/human.polyA.list.txt``
     - ``human.polyA.list.txt``

``genome`` must be the primary assembly FASTA (chromosomes + unplaced/unlocalised scaffolds, not alternate allele contigs). For non-human species, species-matched equivalents of these files are required. See :ref:`notes_non_human` for details.
