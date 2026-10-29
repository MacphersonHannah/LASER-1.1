.. _nextflow_workflow_details:

Nextflow Workflow Details
==========================

This page describes the full configuration options for LASER. Copy the example config that matches your platform from ``example_configs/`` (see `Example Configs`_ below) and edit it before running the pipeline. For a minimal working example, see :doc:`installation_and_requirements`.

----

Platform
--------

.. code-block:: groovy

   params.platform = 'pacbio'   // see platform table in installation_and_requirements

----

PacBio Inputs
-------------

Option A: Starting from Demultiplexed BAMs (``pacbio``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LASER never runs ``lima``/demultiplexing itself — this is always a step you run yourself, before pointing LASER at your data. Once demultiplexing is complete (see :doc:`installation_and_requirements`), point LASER at the output directory in your config:

.. code-block:: groovy

   params.library_dirs = ['/path/to/output_dir']

For multiple libraries, provide one output directory per library — each must contain its own BAMs and primer file. LASER iterates over every entry in ``library_dirs`` automatically, so no extra flag is needed:

.. code-block:: groovy

   params.library_dirs = ['/path/to/lib1_output', '/path/to/lib2_output']

.. note::
   Every library run in parallel increases runtime and required compute. LASER has been tested successfully on 40 ENCODE PacBio RSII samples at once, but the ceiling depends on your cluster and read depth.

Option B: Starting from Refined FLNC Reads (``pacbio_fastqentry``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use this platform if your reads have already been through ``isoseq3 refine`` upstream of LASER — for example, data downloaded from a repository or processed on a different system.

Point ``params.library_dirs`` at a directory containing one FLNC file per sample. Accepted formats: ``.bam``, ``.fastq.gz``, ``.fq.gz``, ``.fastq``, ``.fq``. The sample name is taken from the filename (extension and any trailing ``_flnc``/``.flnc`` stripped):

.. code-block:: text

   /path/to/flnc_data/
       sample1_flnc.bam         → sample name: sample1
       sample2.flnc.fastq.gz    → sample name: sample2

.. code-block:: groovy

   params.platform      = 'pacbio_fastqentry'
   params.library_dirs  = ['/path/to/flnc_data']

If your FLNC files are fastq/fq.gz rather than BAM, ``params.pacbio_template_bam`` must point to a valid PacBio BAM header template (one is provided in ``test_data/``):

.. code-block:: groovy

   params.pacbio_template_bam = "${projectDir}/test_data/pacbio_header_template.bam"

.. important::
   Leave this pointing at the bundled template — do not swap it for your own BAM. Converting fastq back to BAM needs a valid PacBio-style header (the ``@RG``/``@PG`` fields, movie/well identifiers, and other PacBio-specific metadata that ``pbmm2``/``isoseq3`` expect), or downstream steps error out. The bundled template supplies exactly that header without carrying any of your actual read data or sample-specific information.

The underlying conversion approach is based on `Convert2PacBioBAM <https://github.com/jmestret/Convert2PacBioBAM>`_ by Jorge Mestre Tomas — an unofficial tool for converting FASTA/FASTQ reads to a PacBio-format BAM, using a real PacBio BAM as a header template.

PacBio Processing Steps
~~~~~~~~~~~~~~~~~~~~~~~~

The ``pacbio`` platform runs the following steps. The ``pacbio_fastqentry`` platform enters at cluster2, skipping the first two:

.. code-block:: text

   [demux_pacbio.sh — run separately]
     ↓
   Refine (primer removal + poly-A trimming)       ← pacbio only
     ↓
   isoseq3 cluster2 (per sample and bulk)          ← pacbio_fastqentry enters here
     ↓
   pbmm2 align (ISOSEQ preset)
     ↓
   isoseq3 collapse (isoform collapse)
     ↓
   get_abundance (FL count per isoform)
     ↓
   SQANTI3 QC + ORFannotate
     ↓
   DeepVariant variant calling (or Clair3-RNA if configured)
     ↓
   WhatsHap phasing + haplotagging
     ↓
   Transcript summary extraction (read × allele × transcript counts)
     ↓
   Allelic effect testing (abTE / abTU)

.. note::
   ``cluster2`` (isoseq3 v4+) is used — this is significantly faster than the original ``cluster`` command and produces polished HiFi-quality transcripts directly. A separate polish step is not required.

----

ONT Inputs
----------

Option A: Starting from FASTQ (``ont_cdna``, ``ont_directrna``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LASER expects one concatenated FASTQ per sample named ``SAMPLE.fastq.gz``, where the filename becomes the sample name throughout all outputs.

.. code-block:: groovy

   params.ont_input_files = [
       '/path/to/patient1_thalamus.fastq.gz',
       '/path/to/patient2_thalamus.fastq.gz'
   ]

See :doc:`installation_and_requirements` for how to prepare and concatenate ONT FASTQ files.

For ``ont_cdna``, set ``params.ont_chemistry`` to ``r9`` or ``r10`` to select the correct minimap2 preset.

Option B: Starting from Unmapped BAMs (``ont_cdna_minimap_entry``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use this platform when your data is in unmapped BAM format and you want to skip pychopper entirely as it has already been run. Each BAM becomes one sample; the filename (minus ``.bam``) becomes the sample label.

.. code-block:: groovy

   params.platform = 'ont_cdna_minimap_entry'
   params.ont_input_files = [
       '/path/to/sample1.bam',
       '/path/to/sample2.bam'
   ]

All reads go through a single FLAIR stream for both isoform discovery and quantification/variant calling. No FASTQ concatenation is needed.

``ont_cdna`` Only: Pychopper Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``ont_directrna`` and ``ont_cdna_minimap_entry`` platforms skip pychopper entirely. For ``ont_cdna``, set the appropriate kit matching your library preparation:

.. code-block:: groovy

   params.pychopper_kit     = 'PCS114'   // R10.4.1 PCR cDNA barcoding
                                          // 'PCB114' for PCR cDNA barcoding 24-plex
                                          // 'PCS111' for older R9 PCR cDNA
   params.pychopper_rescued = false       // set true to include pychopper-rescued reads

----

References
----------

.. code-block:: groovy

   params.genome      = '/path/to/GRCh38_noY_noMT.fa'
   params.gencode_gtf = '/path/to/gencode.v47_noY_noMT.annotation.gtf'
   params.cage        = '/path/to/refTSS_v3.3_human_coordinate.hg38.bed'
   params.intropolis  = '/path/to/intropolis.v1.hg19_with_liftover_to_hg38.min_count_10.tsv'
   params.polya_atlas = '/path/to/atlas.clusters.2.0.GRCh38.96.mod.bed'
   params.poly_motifs = '${projectDir}/assets/human.polyA.list.txt'

``params.genome`` must be the **primary assembly** FASTA (chromosomes + unplaced/unlocalised scaffolds, not alternate allele contigs) — this is what GENCODE/Ensembl label ``*.primary_assembly.fa``, and is what pbmm2/minimap2 and downstream variant calling expect.

See :doc:`key_outputs` for where to obtain each of these files and typical filenames.

Filtering Y and MT chromosomes out of the genome/GTF is optional but speeds things up:

.. code-block:: bash

   bash scripts/setup/strip_sex_mito_chroms.sh /path/to/genome.fa /path/to/annotation.gtf /path/to/output_prefix

This strips ``chrY``/``Y`` and ``chrM``/``MT`` from both files, and runs ``samtools faidx`` on the resulting FASTA automatically. If you're using a reference you didn't strip with this script, run ``samtools faidx`` on it yourself before starting a LASER run — several downstream steps assume a valid ``.fai`` index already exists alongside the FASTA.

----

Species
-------

LASER is species-agnostic at the pipeline level but requires the correct `CPAT <https://cpat.readthedocs.io/en/latest/>`_ model — so in practice, LASER can only be run on species CPAT supports. One param controls this for both ORFannotate and SQANTI3's species-specific behaviour:

.. code-block:: groovy

   params.species = 'human'   // default; or 'mouse' for mouse data

Supported values correspond to the species available in CPAT — see the `CPAT documentation <https://cpat.readthedocs.io/en/latest/>`_ for the full list. Only human and mouse data have been tested.

.. warning::
   Forgetting to set ``params.species`` for non-human data will score transcripts against the wrong hexamer model, and silently skip SQANTI3's junction-support annotation for human runs where it was never set.

For non-human data, see :ref:`notes_non_human`.

----

Variant Calling
---------------

PacBio
~~~~~~

DeepVariant is used by default and is strongly recommended for PacBio targeted runs. Clair3-RNA is available as an alternative.

.. code-block:: groovy

   params.variant_caller = 'deepvariant'   // default; or 'clair3' to override
   params.dv_container   = "${projectDir}/containers/deepvariant_1.9.0_sandbox"
   params.dv_model       = 'MASSEQ'        // for Iso-Seq / MAS-Seq
   params.clair3_model   = null            // only needed if variant_caller = 'clair3'
                                           // e.g. 'hifi_sequel2_pbmm2'

ONT
~~~

Clair3-RNA is used for all ONT platforms. Select the model matching your chemistry and basecaller:

.. code-block:: groovy

   params.clair3_rna_model = 'ont_r10_dorado_cdna'    // R10.4.1 cDNA with Dorado
                           // 'ont_dorado_drna004'      // direct RNA R10 with Dorado
                           // 'ont_r9_guppy_cdna'       // R9.4.1 cDNA with Guppy
                           // 'ont_guppy_drna002'        // direct RNA R9 with Guppy

To list all available models:

.. code-block:: bash

   singularity exec containers/clair3_rna_latest.sif \
       ls /opt/conda/envs/clair3_rna/bin/clair3_rna_models/

**Clair3 phasing**

.. code-block:: groovy

   params.clair3_enable_phasing = false   // default

Leave this ``false``. It controls Clair3-RNA's ``--enable_phasing_model``, which phases during variant calling. LASER phases separately with WhatsHap after calling, so enabling it duplicates the work.

**Minimum coverage**

.. code-block:: groovy

   params.clair3_min_coverage = 2   // default

Minimum read depth for Clair3-RNA to attempt a call. This sets how many heterozygous sites are available as phase anchors, so it propagates into how many reads WhatsHap can haplotag and therefore how many loci pass ``min_phased_reads_per_locus``. RNA coverage is very uneven across genes, so a low value is deliberate here; the locus-level read filters remove thinly-covered loci afterwards. Raise it if you suspect false heterozygous calls are mis-assigning reads between haplotypes.

----

WhatsHap Phasing Options
------------------------

.. code-block:: groovy

   params.phase_indels   = true    // include indels in phasing (recommended)
   params.phase_min_qual = 0       // recommended: keep at 0 -- raising this can filter
                                    // variants unevenly between alleles, biasing ASE calls
   params.hp_regions     = null    // optional BED file to restrict haplotagging region

----

Allelic Effect Testing
-----------------------

After transcript summary extraction, LASER runs two statistical tests per locus, per sample, on the read × allele × transcript count table:

* **abTE** (allelic bias in total expression): binomial test on allele-1 vs allele-2 total read counts.
* **abTU** (allelic bias in transcript usage): permutation Dirichlet-Multinomial test on the per-transcript allele count distribution — detects a shift in transcript repertoire between alleles, even when total expression is balanced.

Each locus in each sample is classified into one of four categories: no detectable allelic effect, abTE only, abTU without detectable abTE, or both. Loci that don't pass the eligibility thresholds (see :doc:`index`) are kept in the output with a ``tested_status`` explaining why, rather than silently dropped.

The Theta Parameter, and Which Mode to Use
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

abTU's Dirichlet-Multinomial test needs an overdispersion parameter, ``theta``, reflecting how noisy transcript-usage counts are beyond simple multinomial sampling. ``params.theta_mode`` controls how it's set:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Mode
     - When to use it
   * - ``calibrate_pooled`` (default)
     - You have at least 3 samples (2 works but is not ideal — see below) and want a theta tailored to your own data's actual noise characteristics. Runs a permutation-based FPR calibration pooled across every sample in the run, once, and applies the result to every sample.
   * - ``calibrate_pooled_by_group``
     - You have distinct sample groups (e.g. tissues) and want a separate theta calibrated for each, rather than one theta pooled across the whole run or reusing an old value. Runs the exact same calibration as ``calibrate_pooled``, once per group, using only that group's own samples. Requires ``theta_group_table`` — see `Defining Groups for calibrate_pooled_by_group`_ below.
   * - ``fixed``
     - You already know the right theta — from your own prior calibration, for a single-sample run where pooling isn't possible, or you just want to test a specific value. Supply ``theta_fixed`` (one value for every sample) or ``theta_table`` (CSV: ``sample,theta`` for per-sample values).
   * - ``lift``
     - You want to reuse an already-calibrated reference value instead of recalibrating — e.g. a single-sample run, a quick pilot, or matching the noise model of a companion dataset. Looks up a value from ``assets/calibrated_thetas_reference.csv`` by tissue/dataset name (``theta_lift_source``, default ``all_encode``), with optional per-sample override via ``theta_lift_table`` (CSV: ``sample,source``) if different samples in one run should lift from different sources (e.g. brain samples from ``brain_dlpfc``, heart samples from ``heart``).

``assets/calibrated_thetas_reference.csv`` currently contains:

.. code-block:: text

   source,theta,achieved_fpr,n_loci,n_samples_pooled,calibration_date,laser_version,notes
   brain_dlpfc,0.0119,0.0449,,,,1.1,From prior interactive DLPFC calibration (N=302650 permutations total). n_loci/n_samples_pooled/date not recorded from that run -- fill in if known.
   heart,PLACEHOLDER,,,,,1.1,Not yet calibrated in pipeline form. Run CALIBRATE_THETA_POOLED on heart-only samples (or lift from an existing heart run) and fill this row in.
   all_encode,PLACEHOLDER,,,,,1.1,Not yet calibrated in pipeline form. Run CALIBRATE_THETA_POOLED on the full 9-sample ENCODE set and fill this row in.

Anyone can append their own row. See ``assets/README_calibrated_thetas.md``.

.. note::
   **Why** ``calibrate_pooled`` **needs more than one sample:** the calibration pools observed loci across every sample in the run to estimate the empirical false-positive rate at each candidate theta. With a single sample there's nothing to pool — the pipeline refuses to run (``theta_mode = 'calibrate_pooled'`` errors out immediately if fewer than 2 samples are present, before any calibration work starts). With exactly 2 samples it proceeds, but with a warning: the resulting theta rests on a very small basis and may not generalise well. For a single sample, or a small handful you're not confident represents your population well, use ``fixed`` (if you already have a value) or ``lift`` (to borrow an already-validated one) instead.

Defining Groups for calibrate_pooled_by_group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Groups are defined by a two-column CSV, ``sample,group``, nothing more. Every sample in the run needs exactly one row. If any are missing, the pipeline fails immediately, before any calibration runs, rather than silently misassigning one.

.. code-block:: text

   sample,group
   ENCFF156TTD,brain
   ENCFF838DFB,brain
   ENCFF446EFU,heart
   ENCFF708BOP,other_peripheral

The simplest way to build this table is to keep your fastqs/BAMs organised one directory per group (which ``library_dirs`` already expects if you're running multiple libraries) and generate the table straight from those directories. This one-liner uses the exact sample-naming logic LASER itself applies, so the names are guaranteed to match:

.. code-block:: bash

   {
     echo "sample,group"
     for pair in \
       "/path/to/brain_dir:brain" \
       "/path/to/heart_dir:heart" \
       "/path/to/other_dir:other_peripheral"
     do
       dir="${pair%%:*}"; grp="${pair##*:}"
       for f in "$dir"/*; do
         sample=$(basename "$f" | sed -E 's/\.(bam|fastq\.gz|fq\.gz|fastq|fq)$//; s/[._-]?flnc$//')
         echo "${sample},${grp}"
       done
     done
   } > assets/theta_group_table_<run_name>.csv

A worked example from a real run is already in the repo at ``assets/theta_group_table_all_samples.csv``. Then point the config at it:

.. code-block:: groovy

   theta_mode        = 'calibrate_pooled_by_group'
   theta_group_table = "${projectDir}/assets/theta_group_table_<run_name>.csv"

.. note::
   Group size follows the same rule as ``calibrate_pooled``: each group needs at least 2 samples to pool a calibration across — a single sample has nothing to pool against. Groups with exactly 2 get a warning, not an error. Groups with fewer than 2 fail the run before any calibration starts.

abTE / abTU Test Parameters, Explained
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 15 55
   :header-rows: 1

   * - Param
     - Default
     - Meaning
   * - ``abte_fdr_threshold``
     - ``0.05``
     - BH-adjusted p-value cutoff for calling a locus abTE-significant.
   * - ``abte_effect_guard``
     - ``true``
     - If ``true``, statistical significance alone isn't enough: the locus must also clear ``abte_major_frac_threshold`` to be called significant. If ``false``, only the FDR threshold matters.
   * - ``abte_major_frac_threshold``
     - ``0.55``
     - The effect-size guard itself, only enforced when ``abte_effect_guard = true``: the major allele must carry at least this fraction of phased reads (e.g. ``0.55`` means at least a 55/45 split).
   * - ``abte_fdr_scope``
     - ``'sample'``
     - BH correction is computed within each sample's own p-values. This is the only scope the pipeline supports (each ``RUN_ALLELIC_EFFECT_TEST`` process only ever sees one sample's loci at a time), and there are no plans to add a run-wide alternative.
   * - ``abtu_fdr_threshold``
     - ``0.05``
     - Same as ``abte_fdr_threshold``, for abTU.
   * - ``abtu_effect_guard``
     - ``false``
     - Same concept as ``abte_effect_guard``, for abTU. Note the default is the opposite (off by default here, on by default for abTE). If ``true``, requires ``abtu_tvd_threshold`` to also be cleared.
   * - ``abtu_tvd_threshold``
     - ``0.20``
     - The effect-size guard for abTU, only enforced when ``abtu_effect_guard = true``: total variation distance between the two alleles' transcript-usage proportions must be at least this.
   * - ``abtu_perm_b``
     - ``1000``
     - Number of permutation draws used to build the null distribution for the abTU Dirichlet-Multinomial test, per locus. Higher = more precise p-values at the cost of runtime; the minimum resolvable p-value is roughly ``1/(B+1)``.
   * - ``dm_perm_seed``
     - ``123``
     - Random seed for the abTU permutation draws. Set once per run for reproducibility; not something you normally need to change.

Key Params
~~~~~~~~~~

.. code-block:: groovy

   theta_mode = 'calibrate_pooled'   // 'calibrate_pooled' | 'calibrate_pooled_by_group' | 'fixed' | 'lift'

   // calibrate_pooled_by_group
   theta_group_table = null          // CSV: sample,group -- required if theta_mode = 'calibrate_pooled_by_group'

   // fixed
   theta_fixed = null
   theta_table = null                // CSV: sample,theta

   // lift
   theta_reference_csv = "${projectDir}/assets/calibrated_thetas_reference.csv"
   theta_lift_source   = 'all_encode'
   theta_lift_table    = null        // CSV: sample,source

   // calibrate_pooled
   theta_calib_max_loci    = 2000
   theta_calib_n_null_ref  = 500
   theta_calib_n_test      = 50
   theta_calib_nominal_fpr = 0.05
   theta_calib_seed        = 456

   // locus filtering (applied at the read-count-summary stage)
   min_transcripts_per_locus  = 2
   min_reads_per_transcript   = 1
   min_total_reads_per_locus  = 30
   min_phased_reads_per_locus = 30
   min_reads_per_haplotype    = 5
   max_unphased_prop          = 0.30

   // abTE / abTU test parameters
   abte_fdr_threshold        = 0.05
   abte_effect_guard         = true
   abte_major_frac_threshold = 0.55
   abte_fdr_scope            = 'sample'   // the only scope the pipeline supports; see explanation above
   abtu_fdr_threshold = 0.05
   abtu_effect_guard  = false
   abtu_tvd_threshold = 0.20
   abtu_perm_b        = 1000
   dm_perm_seed       = 123

The consequence module also has a handful of optional override params (ORFannotate column names, permutation count, FDR threshold) not listed above, because they already default correctly inside the R scripts themselves — a normal config doesn't need to set any of them. See ``run_consequence_test.R`` and ``build_transcript_features.R`` for the full list, if you ever need to override one (e.g. your ORFannotate version names a column differently).

.. _example_configs:

Example Configs
~~~~~~~~~~~~~~~

See ``example_configs/`` for complete filled-in examples covering each platform:

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Example config
     - Platform
   * - ``example_cDNA_ONT.config``
     - ``ont_cdna``
   * - ``example_ont_cdna_minimap_entry.config``
     - ``ont_cdna_minimap_entry``
   * - ``example_ONT_directRNA.config``
     - ``ont_directrna``
   * - ``example_pacbio_not_fastq.config``
     - ``pacbio``
   * - ``example_pacbiofastq_human.config``
     - ``pacbio_fastqentry`` (human)
   * - ``example_pacbiofastq_mouse.config``
     - ``pacbio_fastqentry`` (mouse)

----

Resources
---------

.. code-block:: groovy

   process {
     cpus   = 16
     memory = '120 GB'
     time   = '100h'

     withName: 'RUN_DEEPVARIANT' {
       maxForks  = 7
       cpus      = 32
     }
     withName: 'RUN_CLAIR3_RNA' {
       container = "${projectDir}/containers/clair3_rna_latest.sif"
       maxForks  = 4
       cpus      = 16
       memory    = '32 GB'
     }

     // Allelic effect testing processes are single-threaded R
     withName: 'PREPARE_LOCI_TABLES|VALIDATE_THETA_SAMPLE_COUNT|VALIDATE_THETA_GROUP_SAMPLE_COUNT|CALIBRATE_THETA_POOLED|CALIBRATE_THETA_POOLED_BY_GROUP|RESOLVE_FIXED_THETA|LIFT_THETA|RUN_ALLELIC_EFFECT_TEST|MERGE_RUN_MANIFEST' {
       cpus   = 1
       memory = '8 GB'
       time   = '12h'
     }
   }

   executor {
     name      = 'local'
     cpus      = 250
     memory    = '900 GB'
     queueSize = 200
   }
