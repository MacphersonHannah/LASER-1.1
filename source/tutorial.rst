.. _tutorial:

Tutorial: Mouse F1 Hybrid
==========================

This tutorial runs LASER end to end on a small, real dataset, then walks through a Jupyter notebook that checks each stage of the run and interprets the allelic results. It takes about an hour of hands-on time plus the pipeline run (up to 6 hours on a 16-CPU node).

All tutorial files live in the ``tutorial/`` folder of the `LASER repository <https://github.com/MacphersonHannah/LASER-1.1>`_:

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - File
     - Purpose
   * - ``tutorial/scripts/1_download_dataset.sh``
     - Downloads and unpacks the tutorial FASTQs and chr17 reference files into ``tutorial/data/``
   * - ``tutorial/scripts/tutorial.config``
     - Ready-made LASER config for the tutorial. All paths are relative to the LASER folder, so nothing needs editing after cloning
   * - ``tutorial/scripts/2_run_tutorial.sh``
     - SLURM submission script that runs LASER with ``tutorial.config``
   * - ``tutorial/scripts/0_tutorial.ipynb``
     - The analysis notebook: QC of every pipeline stage, then interpretation of the abTE/abTU results

----

The Dataset
-----------

The tutorial uses 8 ENCODE PacBio long-read libraries from the **B6 × CAST F1 hybrid mouse** (C57BL/6J × CAST/EiJ), provided as refined FLNC reads (``fastq.gz``) and subset to **chromosome 17**. References are mm10/GRCm38 with GENCODE vM10, also cut down to chr17.

The F1 hybrid is a standard model for allele-specific expression, and it makes the results easy to check:

* **Dense heterozygous sites.** Both parents are inbred, so every site where B6 and CAST differ is heterozygous in the F1. This gives plenty of SNVs to phase on.
* **Known allele identity.** The reference genome is B6, so within any phase block the haplotype carrying the reference alleles is B6 and the other is CAST. The notebook uses this to turn WhatsHap's arbitrary "haplotype 1/2" into "B6/CAST".
* **An interesting chromosome.** Chr17 is highly diverged between laboratory and wild-derived strains. It carries the polymorphic MHC region, the imprinted *Igf2r*/*Airn* cluster, and many cis-eQTLs.

.. note::
   The data are deliberately small (one chromosome, subsampled reads) so the run fits on a single node. Treat the results as a demonstration that the pipeline works, not as biological findings — see `Known Limitations of the Tutorial Data`_.

----

Before You Start
----------------

You will need:

* LASER installed, with ``setup_tools.sh`` and ``setup_envs.sh`` already run — see :doc:`installation_and_requirements`.
* Access to a SLURM compute node with **16 CPUs and ~54 GB of memory** for up to 6 hours. If your cluster isn't SLURM, you can run ``2_run_tutorial.sh`` directly with ``bash`` on a node of that size.
* A Python environment with Jupyter and the notebook's packages: ``pandas``, ``numpy``, ``matplotlib``, ``scipy`` and ``gtfparse``.

----

Step 1: Download the Data
--------------------------

.. code-block:: bash

   bash tutorial/scripts/1_download_dataset.sh

The script can be called from anywhere — it always unpacks into ``tutorial/data/``. It downloads ``laser_tutorial_data.zip`` from the LASER GitHub release, extracts it, deletes the zip, and verifies the files against ``data/checksums.md5`` if that file is present.

After it finishes you should have:

.. code-block:: text

   tutorial/data/
       fastq/          ← 8 FLNC fastq.gz files, one per sample
       reference/      ← chr17 genome (+ .fai), GENCODE vM10 GTF, refTSS, PolyA atlas, polyA motifs

.. note::
   The script uses ``wget``. If your system only has ``curl``, swap in the commented-out ``curl -L -o ...`` line in the script.

----

Step 2: Review the Config
--------------------------

``tutorial/scripts/tutorial.config`` is complete and should run as-is. The settings worth knowing about:

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - Setting
     - Why
   * - ``platform = 'pacbio_fastqentry'``
     - The reads are already-refined FLNC ``fastq.gz`` files, so LASER enters at ``cluster2`` and skips demultiplexing and refine. Every FASTQ in ``library_dirs`` becomes one sample (8 here).
   * - ``pacbio_template_bam``
     - Points at the bundled header template so the FASTQs can be converted back to PacBio BAM — leave this as it is.
   * - ``species = 'mouse'``
     - Selects the mouse CPAT model for ORFannotate and tells SQANTI3 not to apply human-only junction support.
   * - ``intropolis = ''``
     - There is no mouse equivalent of Intropolis, so it is left empty. ``cage`` and ``polya_atlas`` point at mm10, chr17-only versions.
   * - ``variant_caller = 'deepvariant'``, ``dv_model = 'MASSEQ'``
     - The standard PacBio choice; the MASSEQ model is not human-specific, so it is correct for mouse.
   * - ``theta_mode = 'calibrate_pooled'``
     - Mirrors a real multi-sample run. If the chr17 subset leaves too few loci to calibrate, switch to ``theta_mode = 'lift'``.

**Resources.** The ``process`` and ``executor`` blocks are sized to fit inside one 16-CPU / 54 GB SLURM job: the executor is capped at 16 CPUs and 48 GB, which leaves headroom for Nextflow itself. The largest single tasks (``CLUSTER``, ``SQANTI_QC``, ``ORFANNOTATE``) take 8 CPUs and 24 GB, and at most two DeepVariant jobs run at once. If you change the SLURM request in ``2_run_tutorial.sh``, change the two ``executor`` numbers to match — the executor must never promise more than the job holds.

----

Step 3: Run LASER
------------------

Before submitting, open ``tutorial/scripts/2_run_tutorial.sh`` and adapt the environment lines to your cluster:

.. code-block:: bash

   # Point this at YOUR conda install's conda.sh
   source /path/to/miniforge3/etc/profile.d/conda.sh
   conda activate base

   module load singularity      # or however Singularity/Apptainer is provided on your system
   export NXF_VER=25.04.8

.. important::
   Load conda through its ``conda.sh`` hook script, as above. Adding ``miniforge3/bin`` to ``PATH`` instead breaks the ORFannotate step.

Then submit **from the** ``tutorial/scripts/`` **folder** — the script finds ``run_LASER.sh`` and ``tutorial.config`` relative to where it is launched:

.. code-block:: bash

   cd tutorial/scripts
   sbatch 2_run_tutorial.sh

Outputs are written to ``run_tutorial/tutorial/`` in the LASER folder, with one folder per pipeline stage:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Stage folder
     - What LASER did there
   * - ``preprocessing/``
     - Clustering, alignment (pbmm2), isoform collapse, DeepVariant variant calling
   * - ``identify_loci/``
     - WhatsHap phasing and haplotagging, gene counts
   * - ``annotate_isoforms/``
     - SQANTI3, ORFannotate/CPAT, allele-resolved read counts
   * - ``allelic_effect_testing/``
     - abTE and abTU tests, theta calibration, consequence analysis
   * - ``report_*.html``, ``timeline_*.html``, ``trace_*.txt``
     - Nextflow run reports — open in a browser to see runtime and memory per process

See :doc:`key_outputs` for a description of the main result files.

----

Step 4: Open the Notebook
--------------------------

Once the run has finished, start Jupyter from ``tutorial/scripts/`` and open ``0_tutorial.ipynb``:

.. code-block:: bash

   cd tutorial/scripts
   jupyter lab 0_tutorial.ipynb

The notebook only reads files the tutorial run writes, so it should run top to bottom without changes. You can also read it on this site without running anything: see :doc:`tutorial_notebook`. The settings cell at the top controls where it looks:

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Setting
     - Meaning
   * - ``OUT_DIR``
     - The run's output folder. Defaults to ``../../run_tutorial/tutorial``; override with the ``LASER_OUT`` environment variable.
   * - ``REF_DIR``
     - The tutorial reference folder. Defaults to ``../../tutorial/data/reference``; override with ``LASER_REF``.
   * - ``CROSS_MATERNAL``
     - Which strain is the mother (``"B6"`` for ENCODE's B6CASTF1/J mice). Only used to predict the direction of imprinting controls.
   * - ``SAVE_OUTPUTS``
     - If ``True``, every figure and a per-library report card are written to ``notebook_outputs/``.

The thresholds the notebook uses (``min_phased_reads_per_locus``, ``abte_fdr_threshold``, ``abte_major_frac_threshold``) are read from the ``config_used_*.config`` file LASER copies into the output folder, so they always match the run.

----

What the Notebook Covers
-------------------------

The notebook is split into two halves: sections 0–5 check that each stage of the pipeline worked, and sections 6–11 interpret the allelic results. The full notebook, with the figures from the reference tutorial run, is on the next page — compare your own output against it.

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - Notebook section
     - What it checks
   * - 0. Did the run finish?
     - Every sample has its DeepVariant VCF, phased VCF, gene counts and ``allelic_effect_results.tsv``. A missing file usually means a process failed; the Nextflow ``report_*.html`` shows which one.
   * - 1. Variant calling (DeepVariant)
     - Mostly heterozygous calls, a transition/transversion (Ti/Tv) ratio well above the ~0.5 expected from random errors, and the same het sites shared across samples.
   * - 2. Phasing (WhatsHap)
     - Fraction of het SNVs phased, phase block N50 (should be longer than a typical gene), and block consistency — in an inbred F1 every ALT allele is CAST, so all ALT alleles in a correct block sit on one haplotype.
   * - 3. Read phasing and locus eligibility
     - Proportion of reads assigned to haplotype 1, haplotype 2 or unphased, and which loci pass the filters for testing.
   * - 4. Haplotype 1/2 to B6/CAST
     - Orients each gene's alleles as B6 or CAST (the haplotype carrying the reference alleles is B6), then checks for reference mapping bias.
   * - 5. Isoform QC (SQANTI3 + CPAT)
     - SQANTI3 structural categories (FSM, ISM, NIC, NNC and others) and CPAT coding probability.
   * - 6. abTE
     - LASER's allelic bias in total expression calls, plus the *Igf2r* imprinting control.
   * - 7. abTU
     - LASER's allelic bias in transcript usage calls, with total variation distance (TVD) as the effect size.
   * - 8. Why theta matters
     - A QQ plot of abTU p-values to check the test is calibrated. **Always check this before interpreting abTU.**
   * - 9. Four-way classification
     - How many loci are abTE only, abTU without abTE, both, or neither.
   * - 10. Functional consequence
     - For abTU loci, how the alleles' dominant transcripts differ: coding status, NMD, TSS, TES or ORF length.
   * - 11. Report card
     - One row per library combining all of the above, saved to ``notebook_outputs/report_card.tsv``.

.. toctree::
   :maxdepth: 1

   tutorial_notebook

----

Known Limitations of the Tutorial Data
---------------------------------------

The tutorial data are shallow by design. Expect the following, and don't read them as problems with your installation:

* **Two libraries fail QC.** ENCFF476MLA and ENCFF683BZL have very few variants, short phase blocks and almost no phased reads — consistent with failed or shallow runs. The notebook excludes them from the allelic analyses.
* **The abTU test is miscalibrated** on this subset (section 8), which inflates abTU calls and the abTU-containing classes.
* **The** *Igf2r* **imprinting control is inconclusive**, because it has only a handful of phased reads per sample.
* **Mild reference bias** towards B6 is visible after orientation (section 4).

----

Next Steps
----------

* **Inspect reads in IGV.** Load the haplotagged BAMs (``identify_loci/whatshap/<sample>/<sample>.flnc.hp.bam``) and group or colour reads by the ``HP`` tag for a haplotype-split view of any locus.
* **Scale up.** On full-depth libraries the imprinting controls should pass in each library on its own, and many more loci clear ``min_phased_reads_per_locus``.
* **Change theta.** Rerun with ``theta_mode = 'lift'`` or ``'fixed'`` and compare the QQ plot in section 8. With fewer than 3 samples, pooled calibration isn't recommended — see :doc:`nextflow_workflow_details`.
* **Run your own data.** Copy ``tutorial.config`` as a starting point, then change ``platform``, inputs, references, ``species`` and resources — see :doc:`nextflow_workflow_details` and the example configs.
