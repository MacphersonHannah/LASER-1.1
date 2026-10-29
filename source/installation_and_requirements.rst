.. _installation_and_requirements:

Installation and Requirements
==============================

This page covers system requirements, how to install LASER's dependencies, and how to run the pipeline for the first time.

Prerequisites
-------------

* **Linux** (tested on HPC cluster nodes; not tested on macOS/Windows).
* **conda or mamba, via** `Miniforge <https://github.com/conda-forge/miniforge>`_. This is the recommended way to get both. Installing conda on its own (e.g. via Anaconda/Miniconda without mamba) will work but environment solves will be noticeably slower; LASER's setup scripts use ``mamba`` for this reason. If you already have Miniconda, install mamba into it (``conda install -n base -c conda-forge mamba``) rather than starting over.
* **git, wget and a system python3** — used by ``setup_tools.sh`` to fetch and unpack SQANTI3/ORFannotate.
* ``setup_envs.sh`` will install ``jq`` and ``curl`` into your conda **base** environment if they aren't already on PATH (a SQANTI3 post-link script needs them).
* **Singularity or Apptainer** (used for DeepVariant/Clair3-RNA containers). Check with ``singularity --version`` or ``apptainer --version`` / ``module load singularity`` on HPC systems. Ask your HPC administrator if neither is available.
* **Sufficient compute for setup, not just for runs** — see the warning below before running the setup scripts (Steps 2 & 3 below).

.. warning::
   **Do not run** ``setup_tools.sh`` **/** ``setup_envs.sh`` **on an HPC login node.** These scripts do real work (cloning/downloading tools, solving and building multiple conda environments) and commonly exceed login-node limits — many HPC systems cap login nodes at a handful of CPUs and a small amount of memory per user, and setup will simply be killed partway through. Submit setup as a batch job or run it on an interactive compute node instead. A generic SLURM example (adjust partition/account/time/resources to your cluster):

   .. code-block:: bash

      sbatch --job-name=laser_setup \
             --cpus-per-task=4 \
             --mem=32G \
             --time=02:00:00 \
             --wrap="bash scripts/setup/setup_tools.sh && bash scripts/setup/setup_envs.sh"

   The same resource limits apply to the actual pipeline run, which is far more resource-hungry than setup — always run LASER itself on a compute node too, not a login node.

----

Quickstart
----------

This guide assumes a Linux system with conda installed.

**1. Clone the repository**

.. code-block:: bash

   git clone https://github.com/MacphersonHannah/LASER-1.1.git
   cd LASER-1.1

**2. Install required tools**

Run this from the LASER-1.1 root directory (the one containing ``main.nf``); both setup scripts use paths relative to it.

.. code-block:: bash

   bash scripts/setup/setup_tools.sh

This checks for Singularity, downloads ORFannotate, SQANTI3, and DeepVariant, and creates a DeepVariant sandbox. Speak to your HPC administrator if Singularity is not available.

**3. Create required environments**

Again from the LASER-1.1 root directory.

.. code-block:: bash

   bash scripts/setup/setup_envs.sh

Creates all conda environments under ``envs_prebuilt/``. Existing environments are skipped. If you need to remake an environment, delete its folder inside ``envs_prebuilt/`` before rerunning this script.

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Env
     - Covers
   * - ``pacbio_prep_env``
     - Mapping, refinement (FLNC), demultiplexing (pbmm2, minimap2, isoseq, lima)
   * - ``core_utils_env``
     - Orchestration + phasing/VCF utilities (nextflow, whatshap, bcftools). For manual use: ``conda activate envs_prebuilt/core_utils_env``
   * - ``isoform_collapse_env``
     - Isoform collapse + basecall trimming (cdna_cupcake, pychopper)
   * - ``stats_r_env``
     - Allelic effect testing statistics (R)
   * - ``flair_env``
     - Tool-specific, kept separate due to incompatible dependency pins
   * - ``SQANTI3``, ``ORFannotate``
     - Tool-specific, installed by ``setup_tools.sh``

**4. Configure your run**

Copy the example config that matches your platform from ``example_configs/`` and edit the paths to match your data and references. See :doc:`nextflow_workflow_details` and the comments inside each example for guidance.

.. code-block:: bash

   cp example_configs/example_ont_cdna_minimap_entry.config my_run.config
   # edit my_run.config with your paths, then pass it to run_LASER.sh (step 5)

See ``example_configs/``: the :ref:`example configs table <example_configs>` lists which file matches each platform.

**5. Run LASER**

.. code-block:: bash

   bash run_LASER.sh {RUN_NAME} {config_file}

This activates the required environment (``core_utils_env``) automatically.

Working files are written to ``run_{RUN_NAME}/work_{RUN_NAME}/``. Results are written to ``run_{RUN_NAME}/{RUN_NAME}/``. The config file used is automatically copied into the results directory for reproducibility, along with the Nextflow run report, timeline and trace.

Results are hard-linked to the work directory rather than copied, so they cost almost no extra disk and survive deletion of the work directory.

``RUN_NAME`` decides where output goes — ``run_LASER.sh`` always passes ``--outdir`` on the Nextflow command line. Edit ``RUN_NAME`` when starting a new run. To re-run, call it the same ``RUN_NAME`` and the ``-resume`` flag in ``run_LASER.sh`` will resume from whatever has changed. To bypass this, delete ``-resume`` in ``run_LASER.sh`` or delete the folder of the original run output.

.. note::
   Run this on a compute node (or submit as a job), not a login node — see the warning under Requirements above; the pipeline run is far more resource-hungry than setup.

**6. (Optional) Run on test data first**

Each platform has a small test dataset bundled under ``test_data/``, with a ready-to-use config file. These are subsampled to 50,000 reads (10,000 for ``pacbio``) per sample and should complete in under 30 minutes on a node with 16+ cores and 120 GB RAM. Running a test first is a good way to check that setup completed correctly and your cluster environment works before committing to a full-sized run.

.. code-block:: bash

   # Example: test the ont_cdna_minimap_entry platform
   bash run_LASER.sh test_minimap_entry test_data/ONT_cDNA_minimap_input/ONT_cDNA_test_data_minimap_input.config

.. list-table::
   :widths: 25 40 15 20
   :header-rows: 1

   * - Platform
     - Test data directory
     - Test config
     - Samples
   * - ``pacbio``
     - ``test_data/PACBIO_untargeted_not_fastq_entry/``
     - Included
     - 2 samples, 10k reads each
   * - ``pacbio_fastqentry``
     - ``test_data/PACBIO_untargeted_fastq_entry/``
     - Included
     - 2 ENCODE samples, 50k reads each
   * - ``ont_cdna``
     - ``test_data/ONT_cDNA/``
     - Included
     - 2 control brain samples, 50k reads each
   * - ``ont_cdna_minimap_entry``
     - ``test_data/ONT_cDNA_minimap_input/``
     - Included
     - 2 samples, 50k reads each
   * - ``ont_directrna``
     - ``test_data/ONT_directRNA/``
     - Included
     - 3 ENCODE samples, 50k reads each

You will need to update the reference file paths in the test config to match your local installation before running — the test configs use placeholder paths.

----

Platform Selection
------------------

Set ``params.platform`` in your config to one of:

.. list-table::
   :widths: 16 40 14 14 16
   :header-rows: 1

   * - Platform
     - Description
     - Isoform discovery
     - Variant calling
     - Test data
   * - ``pacbio``
     - PacBio Iso-Seq from raw demultiplexed BAMs. Runs the full IsoSeq3 pipeline: refine, cluster2, collapse, and abundance. Standard or MAS-Iso-Seq (run Skera first for MAS).
     - isoseq3 cluster2 + collapse
     - DeepVariant (default) or Clair3-RNA
     - ``test_data/PACBIO_untargeted_not_fastq_entry/``
   * - ``pacbio_fastqentry``
     - PacBio Iso-Seq from already-refined FLNC reads (BAM or fastq.gz). Enters at cluster2, skipping demux and refine. Use this for data downloaded from repositories, processed on another system, or MAS-Iso-Seq that has already been through Skera + lima + isoseq refine.
     - isoseq3 cluster2 + collapse
     - DeepVariant (default) or Clair3-RNA
     - ``test_data/PACBIO_untargeted_fastq_entry/``
   * - ``ont_cdna``
     - ONT cDNA sequencing from FASTQ. Includes pychopper for adapter trimming before alignment. Full-length reads go to FLAIR for isoform discovery; if ``pychopper_rescued = true``, full-length plus rescued reads are used for everything else (variant calling, phasing and quantification). Set ``params.ont_chemistry`` to ``r9`` or ``r10`` to select the correct minimap2 preset.
     - FLAIR
     - Clair3-RNA
     - ``test_data/ONT_cDNA/``
   * - ``ont_cdna_minimap_entry``
     - ONT cDNA sequencing from unmapped BAMs which have already been processed with Pychopper. Skips pychopper entirely; all reads go through a single FLAIR stream for both isoform discovery and quantification. Use this when your data has not been through pychopper, or when you want to skip adapter trimming.
     - FLAIR
     - Clair3-RNA
     - ``test_data/ONT_cDNA_minimap_input/``
   * - ``ont_directrna``
     - ONT direct RNA sequencing. No adapter trimming (pychopper is not used). Reads go straight to minimap2 alignment and FLAIR for isoform discovery.
     - FLAIR
     - Clair3-RNA
     - ``test_data/ONT_directRNA/``

----

Preparing Input Data
--------------------

PacBio: Demultiplexing
~~~~~~~~~~~~~~~~~~~~~~

LASER never runs ``lima``/demultiplexing itself — this is always a step you run yourself, before pointing LASER at your data. Your PacBio BAMs must already be demultiplexed into one BAM per sample by the time you configure a run. The required steps depend on your library type:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Library type
     - Steps before LASER
   * - Standard Iso-Seq
     - Run ``demux_pacbio.sh``
   * - MAS-Iso-Seq
     - Run Skera first, then run ``demux_pacbio.sh`` on the output

**Standard Iso-Seq:**

.. code-block:: bash

   bash scripts/setup/demux_pacbio.sh input.bam primers.fa output_dir

**MAS-Iso-Seq** — first run Skera to segment the concatenated reads, then demultiplex:

.. code-block:: bash

   # Step 1: Skera (must be run separately — not bundled with LASER)
   skera split input.bam mas_primers.fasta skera_output.bam

   # Step 2: Demultiplex with Lima
   bash scripts/setup/demux_pacbio.sh skera_output.bam primers.fa output_dir

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Description
   * - ``input.bam``
     - Raw CCS BAM from the PacBio sequencer (or Skera output for MAS-Iso-Seq)
   * - ``primers.fa``
     - FASTA file containing your 5' and 3' barcode sequences
   * - ``output_dir``
     - Directory to write demultiplexed BAMs — created if it doesn't exist

The primer FASTA is automatically copied into ``output_dir`` alongside the BAMs, as LASER requires them to be in the same directory. Once complete, point LASER at the output directory:

.. code-block:: groovy

   params.library_dirs = ['/path/to/output_dir']

For multiple libraries, provide one output directory per library — each must contain its own BAMs and primer file. LASER iterates over every entry in ``library_dirs`` automatically, so no extra flag is needed. See :doc:`nextflow_workflow_details` for details.

.. note::
   MAS-Iso-Seq libraries must be processed through PacBio's Skera tool before passing to LASER, followed by demultiplexing with Lima. Use ``params.dv_model = 'MASSEQ'`` for variant calling. See :ref:`notes_mas_isoseq` for details.

ONT: Concatenating FASTQ Files (``ont_cdna``, ``ont_directrna``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LASER expects one concatenated FASTQ per sample named ``SAMPLE.fastq.gz``, where the filename becomes the sample name throughout all outputs — so name files clearly (e.g. ``patient1_thalamus.fastq.gz``).

Most basecallers produce multiple FASTQ files per sample (e.g. one per batch in ``fastq_pass/``). Concatenate them using the provided helper script:

.. code-block:: bash

   bash scripts/setup/concat_ONT_fastq.sh fastq_pass/ merged_fastq/ patient1_thalamus

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Description
   * - ``fastq_pass/``
     - Directory containing your raw FASTQ files
   * - ``merged_fastq/``
     - Output directory — created if it doesn't exist
   * - ``patient1_thalamus``
     - Sample name — output will be ``merged_fastq/patient1_thalamus.fastq.gz``

Run once per sample. Once all samples are concatenated:

.. code-block:: groovy

   params.ont_input_files = [
       '/path/to/patient1_thalamus.fastq.gz',
       '/path/to/patient2_thalamus.fastq.gz'
   ]

ONT: Unmapped BAMs (``ont_cdna_minimap_entry``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No FASTQ concatenation is needed for this platform. Each unmapped BAM becomes one sample, and the filename (minus ``.bam``) becomes the sample label. See :doc:`nextflow_workflow_details` for the config.
