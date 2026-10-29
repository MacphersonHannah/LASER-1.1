.. _troubleshooting:

Troubleshooting and Notes
==========================

This page covers common errors and platform-specific notes for running LASER.

----

Troubleshooting
---------------

``Variable declarations cannot be mixed with config statements``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run ``bash scripts/setup/validate_config.sh {config_file} --fix`` — this catches the usual cause (a stray byte-order mark or Windows line endings left behind after a config file was re-saved by another editor). See `nextflow-io/nextflow#5997 <https://github.com/nextflow-io/nextflow/issues/5997>`_ for background.

``No space left on device`` from Singularity/Apptainer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Singularity stages container images to local ``/tmp`` during a run. This error usually just means ``/tmp`` on that node is full. Check available space and either free some up or run on a node with more local disk.

WhatsHap haplotagging gives very low phasing rates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ensure ``--no-reference`` is set in the haplotag command in ``whatshap_phase_tag.nf``. Without it, WhatsHap's realignment-based allele detection fails on reads with multiple indels near phase anchors, which is common in PacBio data. The phasing step should still use ``--reference`` for best phasing quality. This can also be a biological issue — loci must have expressed SNV or indel variants that differ between alleles to be phased.

Pipeline resumes from cache
~~~~~~~~~~~~~~~~~~~~~~~~~~~

LASER supports Nextflow's ``-resume`` flag. All upstream steps that haven't changed will be cached and reused. ``run_LASER.sh`` passes ``-resume`` automatically; to force a fully fresh run, remove ``-resume`` from the ``NF_CMD`` array inside ``run_LASER.sh`` or delete the original output folder.

Testing against an individual's own genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only the human/mouse genome has been used as a reference in testing, but testing RNA/cDNA against the individual's own genome should theoretically work. Feel free to get in touch if you're testing this and have any trouble.

Running single samples
~~~~~~~~~~~~~~~~~~~~~~

Whilst LASER is equipped to do this, you cannot derive a theta value for the statistical tests from only one sample, as it requires multiple samples to get a measure of noisiness. You can, however, import a previously generated theta (``theta_mode = 'fixed'`` or ``'lift'`` — see :doc:`nextflow_workflow_details`). You will likely need to run multiple samples together initially to generate this number before running your samples of interest separately. We haven't fully investigated what causes the greatest theta changes, so the guidance is to generate a theta from the libraries most similar to your libraries of interest — i.e. same platform, same tissue type, ideally same library prep.

----

.. _notes_mas_isoseq:

Notes on MAS-Iso-Seq
---------------------

MAS-Iso-Seq libraries must be processed through PacBio's Skera tool before passing to LASER, followed by demultiplexing with Lima. Pass the resulting BAMs via ``params.library_dirs`` as described in :doc:`nextflow_workflow_details`. Use ``params.dv_model = 'MASSEQ'`` for variant calling.

----

.. _notes_non_human:

Notes on Non-Human Data
-----------------------

LASER is species-agnostic at the pipeline level. The following points require attention when running on non-human data.

``params.species`` must be set (required for non-human data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both ORFannotate (coding-potential model) and SQANTI3 (whether to apply the human-specific ``--coverage``/intropolis junction-support flag) read this single param. The default is ``'human'``. Set this to match your organism — for mouse:

.. code-block:: groovy

   params.species = 'mouse'

Forgetting this will score transcripts against the wrong hexamer model, and silently skip SQANTI3's junction-support annotation for human runs where it was never set.

SQANTI3 Reference Files Need Species-Matched Equivalents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``cage``, ``intropolis``, ``polya_atlas``, and ``poly_motifs`` files listed in :doc:`key_outputs` are human-specific. For mouse:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Mouse equivalent
   * - ``cage``
     - refTSS_v3.3 mouse coordinate BED (download from the refTSS database)
   * - ``polya_atlas``
     - Mouse PolyA_DB or APADB atlas BED
   * - ``poly_motifs``
     - Human list is appropriate — polyA hexamers are conserved
   * - ``intropolis``
     - No mouse equivalent exists. Pass an empty file: ``touch empty_intropolis.tsv``

Passing an empty file for ``intropolis`` is safe — SQANTI3 will run without junction support evidence and simply leave that column unpopulated.

DeepVariant Model
~~~~~~~~~~~~~~~~~

``dv_model = 'MASSEQ'`` is correct for mouse PacBio Iso-Seq data. The MASSEQ model is trained on PacBio read characteristics, not human-specific sequence content.
