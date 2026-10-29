.. LASER documentation master file

LASER: Long-read Allele-Specific Expression Resolver
======================================================

.. figure:: ../img/LASER-logo.png
   :alt: LASER Logo
   :align: center
   :width: 60%

LASER is a Nextflow pipeline for long-read RNA-seq analysis, designed to identify allele-specific expression (ASE) and transcript usage from PacBio Iso-Seq and Oxford Nanopore Technology (ONT) data.

----

Overview
--------

Allele-specific expression (ASE) — where the two copies of a gene are expressed at different levels — is a key mechanism in disease, dosage compensation, and imprinting. Long-read sequencing offers a direct route to resolving ASE, because individual reads are long enough to span both heterozygous SNVs and splice junctions, linking variant identity to transcript structure in a single molecule.

LASER is built to exploit this. Given long-read RNA-seq data from either PacBio Iso-Seq or ONT, it performs end-to-end processing: read refinement and alignment, transcript-level isoform collapse and quantification, variant calling and phasing, and finally statistical assessment of allelic bias — both in total expression (abTE) and in transcript usage (abTU).

LASER is designed for diploid organisms, has been tested on human and mouse tissue and iPSC samples, and is species-agnostic at the pipeline level.

New to LASER? The :doc:`tutorial` runs the full pipeline on a small mouse dataset and walks through the results.


I. The LASER Nextflow Pipeline
-------------------------------

LASER is a scalable, automated Nextflow pipeline that takes raw long-read RNA-seq data from demultiplexing to per-locus ASE calls.

**Key features:**

* **Multi-platform support**: Compatible with PacBio Iso-Seq (Sequel IIe, Revio) and ONT (cDNA and direct RNA) data, with multiple entry points depending on upstream processing — demultiplexed BAMs or refined FLNC reads for PacBio, FASTQ or unmapped BAMs for ONT.

* **Isoform-aware quantification**: Uses isoseq3 (PacBio) or FLAIR (ONT) to resolve and quantify full-length isoforms, with SQANTI3 QC and ORFannotate annotation.

* **Integrated variant calling and phasing**: Calls variants with DeepVariant (PacBio default; Clair3-RNA optional) or Clair3-RNA (ONT), phases with WhatsHap, and haplotags reads for haplotype-resolved expression.

* **ASE and abTU assessment**: Identifies loci with allelic bias in total expression and transcript usage, with configurable thresholds for read depth, phasing quality, and isoform diversity.

Pipeline Overview
~~~~~~~~~~~~~~~~~

.. mermaid::
   :alt: LASER pipeline flowchart, from input reads through preprocessing, locus identification, isoform annotation and transcript summary extraction to allelic effect testing

   flowchart TD
       A["Input reads<br/>(PacBio FLNC / ONT FASTQ / ONT BAM)"] --> B["Platform-specific preprocessing<br/>mapping, isoform discovery/collapse,<br/>variant calling"]
       B --> C["identify_loci<br/>read-sample map,<br/>WhatsHap phasing + haplotagging"]
       C --> D["annotate_isoforms<br/>SQANTI3 QC, ORF annotation,<br/>haplotype tag extraction,<br/>read-to-isoform lookup"]
       D --> E["Transcript summary extraction<br/>read x allele x transcript counts"]
       E --> F{"Allelic effect testing"}
       F --> G["theta_mode:<br/>calibrate_pooled / calibrate_pooled_by_group / fixed / lift"]
       G --> H["abTE test<br/>binomial, allele totals"]
       G --> I["abTU test<br/>permutation Dirichlet-Multinomial"]
       H --> J["Four-way classification<br/>per locus, per sample"]
       I --> J
       J --> K["allelic_effect_testing/&lt;sample&gt;/<br/>allelic_effect_results.tsv"]



II. Supported Input Types
--------------------------

**PacBio (Iso-Seq)**

* Any PacBio machine (Sequel IIe, Revio)
* Single-sample or multiplexed libraries
* MAS-Iso-Seq or classic Iso-Seq
* Targeted or transcriptome-wide

**ONT (Oxford Nanopore)**

* Direct RNA sequencing (``ont_directrna`` platform)
* cDNA sequencing (``ont_cdna`` platform) — from FASTQ, includes pychopper adapter trimming
* cDNA sequencing (``ont_cdna_minimap_entry`` platform) — from unmapped BAMs, skips pychopper
* One FASTQ file per sample (named ``SAMPLE.fastq.gz``) for ``ont_cdna`` and ``ont_directrna``; one unmapped BAM per sample for ``ont_cdna_minimap_entry``


III. Sample and Locus Eligibility
----------------------------------

**Sample:** LASER is designed for diploid organisms and has been tested on human and mouse tissue and iPSC samples.

**Loci:** For a locus to be eligible it must be diploid — targets on the Y chromosome and mitochondrial genome are excluded. Loci must contain informative positions: heterozygous, expressed SNVs that differ between the two alleles, typically exonic or in UTRs. LASER identifies these automatically.

Assessment of allelic bias in total expression (abTE) and allelic bias in transcript usage (abTU) is performed by the allelic effect testing stage, and requires total expression above a configurable threshold (default: 30 reads, ``min_total_reads_per_locus``), a minimum number of phased reads (default: 30, ``min_phased_reads_per_locus``), and a bound on the unphased proportion (default: 30%, ``max_unphased_prop``). Assessment of abTU additionally requires more than one detected transcript at the locus (default: 2, ``min_transcripts_per_locus``). See :doc:`nextflow_workflow_details` for the full parameter list.

----

.. toctree::
   :maxdepth: 1
   :caption: Pipeline Documentation:

   installation_and_requirements
   nextflow_workflow_details
   key_outputs
   troubleshooting

.. toctree::
   :maxdepth: 1
   :caption: Tutorial:

   tutorial
