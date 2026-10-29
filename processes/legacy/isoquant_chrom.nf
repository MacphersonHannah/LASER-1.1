// processes/identify_loci/ont_specific/isoquant_chrom.nf
nextflow.enable.dsl = 2

// Subset each BAM to one chromosome 
process SPLIT_BAM_BY_CHROM {
    tag "${sample} - ${chrom}"

    conda "${projectDir}/envs_prebuilt/isoquant"

    input:
    tuple val(sample), path(bam), path(bai)
    each  chrom

    output:
    tuple val(chrom), val(sample),
          path("${sample}.${chrom}.bam"),
          path("${sample}.${chrom}.bam.bai")

    script:
    """
    samtools view -b ${bam} ${chrom} > ${sample}.${chrom}.bam
    samtools index ${sample}.${chrom}.bam
    """
}

// IsoQuant on one chromosome across all samples
// For targeted libraries, many chromosomes will have zero or very few reads.
// isoQuant exits 0 but silently skips writing some output files
// All cp commands are therefore guarded

process ISOQUANT_CHROM {
    tag "IsoQuant - ${chrom} (${labels.size()} samples)"


    conda "${projectDir}/envs_prebuilt/isoquant"

    input:
    tuple val(chrom),
          val(labels),
          path(bams),
          path(bais)
    path  reference
    path  annotation
    val   stranded      // reverse (cDNA ONT) or forward (direct RNA / PacBio cDNA)
    val   data_type     // nanopore or pacbio

    output:
    tuple val(chrom), path("bulk_${chrom}.transcript_models.gtf"),         emit: gtf
    tuple val(chrom), path("bulk_${chrom}.transcript_counts.tsv"),         emit: counts
    tuple val(chrom), path("bulk_${chrom}.gene_counts.tsv"),               emit: gene_counts
    tuple val(chrom), path("bulk_${chrom}.model_reads.tsv.gz"),            emit: model_reads
    tuple val(chrom), path("bulk_${chrom}.sqanti_classification.tsv"),     emit: sqanti_classification
    tuple val(chrom), path("bulk_${chrom}.read_assignments.tsv.gz"),       emit: read_assignments, optional: true

    script:
    def bam_str   = (bams   instanceof List ? bams   : [bams]).collect { it.toString() }.join(' ')
    def label_str = (labels instanceof List ? labels : [labels]).collect { it.toString() }.join(' ')

    """
    for bam in ${bam_str}; do
        [ -f "\${bam}.bai" ] || samtools index -@ ${task.cpus} "\${bam}"
    done

    isoquant \\
        --data_type                  ${data_type}        \\
        --reference                  ${reference}        \\
        --genedb                     ${annotation}       \\
        --complete_genedb                                \\
        --bam                        ${bam_str}          \\
        --labels                     ${label_str}        \\
        --stranded                   ${stranded}         \\
        --output                     out_${chrom}        \\
        --prefix                     bulk_${chrom}       \\
        --threads                    ${task.cpus}        \\
        --transcript_quantification  all                 \\
        --gene_quantification        all                 \\
        --matching_strategy          loose               \\
        --normalization_method       usable_reads        \\
        --count_exons                                    \\
        --sqanti_output                                  \\
        --min_mapq                   0                   \\
        --inconsistent_mapq_cutoff   0

    # GTF
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_models.gtf ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_models.gtf \
           bulk_${chrom}.transcript_models.gtf
    else
        touch bulk_${chrom}.transcript_models.gtf
    fi

    # transcript counts 
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_counts.tsv ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_counts.tsv \
           bulk_${chrom}.transcript_counts.tsv
    else
        touch bulk_${chrom}.transcript_counts.tsv
    fi

    # gene counts 
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.gene_counts.tsv ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.gene_counts.tsv \
           bulk_${chrom}.gene_counts.tsv
    else
        touch bulk_${chrom}.gene_counts.tsv
    fi

    # model reads (must be valid gzip - zcat is used downstream)
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_model_reads.tsv.gz ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_model_reads.tsv.gz \
           bulk_${chrom}.model_reads.tsv.gz
    else
        printf '' | gzip > bulk_${chrom}.model_reads.tsv.gz
    fi

    # SQANTI classification ─
    # Most commonly absent for chromosomes with few reads in targeted libraries
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_models_SQANTI.tsv ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_models_SQANTI.tsv \
           bulk_${chrom}.sqanti_classification.tsv
    else
        touch bulk_${chrom}.sqanti_classification.tsv
    fi

    # read assignments (already optional, keep existing guard) 
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.read_assignments.tsv.gz ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.read_assignments.tsv.gz \
           bulk_${chrom}.read_assignments.tsv.gz
    fi
    """

    stub:
    """
    touch bulk_${chrom}.transcript_models.gtf
    touch bulk_${chrom}.transcript_counts.tsv
    touch bulk_${chrom}.gene_counts.tsv
    printf '' | gzip > bulk_${chrom}.model_reads.tsv.gz
    touch bulk_${chrom}.sqanti_classification.tsv
    touch bulk_${chrom}.read_assignments.tsv.gz
    """
}


// Merge per-chromosome outputs
// The header is taken from the first non empty file where possible.

process MERGE_ISOQUANT {
    tag "Merge IsoQuant outputs"

    publishDir "${params.outdir}/isoquant", mode: 'copy'

    conda "${projectDir}/envs_prebuilt/isoquant"

    input:
    path gtfs
    path counts
    path gene_counts
    path model_reads
    path sqanti_classifications

    output:
    path "out/bulk/bulk.transcript_models.gtf",            emit: transcript_models_gtf
    path "out/bulk/bulk.transcript_counts.tsv",            emit: combined_counts
    path "out/bulk/bulk.gene_counts.tsv",                  emit: combined_gene_counts
    path "out/bulk/bulk.transcript_model_reads.tsv.gz",    emit: model_reads
    path "out/bulk/bulk.sqanti_classification.tsv",        emit: sqanti_classification

    script:
    """
    mkdir -p out/bulk

    # GTF 
    # Take header from first non-empty GTF; concat all data lines
    first_gtf=\$(cat <(ls *.transcript_models.gtf) | xargs grep -l "." 2>/dev/null | head -1 || ls *.transcript_models.gtf | head -1)
    grep "^#" "\$first_gtf" > out/bulk/bulk.transcript_models.gtf || true
    cat *.transcript_models.gtf | grep -v "^#" | sort -k1,1 -k4,4n \
        >> out/bulk/bulk.transcript_models.gtf

    # transcript counts 
    first_counts=\$(grep -l "." *.transcript_counts.tsv 2>/dev/null | head -1 || ls *.transcript_counts.tsv | head -1)
    head -1 "\$first_counts" > out/bulk/bulk.transcript_counts.tsv || true
    cat *.transcript_counts.tsv | grep -v "^#\\|feature_id" \
        >> out/bulk/bulk.transcript_counts.tsv

    # gene counts 
    first_gene=\$(grep -l "." *.gene_counts.tsv 2>/dev/null | head -1 || ls *.gene_counts.tsv | head -1)
    head -1 "\$first_gene" > out/bulk/bulk.gene_counts.tsv || true
    cat *.gene_counts.tsv | grep -v "^#\\|feature_id" \
        >> out/bulk/bulk.gene_counts.tsv

    # model reads 
    # Find first non-empty gzip to extract header; skip empty gzips silently
    header_written=false
    for f in *.model_reads.tsv.gz; do
        first_line=\$(zcat "\$f" 2>/dev/null | head -1)
        if [ -n "\$first_line" ]; then
            echo "\$first_line" | gzip > header.tsv.gz
            header_written=true
            break
        fi
    done
    if [ "\$header_written" = false ]; then
        printf '' | gzip > header.tsv.gz
    fi

    printf '' | gzip > data.tsv.gz
    for f in *.model_reads.tsv.gz; do
        content=\$(zcat "\$f" 2>/dev/null | tail -n +2)
        if [ -n "\$content" ]; then
            echo "\$content" | gzip >> data.tsv.gz
        fi
    done
    cat header.tsv.gz data.tsv.gz > out/bulk/bulk.transcript_model_reads.tsv.gz
    rm header.tsv.gz data.tsv.gz

    # SQANTI 
    first_sq=\$(grep -l "." *.sqanti_classification.tsv 2>/dev/null | head -1 || ls *.sqanti_classification.tsv | head -1)
    head -1 "\$first_sq" > out/bulk/bulk.sqanti_classification.tsv || true
    cat *.sqanti_classification.tsv | grep -v "^#\\|isoform" \
        >> out/bulk/bulk.sqanti_classification.tsv
    """

    stub:
    """
    mkdir -p out/bulk
    touch out/bulk/bulk.transcript_models.gtf
    touch out/bulk/bulk.transcript_counts.tsv
    touch out/bulk/bulk.gene_counts.tsv
    printf '' | gzip > out/bulk/bulk.transcript_model_reads.tsv.gz
    touch out/bulk/bulk.sqanti_classification.tsv
    """
}
