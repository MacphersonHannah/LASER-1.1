nextflow.enable.dsl = 2
// Speed up IsoQuant by parallelising

//  Split BAMs by chromosome
process SPLIT_BAM_BY_CHROM {
    tag "${sample} - ${chrom}"
    conda "${projectDir}/envs_prebuilt/isoquant_v2"
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

//  Run IsoQuant per chromosome
// Captures both reference-based and discovery outputs:
//   - transcript_counts.tsv: known transcript counts
//   - discovered_transcript_counts.tsv : novel + known transcript model counts
//   - read_assignments.tsv.gz: read > isoform mappings (used downstream for bulk_pbid lookup)
//   - transcript_model_reads.tsv.gz: reads > novel models only (may be empty)

process ISOQUANT_CHROM {
    tag "IsoQuant - ${chrom} (${labels.size()} samples)"
    
    conda "${projectDir}/envs_prebuilt/isoquant_v2"
    
    input:
    tuple val(chrom),
          val(labels),
          path(bams),
          path(bais)
    path  reference
    path  annotation
    val   stranded
    val   data_type
    
    output:
    tuple val(chrom), path("bulk_${chrom}.transcript_models.gtf"),              emit: gtf
    tuple val(chrom), path("bulk_${chrom}.transcript_counts.tsv"),              emit: counts
    tuple val(chrom), path("bulk_${chrom}.gene_counts.tsv"),                    emit: gene_counts
    tuple val(chrom), path("bulk_${chrom}.model_reads.tsv.gz"),                 emit: model_reads
    tuple val(chrom), path("bulk_${chrom}.sqanti_classification.tsv"),          emit: sqanti_classification
    tuple val(chrom), path("bulk_${chrom}.read_assignments.tsv.gz"),            emit: read_assignments
    tuple val(chrom), path("bulk_${chrom}.discovered_transcript_counts.tsv"),   emit: discovered_counts
    
    script:
    def bam_str   = (bams   instanceof List ? bams   : [bams]).collect { it.toString() }.join(' ')
    def label_str = (labels instanceof List ? labels : [labels]).collect { it.toString() }.join(' ')
    
    """
    # Stage references to local disk

    LOCAL_DIR=""
    trap '[ -n "\$LOCAL_DIR" ] && rm -rf "\$LOCAL_DIR"' EXIT
    for dir in /tmp /scratch /local; do
        if [ -d "\$dir" ] && [ -w "\$dir" ]; then
            LOCAL_DIR="\$dir/isoquant_refs_\${HOSTNAME}_\$\$"
            mkdir -p "\$LOCAL_DIR"
            break
        fi
    done
    
    if [ -n "\$LOCAL_DIR" ]; then
        echo "Staging references to local disk: \$LOCAL_DIR"
        cp ${reference} "\$LOCAL_DIR/" 2>/dev/null || true
        LOCAL_REF="\$LOCAL_DIR/${reference.name}"
        [ -f "\$LOCAL_REF" ] && REF="\$LOCAL_REF" || REF="${reference}"
        ANN="${annotation}"
    else
        echo "No local disk available, using NFS paths"
        REF="${reference}"
        ANN="${annotation}"
    fi
    
    # Run IsoQuant

    for bam in ${bam_str}; do
        [ -f "\${bam}.bai" ] || samtools index -@ ${task.cpus} "\${bam}"
    done
    
    isoquant \
        --data_type                  ${data_type}        \
        --reference                  \$REF               \
        --genedb                     \$ANN               \
        --bam                        ${bam_str}          \
        --labels                     ${label_str}        \
        --stranded                   ${stranded}         \
        --output                     out_${chrom}        \
        --prefix                     bulk_${chrom}       \
        --threads                    ${task.cpus}        \
        --transcript_quantification  all                 \
        --gene_quantification        all                 \
        --matching_strategy          loose               \
        --model_construction_strategy all                \
        --normalization_method       usable_reads        \
        --large_output               read2transcripts    \
        --count_exons                                    \
        --sqanti_output                                  \
        --min_mapq                   0                   \
        --inconsistent_mapq_cutoff   0

    # Clean local staging
    [ -n "\$LOCAL_DIR" ] && rm -rf "\$LOCAL_DIR"
    
    # Copy outputs - all have fallbacks so .collect() never hangs

    # GTF
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_models.gtf ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_models.gtf bulk_${chrom}.transcript_models.gtf
    else
        touch bulk_${chrom}.transcript_models.gtf
    fi

    # Reference transcript counts
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_counts.tsv ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_counts.tsv bulk_${chrom}.transcript_counts.tsv
    else
        touch bulk_${chrom}.transcript_counts.tsv
    fi

    # Gene counts
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.gene_counts.tsv ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.gene_counts.tsv bulk_${chrom}.gene_counts.tsv
    else
        touch bulk_${chrom}.gene_counts.tsv
    fi

    # Model reads (read > transcript model, known + novel via --large_output read2transcripts)
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_model_reads.tsv.gz ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_model_reads.tsv.gz bulk_${chrom}.model_reads.tsv.gz
    else
        printf '' | gzip > bulk_${chrom}.model_reads.tsv.gz
    fi

    # SQANTI classification
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_models_SQANTI.tsv ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.transcript_models_SQANTI.tsv bulk_${chrom}.sqanti_classification.tsv
    else
        touch bulk_${chrom}.sqanti_classification.tsv
    fi

    # Read assignments (read > known isoform)
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.read_assignments.tsv.gz ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.read_assignments.tsv.gz bulk_${chrom}.read_assignments.tsv.gz
    else
        printf '' | gzip > bulk_${chrom}.read_assignments.tsv.gz
    fi

    # Discovered transcript counts (novel + known)
    if [ -f out_${chrom}/bulk_${chrom}/bulk_${chrom}.discovered_transcript_counts.tsv ]; then
        cp out_${chrom}/bulk_${chrom}/bulk_${chrom}.discovered_transcript_counts.tsv \
           bulk_${chrom}.discovered_transcript_counts.tsv
    else
        touch bulk_${chrom}.discovered_transcript_counts.tsv
    fi
    """
    
    stub:
    """
    touch bulk_${chrom}.transcript_models.gtf
    touch bulk_${chrom}.transcript_counts.tsv
    touch bulk_${chrom}.gene_counts.tsv
    printf '' | gzip > bulk_${chrom}.model_reads.tsv.gz
    touch bulk_${chrom}.sqanti_classification.tsv
    printf '' | gzip > bulk_${chrom}.read_assignments.tsv.gz
    touch bulk_${chrom}.discovered_transcript_counts.tsv
    """
}

// Merge chromosome outputs
process MERGE_ISOQUANT {
    tag "Merge IsoQuant outputs"
    
    publishDir "${params.outdir}/preprocessing/isoquant", mode: 'copy'
    
    conda "${projectDir}/envs_prebuilt/isoquant_v2"
    
    input:
    path gtfs
    path counts
    path gene_counts
    path model_reads
    path sqanti_classifications
    path read_assignments
    path discovered_counts
    
    output:
    path "out/bulk/bulk.transcript_models.gtf",               emit: transcript_models_gtf
    path "out/bulk/bulk.transcript_counts.tsv",               emit: combined_counts
    path "out/bulk/bulk.gene_counts.tsv",                     emit: combined_gene_counts
    path "out/bulk/bulk.transcript_model_reads.tsv.gz",       emit: model_reads
    path "out/bulk/bulk.sqanti_classification.tsv",           emit: sqanti_classification
    path "out/bulk/bulk.read_assignments.tsv.gz",             emit: read_assignments
    path "out/bulk/bulk.discovered_transcript_counts.tsv",    emit: discovered_counts
    
    script:
    """
    mkdir -p out/bulk

    #  GTF 
    first_gtf=\$(ls *.transcript_models.gtf | head -1)
    grep "^#" "\$first_gtf" > out/bulk/bulk.transcript_models.gtf || true
    cat *.transcript_models.gtf | grep -v "^#" | sort -k1,1 -k4,4n \
        >> out/bulk/bulk.transcript_models.gtf

    #  reference transcript counts 
    first_counts=\$(ls *.transcript_counts.tsv | head -1)
    head -1 "\$first_counts" > out/bulk/bulk.transcript_counts.tsv || true
    cat *.transcript_counts.tsv | grep -v "^#\\|feature_id" \
        >> out/bulk/bulk.transcript_counts.tsv

    #  gene counts 
    first_gene=\$(ls *.gene_counts.tsv | head -1)
    head -1 "\$first_gene" > out/bulk/bulk.gene_counts.tsv || true
    cat *.gene_counts.tsv | grep -v "^#\\|feature_id" \
        >> out/bulk/bulk.gene_counts.tsv

    #  model reads (may be empty) 
    header_written=false
    for f in *.model_reads.tsv.gz; do
        first_line=\$(zcat "\$f" 2>/dev/null | head -1)
        if [ -n "\$first_line" ]; then
            echo "\$first_line" | gzip > header.tsv.gz
            header_written=true
            break
        fi
    done
    [ "\$header_written" = false ] && printf '' | gzip > header.tsv.gz
    printf '' | gzip > data.tsv.gz
    for f in *.model_reads.tsv.gz; do
        content=\$(zcat "\$f" 2>/dev/null | tail -n +2)
        [ -n "\$content" ] && echo "\$content" | gzip >> data.tsv.gz
    done
    cat header.tsv.gz data.tsv.gz > out/bulk/bulk.transcript_model_reads.tsv.gz
    rm header.tsv.gz data.tsv.gz

    #  SQANTI classification 
    first_sq=\$(grep -rl "." *.sqanti_classification.tsv 2>/dev/null | head -1)
    if [ -n "\$first_sq" ]; then
        head -1 "\$first_sq" > out/bulk/bulk.sqanti_classification.tsv
        cat *.sqanti_classification.tsv | grep -v "^#\\|isoform" \
            >> out/bulk/bulk.sqanti_classification.tsv
    else
        touch out/bulk/bulk.sqanti_classification.tsv
    fi

    #  read assignments (read > isoform, for bulk_pbid lookup) 
    # Merge all per-chromosome read_assignments.tsv.gz files.
    # Header is taken from first non-empty file; comment lines (#) are skipped.
    ra_header_written=false
    for f in *.read_assignments.tsv.gz; do
        [ -f "\$f" ] || continue
        first_line=\$(zcat "\$f" 2>/dev/null | grep -v "^#" | head -1)
        if [ -n "\$first_line" ]; then
            echo "\$first_line" | gzip > ra_header.tsv.gz
            ra_header_written=true
            break
        fi
    done
    [ "\$ra_header_written" = false ] && printf '' | gzip > ra_header.tsv.gz

    printf '' | gzip > ra_data.tsv.gz
    for f in *.read_assignments.tsv.gz; do
        [ -f "\$f" ] || continue
        content=\$(zcat "\$f" 2>/dev/null | grep -v "^#" | tail -n +2)
        [ -n "\$content" ] && echo "\$content" | gzip >> ra_data.tsv.gz
    done
    cat ra_header.tsv.gz ra_data.tsv.gz > out/bulk/bulk.read_assignments.tsv.gz
    rm ra_header.tsv.gz ra_data.tsv.gz

    # discovered transcript counts (novel + known transcript models)
    first_disc=\$(ls *.discovered_transcript_counts.tsv 2>/dev/null | head -1)
    if [ -n "\$first_disc" ]; then
        head -1 "\$first_disc" > out/bulk/bulk.discovered_transcript_counts.tsv
        cat *.discovered_transcript_counts.tsv | grep -v "^#\\|feature_id" \
            >> out/bulk/bulk.discovered_transcript_counts.tsv
    else
        touch out/bulk/bulk.discovered_transcript_counts.tsv
    fi
    """
    
    stub:
    """
    mkdir -p out/bulk
    touch out/bulk/bulk.transcript_models.gtf
    touch out/bulk/bulk.transcript_counts.tsv
    touch out/bulk/bulk.gene_counts.tsv
    printf '' | gzip > out/bulk/bulk.transcript_model_reads.tsv.gz
    touch out/bulk/bulk.sqanti_classification.tsv
    printf '' | gzip > out/bulk/bulk.read_assignments.tsv.gz
    touch out/bulk/bulk.discovered_transcript_counts.tsv
    """
}
