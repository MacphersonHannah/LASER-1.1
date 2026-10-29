// StringTie2 long-read transcript discovery and quantification for ONT cDNA.
//   STRINGTIE - per-sample transcript discovery (no -e)
//   STRINGTIE_MERGE - unified transcript catalogue across all samples
//   STRINGTIE_QUANT - per-sample quantification against merged catalogue (-e -b)
// Read-to-transcript info is produced by featureCounts in transcript mode.
nextflow.enable.dsl = 2

// Per-sample transcript discovery 
process STRINGTIE {
    tag { sample }
    publishDir "${params.outdir}/preprocessing/stringtie/discovery/${sample}", mode: 'link'
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    input:
    tuple val(sample), path(bam)
    path annotation
    val  strand_flag   // '--rf' for reverse-stranded cDNA (PCB114), '--fr' for direct RNA

    output:
    path "${sample}.stringtie.gtf", emit: gtf

    script:
    """
    stringtie \
        -L \
        -G ${annotation} \
        ${strand_flag} \
        -o ${sample}.stringtie.gtf \
        -p ${task.cpus} \
        ${bam}
    """
}

// Merge per-sample GTFs into unified catalogue
process STRINGTIE_MERGE {
    tag "StringTie merge"
    publishDir "${params.outdir}/preprocessing/stringtie/merged", mode: 'copy'
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    input:
    path gtfs           // all per-sample GTFs collected
    path annotation

    output:
    path "merged.gtf", emit: merged_gtf

    script:
    """
    ls *.stringtie.gtf > gtf_list.txt
    stringtie --merge \
        -G ${annotation} \
        -o merged.gtf \
        gtf_list.txt
    """
}

// sample quantification against merged catalogue
process STRINGTIE_QUANT {
    tag { sample }
    publishDir "${params.outdir}/preprocessing/stringtie/quant/${sample}", mode: 'copy'
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    input:
    tuple val(sample), path(bam)
    path merged_gtf
    val  strand_flag
    val  fc_strand     // featureCounts strand: 0=unstranded, 1=forward, 2=reverse

    output:
    tuple val(sample), path("${sample}.quant.gtf"),                    emit: quant_gtf
    tuple val(sample), path("ballgown/${sample}"),                      emit: ballgown
    tuple val(sample), path("${sample}.transcript_reads.tsv.gz"),       emit: transcript_reads

    script:
    """
    # StringTie quantification 
    mkdir -p ballgown/${sample}
    stringtie \
        -L \
        -e \
        -b ballgown/${sample} \
        -G ${merged_gtf} \
        ${strand_flag} \
        -o ${sample}.quant.gtf \
        -p ${task.cpus} \
        ${bam}

    # featureCounts: read > transcript assignment 
    # -f = feature level (exon), -O = allow overlap, -s = strandedness
    # -R BAM writes per-read assignment to a BAM file
    featureCounts \
        -a ${merged_gtf} \
        -o ${sample}.fc_counts.txt \
        -f \
        -O \
        -s ${fc_strand} \
        --primary \
        -R BAM \
        -T ${task.cpus} \
        ${bam}

    # Convert featureCounts BAM tags to read>transcript TSV 
    # XS tag = assigned feature (transcript_id from GTF)
    # XT tag = status: Assigned / Unassigned_*
    # Output format matches IsoQuant model_reads: read_id <tab> transcript_id
    samtools view ${bam}.featureCounts.bam | \
        awk 'BEGIN{OFS="\\t"; print "read_id","transcript_id"}
        {
            read_id=\$1; tx="*"
            for(i=12;i<=NF;i++){
                if(\$i~/^XS:Z:/) tx=substr(\$i,6)
                if(\$i~/^XT:Z:Assigned/) assigned=1
            }
            if(tx!="*" && assigned) print read_id, tx
            assigned=0
        }' | gzip > ${sample}.transcript_reads.tsv.gz

    # Clean up large intermediate BAM
    rm -f ${bam}.featureCounts.bam
    """
}

// Merge per-sample transcript_reads into bulk file 
process MERGE_STRINGTIE_READS {
    tag "Merge StringTie read assignments"
    publishDir "${params.outdir}/preprocessing/stringtie/merged", mode: 'copy'
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    input:
    path transcript_reads   // all per-sample transcript_reads.tsv.gz collected

    output:
    path "bulk.transcript_model_reads.tsv.gz", emit: model_reads

    script:
    """
    # Write header once, then append data from all samples
    echo -e "read_id\\ttranscript_id" | gzip > bulk.transcript_model_reads.tsv.gz
    for f in *.transcript_reads.tsv.gz; do
        zcat "\$f" | tail -n +2 | gzip >> bulk.transcript_model_reads.tsv.gz
    done
    """
}

// Build combined count matrix from StringTie ballgown tables
process STRINGTIE_COUNT_MATRIX {
    tag "StringTie count matrix"
    publishDir "${params.outdir}/preprocessing/stringtie/counts", mode: 'copy'
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    input:
    path quant_gtfs     // all per-sample quant GTFs collected

    output:
    path "transcript_counts.tsv", emit: counts

    script:
    """
    #!/usr/bin/env python3
    import os, re
    from collections import defaultdict

    gtf_files = sorted([f for f in os.listdir('.') if f.endswith('.quant.gtf')])
    samples = [re.sub(r'\\.quant\\.gtf$', '', f) for f in gtf_files]

    # Parse coverage (TPM proxy) from StringTie quant GTFs
    # StringTie outputs cov (coverage) and TPM in transcript lines
    tx_data = defaultdict(dict)
    tx_meta = {}

    for gtf, sample in zip(gtf_files, samples):
        with open(gtf) as fh:
            for line in fh:
                if line.startswith('#') or '\\ttranscript\\t' not in line:
                    continue
                fields = line.strip().split('\\t')
                attrs = fields[8]
                tx_id = re.search(r'transcript_id "([^"]+)"', attrs)
                cov   = re.search(r'cov "([^"]+)"', attrs)
                gene  = re.search(r'gene_id "([^"]+)"', attrs)
                ref   = re.search(r'reference_id "([^"]+)"', attrs)
                if tx_id and cov:
                    tid = tx_id.group(1)
                    tx_data[tid][sample] = cov.group(1)
                    if tid not in tx_meta:
                        tx_meta[tid] = {
                            'gene': gene.group(1) if gene else 'NA',
                            'ref':  ref.group(1) if ref else 'novel'
                        }

    with open('transcript_counts.tsv', 'w') as out:
        out.write('transcript_id\\tgene_id\\treference_id\\t' + '\\t'.join(samples) + '\\n')
        for tid in sorted(tx_data.keys()):
            meta = tx_meta[tid]
            vals = [tx_data[tid].get(s, '0') for s in samples]
            out.write(f"{tid}\\t{meta['gene']}\\t{meta['ref']}\\t" + '\\t'.join(vals) + '\\n')

    print(f"Written {len(tx_data)} transcripts x {len(samples)} samples")
    """
}
