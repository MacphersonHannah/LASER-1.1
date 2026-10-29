// FLAIR quantify for direct RNA sequencing.
// Key differences from cDNA version:
//   - No --stringent: reads are fragmented, can't require full-length coverage
//   - No --check_splice: splice site coverage requirements too strict for direct RNA
//   - Keep --trust_ends: 3' end IS reliable in direct RNA (polyadenylated end)
//     This helps with TES definition for isoform distinction
//   - Keep --quality 0: same reasoning as cDNA

process FLAIR_QUANTIFY_DRNA {

    tag { sample }
    publishDir "${params.outdir}/preprocessing/flair/quantify/${sample}", mode: 'link'

    conda "${projectDir}/envs_prebuilt/flair_env"

    input:
        tuple val(sample), path(fastq)
        path(isoforms_fa)
        path(isoforms_bed)

    output:
        tuple val(sample), path("${sample}.flair.counts.tsv"),     emit: counts
        tuple val(sample), path("${sample}.isoform.read.map.txt"), emit: read_map

    script:
    """
    echo -e "${sample}\tsample\tbatch\t${fastq}" > manifest.tsv

    flair quantify \\
        --isoforms ${isoforms_fa} \\
        --reads_manifest manifest.tsv \\
        --isoform_bed ${isoforms_bed} \\
        --output ${sample}.flair \\
        --sample_id_only \\
        --quality 0 \\
        --trust_ends \\
        --generate_map \\
        --threads ${task.cpus}

    mv ${sample}.flair.*.isoform.read.map.txt \\
       ${sample}.isoform.read.map.txt
    """
}
