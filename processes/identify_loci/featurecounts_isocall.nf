// Assigns per-read transcript labels via featureCounts (largestOverlap,
// -g transcript_id) against the isocall merged GTF, then converts CORE
// output to FLAIR read-map format for build_read_isoform_lookup.py.

process FEATURECOUNTS_ISOCALL_READ_MAP {
    tag { "all_samples" }
    publishDir "${params.outdir}/preprocessing/featurecounts_isocall", mode: 'link'

    input:
        path(bams)   // collected list of *.flnc.mapped.bam (all samples)
        path(gtf)    // isocall_merged.isoforms.gtf

    output:
        path("isocall.read_map.txt"), emit: read_map

    conda "${projectDir}/envs_prebuilt/featurecounts_env"

    script:
    def bam_args = (bams instanceof List)
        ? bams.collect { it.toString() }.join(' ')
        : bams.toString()
    """
    set -euo pipefail

    featureCounts \\
        -L \\
        --largestOverlap \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -t exon -g transcript_id \\
        --primary \\
        -s 1 \\
        -R CORE \\
        ${bam_args} \\
        -o fc_transcripts.txt \\
        2> fc_transcripts.log

    # merge all per-BAM CORE files into a single FLAIR read map

    cat *.featureCounts \\
        | awk -F'\\t' '!/^#/ && \$2=="Assigned" && NF>=4 {print \$4"\\t"\$1}' \\
        | sort -k1,1 \\
        | awk -F'\\t' '
            BEGIN { prev=""; reads="" }
            {
                if (\$1 == prev) {
                    reads = reads "," \$2
                } else {
                    if (prev != "") print prev "\\t" reads
                    prev = \$1
                    reads = \$2
                }
            }
            END { if (prev != "") print prev "\\t" reads }
        ' > isocall.read_map.txt

    n=\$(wc -l < isocall.read_map.txt)
    echo "[featurecounts_isocall] Transcripts with >= 1 assigned read: \$n" >&2
    cat fc_transcripts.log >&2
    """
}
