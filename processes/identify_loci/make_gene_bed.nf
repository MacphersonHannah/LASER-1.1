process MAKE_GENE_BED {
    tag { gtf.baseName }
    publishDir "${params.outdir}/preprocessing/make_gene_bed", mode: 'link'

    input:
        path gtf

    output:
        path "gene_loci_named.bed", emit: bed_named

    conda "${projectDir}/envs_prebuilt/core_utils_env"

    script:
    """
    set -euo pipefail

    # Make BED: chr, start, end, gene_name, gene_id, strand
    # Exclude mitochondrial (chrM, MT) and Y chromosome genes (chrY, Y)
    awk 'BEGIN{IGNORECASE=1}
         \$3 == "gene" &&
         \$1 !~ /(^chr)?m(t)?\$/ &&
         \$1 !~ /(^chr)?y\$/
    ' ${gtf} | \\
    awk 'BEGIN{OFS="\\t"}{
        gene_name="NA"; gene_id="NA";
        if (match(\$0, /gene_name "[^"]+"/)) {
            gene_name=substr(\$0, RSTART+11, RLENGTH-12)
        }
        if (match(\$0, /gene_id "[^"]+"/)) {
            gene_id=substr(\$0, RSTART+9, RLENGTH-10)
        }
        strand=\$7;
        if (strand != "+" && strand != "-") strand=".";
        print \$1, \$4-1, \$5, gene_name, gene_id, strand
    }' > gene_loci_named.bed
    """
}
