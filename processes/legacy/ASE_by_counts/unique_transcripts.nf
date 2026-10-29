process UNIQUE_TRANSCRIPTS {

    tag { barcode }
    publishDir "${params.outdir}/ase_by_counts/unique_transcripts/${sample}", mode: 'copy', overwrite: true
    input:
        tuple val(barcode), path(sqlite_db), path(pass_csv)

    output:
        path("${barcode}_haplotype_unique.tsv"), emit: uniq_tsv

    conda "${projectDir}/envs_prebuilt/nf_runner_laser"

    script:
    """
    set -euo pipefail

    echo "[UNIQUE] barcode=${barcode}"
    echo "[UNIQUE] sqlite_db=${sqlite_db}"
    echo "[UNIQUE] pass_csv=${pass_csv}"

    # extract loci that are PASS for *this* barcode
    loci=\$(awk -F',' -v bc=${barcode} '
        NR==1 {
            for(i=1;i<=NF;i++) if(\$i==bc) col=i; next
        }
        toupper(\$col)=="PASS" {printf "%s ", \$1}
        ' ${pass_csv})

    echo "[UNIQUE] loci tokens: \${loci}"

    # if none, write empty header-only file and exit cleanly
    if [ -z "\$loci" ]; then
        echo -e "barcode\tlocus\tisoform_id\thaplotype\tcount\tproportion" > ${barcode}_haplotype_unique.tsv
        echo "[UNIQUE] no PASS loci for ${barcode} -> created empty file"
        exit 0
    fi

    python3 ${projectDir}/scripts/haplotype_unique.py \\
        --db ${sqlite_db} \\
        --barcode ${barcode} \\
        --loci \$loci \\
        --count-thresh ${params.unique_haplotype_count_threshold} \\
        --prop-thresh  ${params.unique_haplotype_proportion_threshold} \\
        --out ${barcode}_haplotype_unique.tsv
    """
}