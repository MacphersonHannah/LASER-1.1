// Fast ONT unmapped-BAM validation - checks file exists, is a valid BAM,
// and reads are unmapped (no reference sequences in header).
//
// Checks existence, BAM magic bytes, samtools quickcheck, unmapped header, and unmapped flag.

process VALIDATE_ONT_BAM_INPUTS {

    tag { sample }
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path(bam)

    script:
    """
    set -euo pipefail

    echo "[VALIDATE_BAM] ${sample} > ${bam}"

    # 1. exists and non-empty
    if [ ! -s "${bam}" ]; then
        echo "ERROR: ${bam} is empty or does not exist" >&2
        exit 1
    fi
    fsize=\$(stat -c%s "${bam}")
    echo "[VALIDATE_BAM] file size: \$(numfmt --to=iec \$fsize)"

    # 2. BAM starts with BAM\\1
    # gzip-decompress first 4 bytes; BAM is bgzf (gzip-compatible)
    magic=\$(head -c 100 "${bam}" | gzip -dc 2>/dev/null | head -c 4 || true)
    if [ "\${magic:0:3}" != "BAM" ]; then
        echo "ERROR: ${bam} does not appear to be a BAM file" >&2
        exit 1
    fi

    # 3. samtools quickcheck 
        qc_msg=\$(samtools quickcheck "${bam}" 2>&1) || {
        if echo "\$qc_msg" | grep -qi "no targets"; then
            echo "[VALIDATE_BAM] quickcheck: no targets in header (expected for unmapped BAM)"
        else
            echo "ERROR: samtools quickcheck failed for ${bam} - file may be truncated or corrupt" >&2
            echo "ERROR: \$qc_msg" >&2
            exit 1
        fi
    }

    # 4. check header has no @SQ lines (unmapped BAM) 
    sq_count=\$(samtools view -H "${bam}" | grep -c '^@SQ' || true)
    if [ "\$sq_count" -gt 0 ]; then
        echo "WARNING: ${bam} has \$sq_count @SQ header lines - this looks like a mapped BAM." >&2
        echo "WARNING: ont_cdna_minimap_entry expects unmapped BAMs. Reads will be re-mapped with minimap2." >&2
        echo "WARNING: Continuing, but check your input if results look unexpected." >&2
    fi

    # 5. first read is unmapped (flag & 0x4)
    first_flag=\$(set +o pipefail; samtools view "${bam}" | head -n 1 | cut -f 2)
    if [ -n "\$first_flag" ]; then
        if [ \$((first_flag & 4)) -eq 0 ]; then
            echo "WARNING: first read in ${bam} appears to be mapped (flag=\$first_flag)." >&2
            echo "WARNING: ont_cdna_minimap_entry expects unmapped reads. Continuing anyway." >&2
        fi
    fi

    # 6. read count (quick estimate from first few reads)
    nreads=\$(samtools view -c "${bam}")
    echo "[VALIDATE_BAM] read count: \$nreads"

    echo "[VALIDATE_BAM] PASSED: ${sample} (${bam})"
    """
}
