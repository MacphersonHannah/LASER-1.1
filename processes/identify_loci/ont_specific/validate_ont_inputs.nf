// Fast ONT fastq validation - checks file exists, is gzipped, and the first
// record is well-formed. Does NOT stream the whole file.
//
// Truncation at the end is extremely rare and all downstream tools (pychopper,
// minimap2) will fail loudly with a clear error if the file is corrupt.
// Streaming a 50GB file just to check the last 4 lines is not worth it.
//
// Checks existence, gzip magic bytes, and first fastq record.

process VALIDATE_ONT_INPUTS {

    tag { sample }
    conda "${projectDir}/envs_prebuilt/pacbio_prep_env"

    input:
        tuple val(sample), path(fastq)

    output:
        tuple val(sample), path(fastq)

    script:
    """
    set -euo pipefail

    echo "[VALIDATE] ${sample} > ${fastq}"

    # 1. exists and non-empty  
    if [ ! -s "${fastq}" ]; then
        echo "ERROR: ${fastq} is empty or does not exist" >&2
        exit 1
    fi
    fsize=\$(stat -c%s "${fastq}")
    echo "[VALIDATE] file size: \$(numfmt --to=iec \$fsize)"

    # 2. gzip magic bytes (reads 2 bytes only)  
    magic=\$(head -c 2 "${fastq}" | od -An -tx1 | tr -d ' \n')
    if [ "\$magic" != "1f8b" ]; then
        echo "ERROR: ${fastq} does not appear to be gzip compressed (magic=\$magic)" >&2
        exit 1
    fi

    # 3. validate first fastq record (decompresses ~1KB then stops)
    first4=\$(zcat "${fastq}" 2>/dev/null | head -n 4 || true)

    if [ -z "\$first4" ]; then
        echo "ERROR: cannot decompress ${fastq} - file may be corrupt" >&2
        exit 1
    fi

    line1=\$(echo "\$first4" | sed -n '1p')
    line3=\$(echo "\$first4" | sed -n '3p')

    if [[ "\$line1" != @* ]]; then
        echo "ERROR: first line does not start with '@' - not a valid fastq" >&2
        echo "  First line: \$line1" >&2
        exit 1
    fi
    if [[ "\$line3" != +* ]]; then
        echo "ERROR: third line does not start with '+' - malformed fastq record" >&2
        echo "  Third line: \$line3" >&2
        exit 1
    fi

    # Extract sample name from read header for sanity check
    read_name=\$(echo "\$line1" | awk '{print substr(\$1,2)}')
    echo "[VALIDATE] first read: \$read_name"
    echo "[VALIDATE] PASSED: ${sample} (${fastq})"
    """
}
