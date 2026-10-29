#!/usr/bin/env bash
# Usage: bash make_read_sample_map.sh <bam_dir> <output_tsv> <bam_glob_pattern>
# eg bash make_read_sample_map.sh pbmm2_flnc/ read_sample_map.tsv '*/*.mapped.bam'
# eg bash make_read_sample_map.sh minimap2_ont/ read_sample_map.tsv '*/*.sorted.bam'

set -euo pipefail

bam_dir="${1}"
out="${2:-read_sample_map.tsv}"
pattern="${3:-*/*.mapped.bam}"

echo -e "read_id\tsample" > "$out"
echo "Extracting read names from ${pattern} ..." >&2

# Validate at least one BAM exists before launching parallel jobs
first_bam=( "${bam_dir}"/${pattern} )
if [ ! -f "${first_bam[0]}" ]; then
    echo "No BAM files matched pattern: ${bam_dir}/${pattern}" >&2
    exit 1
fi

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# Launch all samtools view calls in parallel, each writing to a temp file
pids=()
for bam in "${bam_dir}"/${pattern}; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" .bam | sed 's/\.\(sorted\|flnc\.mapped\|mapped\)$//')
    echo "  sample: ${sample}" >&2
    tmpfile="${tmpdir}/${sample}.tsv"
    samtools view -@ 4 "$bam" | \
        awk -v s="$sample" 'BEGIN{OFS="\t"}{print $1,s}' > "$tmpfile" &
    pids+=($!)
done

# Wait for all parallel jobs
for pid in "${pids[@]}"; do wait "$pid"; done

# Concatenate temp files into output (order doesn't matter for lookup)
cat "${tmpdir}"/*.tsv >> "$out"
echo "Done. $(wc -l < "$out") read entries written." >&2