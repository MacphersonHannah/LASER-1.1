#!/usr/bin/env bash
set -euo pipefail

BED=$1
BAM=$2
OUT=$3

# Sort BED for consistent behaviour
sorted_bed=$(mktemp)
sort -k1,1 -k2,2n "$BED" > "$sorted_bed"

# Ensure BAM index exists
if [ ! -f "${BAM}.bai" ]; then
    samtools index "$BAM"
fi

# All genes
all_genes=$(mktemp)
cut -f4 "$sorted_bed" | sort -u > "$all_genes"

# Genes with >=1 read overlap
expressed=$(mktemp)
bedtools intersect -u -a "$sorted_bed" -b "$BAM" |
    cut -f4 | sort -u > "$expressed"

# Counts
total=$(wc -l < "$all_genes")
exp=$(wc -l < "$expressed")
unexp=$((total - exp))

# Output
{
    echo -e "metric\tvalue"
    echo -e "total_genes\t$total"
    echo -e "expressed_genes\t$exp"
    echo -e "unexpressed_genes\t$unexp"
} > "$OUT"
