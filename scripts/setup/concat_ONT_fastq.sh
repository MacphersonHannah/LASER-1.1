#!/usr/bin/env bash

## make executable chmod +x SCRIPT
## run: ./concat_fastq.sh fastq_pass merged_fastq my_output_name

set -euo pipefail

INPUT_DIR="$1"
OUTPUT_DIR="$2"
BASENAME="$3"

mkdir -p "$OUTPUT_DIR"

OUTFILE="$OUTPUT_DIR/${BASENAME}.fastq.gz"

echo "Concatenating FASTQ files from: $INPUT_DIR"
echo "Output file: $OUTFILE"

cat "$INPUT_DIR"/*.fastq.gz > "$OUTFILE"

echo "Done."