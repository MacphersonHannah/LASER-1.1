#!/usr/bin/env bash

## Demultiplex PacBio Iso-Seq BAM files using Lima.
##
##   bash scripts/setup/demux_pacbio.sh input.bam primers.fa output_dir [--mas-seq]
##
## Arguments:
##   input.bam    raw or CCS BAM from PacBio sequencer
##   primers.fa   FASTA file containing 5' and 3' barcode sequences
##   output_dir   directory to write demultiplexed BAMs (created if needed;
##                  the primers.fa will be copied here so LASER can find it)
##   --mas-seq    optional flag: add if library is MAS-Iso-Seq
##
## Output:
##   output_dir/
##     *.bam              one BAM per barcode combination
##     *.bam.pbi          PacBio BAM indices
##     lima.counts        read counts per barcode pair
##     lima.summary       demux summary
##     primers.fa         copy of the primer file (required by LASER)

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: bash scripts/setup/demux_pacbio.sh input.bam primers.fa output_dir [--mas-seq]"
    exit 1
fi

INPUT_BAM="$1"
PRIMERS_FA="$2"
OUTPUT_DIR="$3"
MAS_SEQ=false

if [ "${4:-}" = "--mas-seq" ]; then
    MAS_SEQ=true
fi

# Validate inputs 
if [ ! -f "$INPUT_BAM" ]; then
    echo "ERROR: Input BAM not found: $INPUT_BAM"
    exit 1
fi

if [ ! -f "$PRIMERS_FA" ]; then
    echo "ERROR: Primers FASTA not found: $PRIMERS_FA"
    exit 1
fi

# Locate Lima 
# Use the LASER-bundled pacbio_prep_env conda env - should be at envs_prebuilt/pacbio_prep_env

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
LIMA_ENV="$PROJECT_ROOT/envs_prebuilt/pacbio_prep_env"

if command -v lima &> /dev/null; then
    LIMA_CMD="lima"
    echo "Using Lima from PATH: $(command -v lima)"
elif [ -f "$LIMA_ENV/bin/lima" ]; then
    LIMA_CMD="$LIMA_ENV/bin/lima"
    echo "Using Lima from LASER env: $LIMA_CMD"
else
    echo "ERROR: Lima not found on PATH or in $LIMA_ENV"
    echo "Run bash scripts/setup/setup_envs.sh first, or install Lima with:"
    echo "  conda install -c bioconda lima"
    exit 1
fi

# Set up output directory 
mkdir -p "$OUTPUT_DIR"

BASENAME=$(basename "$INPUT_BAM" .bam)
OUTPUT_BAM="$OUTPUT_DIR/${BASENAME}.bam"

echo ""
echo "=== PacBio demultiplexing with Lima ==="
echo "  Input BAM:   $INPUT_BAM"
echo "  Primers:     $PRIMERS_FA"
echo "  Output dir:  $OUTPUT_DIR"
echo "  MAS-Seq:     $MAS_SEQ"
echo ""

# Run Lima 
$LIMA_CMD \
    "$INPUT_BAM" \
    "$PRIMERS_FA" \
    "$OUTPUT_BAM" \
    --isoseq \
    --peek-guess \
    --num-threads 8 \
    --log-level INFO \
    --split-named

# Copy primers into output dir so LASER can find them 
# LASER requires the primer .fa to sit alongside the BAMs in the library dir
cp "$PRIMERS_FA" "$OUTPUT_DIR/$(basename "$PRIMERS_FA")"

echo ""
echo "=== Demultiplexing complete ==="
echo ""
echo "Output BAMs written to: $OUTPUT_DIR"
echo "Primer file copied to:  $OUTPUT_DIR/$(basename "$PRIMERS_FA")"
echo ""
echo "BAM files produced:"
ls "$OUTPUT_DIR"/*.bam 2>/dev/null | while read f; do echo "  $f"; done
echo ""
echo "To use in LASER, set in your run config file:"
echo "    params.library_dirs     = ['$(realpath "$OUTPUT_DIR")']"
echo "    params.start_from_demux = true"
