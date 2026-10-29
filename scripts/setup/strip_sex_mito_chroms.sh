#!/usr/bin/env bash
# strips Y and mitochondrial (Y/M/MT) contigs from a genome FASTA and its matching GTF, then reindexes the FASTA with samtools faidx.
#
# This is NOT required for LASER to run correctly as loci on excluded
# chromosomes are already skipped downstream but running this first will speed up some steps because
# of less sequence for aligners andnannotation tools to load and smaller intermediate files
#
# Usage: bash scripts/setup/strip_sex_mito_chroms.sh <in.fa> <in.gtf> <out_prefix>
#
# Matches both UCSC-style (chrY, chrM) and Ensembl-style (Y, MT) naming so should work for any standard refs.

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: bash scripts/setup/strip_sex_mito_chroms.sh <in.fa> <in.gtf> <out_prefix>" >&2
  exit 1
fi

IN_FA="$1"
IN_GTF="$2"
OUT_PREFIX="$3"

OUT_FA="${OUT_PREFIX}.fa"
OUT_GTF="${OUT_PREFIX}.gtf"

if ! command -v samtools &> /dev/null; then
  echo "ERROR: samtools not found on PATH. Activate an env that has it (eg pacbio_prep_env) first." >&2
  exit 1
fi

echo "=== Stripping Y/M/MT contigs ==="
echo "  in fasta : ${IN_FA}"
echo "  in gtf   : ${IN_GTF}"
echo "  out fasta: ${OUT_FA}"
echo "  out gtf  : ${OUT_GTF}"

#  FASTA: keep every contig whose header does NOT match Y/M/MT (either naming style)
awk '
  /^>/ {
    header = substr($1, 2)
    keep = (header !~ /^(chr)?(Y|M|MT)$/)
  }
  keep { print }
' "${IN_FA}" > "${OUT_FA}"

# GTF: drop any line whose first column (chromosome) matches Y/M/MT
awk 'BEGIN{FS="\t"; OFS="\t"} $1 !~ /^(chr)?(Y|M|MT)$/' "${IN_GTF}" > "${OUT_GTF}"

echo "Indexing stripped FASTA with samtools faidx..."
samtools faidx "${OUT_FA}"

removed_contigs=$(grep "^>" "${IN_FA}" | sed 's/^>//' | awk '{print $1}' | grep -E "^(chr)?(Y|M|MT)$" || true)
echo ""
echo "Done. Removed contig(s): $(echo "${removed_contigs}" | tr '\n' ' ')"
echo "Point your config at:"
echo "  params.genome      = \"$(realpath "${OUT_FA}")\""
echo "  params.gencode_gtf = \"$(realpath "${OUT_GTF}")\""
