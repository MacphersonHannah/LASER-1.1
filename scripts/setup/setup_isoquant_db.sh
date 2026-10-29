#!/bin/bash
set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/annotation.gtf"
    exit 1
fi

GTF_PATH="$1"
DB_PATH="${GTF_PATH%.gtf}.db"

echo "=== Setting up IsoQuant gene database ==="

# Find isoquant conda environment
ISOQUANT_ENV=$(find envs_prebuilt -name "isoquant_v2" -type d 2>/dev/null | head -1)
if [ -z "$ISOQUANT_ENV" ]; then
    echo "ERROR: IsoQuant environment not found. Run setup_envs.sh first."
    exit 1
fi

echo "Found IsoQuant environment: $ISOQUANT_ENV"

echo "Converting GTF to database..."
echo "  Input:  $GTF_PATH"
echo "  Output: $DB_PATH"
echo "This may take 30-60 minutes for a full genome annotation..."

# Create database with optimisations
"$ISOQUANT_ENV/bin/python" << PYEOF
import gffutils
print("Starting database creation...")
db = gffutils.create_db(
    "$GTF_PATH",
    "$DB_PATH",
    force=True,
    merge_strategy='merge',
    disable_infer_genes=True,
    disable_infer_transcripts=True
)
print("Database created successfully!")
PYEOF

echo "Database ready: $DB_PATH"
