#!/usr/bin/env bash
set -euo pipefail

echo "========================================"
echo "   LASER Pipeline - Complete Setup"
echo "========================================"
echo ""

# Check if GTF path provided
if [ $# -eq 0 ]; then
    echo "Usage: bash setup_all.sh /path/to/annotation.gtf"
    echo ""
    echo "This script will:"
    echo "  1. Download tools and containers (DeepVariant, Clair3-RNA, SQANTI3, ORFannotate)"
    echo "  2. Create all conda environments (~20 minutes)"
    echo "  3. Build IsoQuant gene database from your GTF (~30-60 minutes)"
    echo ""
    echo "Example:"
    echo "  bash setup_all.sh /home/user/references/gencode.v47_noY_noMT.annotation.gtf"
    echo ""
    echo "Requirements:"
    echo "  - conda/mamba installed and in PATH"
    echo "  - Singularity or Apptainer installed"
    echo "  - Internet connection"
    echo "  - ~50GB free disk space"
    exit 1
fi

GTF_PATH="$1"

# Check GTF exists
if [ ! -f "$GTF_PATH" ]; then
    echo "ERROR: GTF file not found: $GTF_PATH"
    exit 1
fi

# Check conda/mamba available
if ! command -v mamba &> /dev/null && ! command -v micromamba &> /dev/null && ! command -v conda &> /dev/null; then
    echo "ERROR: conda, mamba, or micromamba not found."
    echo "Please install conda first: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

START_TIME=$(date +%s)


# STEP 1: Setup Tools (containers, SQANTI3, ORFannotate)

echo ""
echo "========================================"
echo "STEP 1/3: Setting up tools and containers"
echo "========================================"
echo ""

if [ ! -f scripts/setup/setup_tools.sh ]; then
    echo "ERROR: setup_tools.sh not found. Run from LASER-1.1 root directory."
    exit 1
fi

bash scripts/setup/setup_tools.sh
if [ $? -ne 0 ]; then
    echo "ERROR: Tool setup failed. Check output above."
    exit 1
fi


# STEP 2: Setup Conda Environments

echo ""
echo "========================================"
echo "STEP 2/3: Creating conda environments"
echo "========================================"
echo "This will take ~20 minutes..."
echo ""

if [ ! -f scripts/setup/setup_envs.sh ]; then
    echo "ERROR: setup_envs.sh not found."
    exit 1
fi

bash scripts/setup/setup_envs.sh
if [ $? -ne 0 ]; then
    echo "ERROR: Environment setup failed. Check output above."
    exit 1
fi



# FINISH

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))

DB_PATH="${GTF_PATH%.gtf}.db"

echo ""
echo "========================================"
echo "  Setup Complete!"
echo "========================================"
echo ""
echo "Total time: ${HOURS}h ${MINUTES}m"
echo ""
echo "Created:"
echo " Tools and containers in containers/"
echo " Conda environments in envs_prebuilt/"
echo " IsoQuant database at:"
echo " $DB_PATH"
echo ""
echo "Next steps:"
echo "  1. Copy and edit a config from example_configs/:"
echo "     - Set params.gencode_gtf to your GENCODE GTF"
echo "     - Configure other reference paths (genome, cage, etc.)"
echo "     - Set params.library_dirs or params.ont_input_files"
echo ""
echo "  2. Activate the runner environment:"
echo "     conda activate envs_prebuilt/core_utils_env"
echo ""
echo "  3. Run LASER:"
echo "     bash run_LASER.sh MY_RUN_NAME"
echo ""