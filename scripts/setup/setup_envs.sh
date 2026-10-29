#!/usr/bin/env bash
set -uo pipefail  # Removed -e so it doesn't exit on errors

echo "=== Setting up LASER environments ==="

echo "Run this from the LASER-1.1 root directory (the one containing main.nf)." 


ENV_DIR="$PWD/envs_prebuilt"
mkdir -p "$ENV_DIR"

# Detect mamba or micromamba
if command -v mamba &> /dev/null; then
    MAMBA_CMD="mamba"
elif command -v micromamba &> /dev/null; then
    MAMBA_CMD="micromamba"
else
    echo "ERROR: mamba or micromamba not found. Please install one first."
    exit 1
fi

for tool in jq curl; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Installing $tool into base env (harmless general utility, not a fix for SQANTI3 specifically)..."
        $MAMBA_CMD install -n base -c conda-forge "$tool" -y
    fi
done

# Track results
FAILED_ENVS=()
SUCCESS_ENVS=()

create_env() {
    local name=$1
    local yaml=$2
    local path="$ENV_DIR/$name"

    if [ -d "$path" ]; then
        echo "Environment $name already exists - skipping."
        SUCCESS_ENVS+=("$name (existed)")
    else
        echo ""
        echo "=========================================="
        echo "Creating $name..."
        echo "=========================================="
        
        if $MAMBA_CMD env create -p "$path" -f "$yaml"; then
            echo "$name created successfully"
            SUCCESS_ENVS+=("$name")
        else
            echo "$name FAILED"
            FAILED_ENVS+=("$name")
        fi
    fi
}

# Core environments
create_env "pacbio_prep_env"      envs/pacbio_prep_env.yaml
create_env "core_utils_env"       envs/core_utils_env.yaml
create_env "isoform_collapse_env" envs/isoform_collapse_env.yaml
create_env "featurecounts_env"    envs/featurecounts_env.yaml
create_env "flair_env"            envs/flair_env.yaml
create_env "stats_r_env"          envs/stats_r_env.yaml


# Tool-specific environments (installed by setup_tools.sh)
SQANTI3_YAML=$(find tools/SQANTI3v5.2.2 -maxdepth 2 \( -iname "*.yml" -o -iname "*.yaml" \) 2>/dev/null | grep -i conda | head -1)
if [ -n "$SQANTI3_YAML" ]; then
    # bug fix: extra channels caused the solver to hang so kept just the two that matter
    if grep -qE "^\s*-\s*(anaconda|defaults|\br\b)" "$SQANTI3_YAML"; then
        echo "Patching ${SQANTI3_YAML}: trimming channels to conda-forge + bioconda (fixes solver hang)"
        sed -i '/^channels:/,/^[^[:space:]-]/{/^channels:/!{/^[^[:space:]-]/!d}}' "$SQANTI3_YAML"
        sed -i '/^channels:/a\  - conda-forge\n  - bioconda\n  - nodefaults' "$SQANTI3_YAML"
    fi
    create_env "SQANTI3" "$SQANTI3_YAML"
else
    echo "SQANTI3 FAILED: could not find a conda env yaml under tools/SQANTI3v5.2.2/ (checked *.yml/*.yaml matching 'conda'). Run setup_tools.sh first, or check the extracted folder manually."
    FAILED_ENVS+=("SQANTI3")
fi

create_env "ORFannotate" tools/ORFannotate/ORFannotate.conda_env.yml

echo ""
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo ""

if [ ${#SUCCESS_ENVS[@]} -gt 0 ]; then
    echo "Successfully created/found (${#SUCCESS_ENVS[@]}):"
    for env in "${SUCCESS_ENVS[@]}"; do
        echo "  - $env"
    done
fi

echo ""

if [ ${#FAILED_ENVS[@]} -gt 0 ]; then
    echo "FAILED (${#FAILED_ENVS[@]}):"
    for env in "${FAILED_ENVS[@]}"; do
        echo "  - $env"
    done
    echo ""
    echo "Fix the YAML files for failed environments and run again."
    exit 1
else
    echo "All environments created successfully!"
fi