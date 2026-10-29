#!/usr/bin/env bash
set -euo pipefail

echo "=== Setting up LASER pipeline dependencies ==="

# Create main directories
mkdir -p tools
mkdir -p containers
mkdir -p containers/cache containers/tmp

# Clone or update tools
cd tools
echo "Checking ORFannotate..."
if [ -d ORFannotate ]; then
    echo "  ORFannotate already present - skipping clone."
else
    git clone --branch v1.1.3 https://github.com/egustavsson/ORFannotate.git
fi

echo "Checking SQANTI3..."
SQANTI_DIR="SQANTI3v5.2.2"
SQANTI_VERSION="v5.2.2"
if [ -d "$SQANTI_DIR" ]; then
    echo "  ${SQANTI_DIR} already present - skipping download."
else
    mkdir -p "$SQANTI_DIR"
    # Release asset naming isn't consistent across every SQANTI3 version so tries the direct release zip first (fast, matches recent versions convention), otherwise tries GitHub's tag archive
    if wget -q "https://github.com/ConesaLab/SQANTI3/releases/download/${SQANTI_VERSION}/SQANTI3_${SQANTI_VERSION}.zip" -O sqanti_dl.zip \
        && [ -s sqanti_dl.zip ]; then
        DL_FORMAT="zip"
    else
        echo "  Release zip asset not found for ${SQANTI_VERSION} - falling back to the tag source archive."
        rm -f sqanti_dl.zip
        wget -q "https://github.com/ConesaLab/SQANTI3/archive/refs/tags/${SQANTI_VERSION}.tar.gz" -O sqanti_dl.tar.gz
        DL_FORMAT="tar"
    fi

    python3 << PYEOF
import zipfile, tarfile, os, shutil

fmt = "${DL_FORMAT}"
target = "${SQANTI_DIR}"

if fmt == "zip":
    with zipfile.ZipFile("sqanti_dl.zip") as z:
        z.extractall(target)
        for info in z.infolist():
            perm = info.external_attr >> 16
            if perm:
                os.chmod(os.path.join(target, info.filename), perm)
else:
    # tarfile.extractall() preserves Unix permissions by default- no
    # manual chmod pass needed here, unlike the zip case above.
    with tarfile.open("sqanti_dl.tar.gz") as t:
        t.extractall(target)

entries = os.listdir(target)
if len(entries) == 1 and os.path.isdir(os.path.join(target, entries[0])):
    nested = os.path.join(target, entries[0])
    for item in os.listdir(nested):
        shutil.move(os.path.join(nested, item), os.path.join(target, item))
    os.rmdir(nested)
PYEOF
    rm -f sqanti_dl.zip sqanti_dl.tar.gz

    if [ ! -e "${SQANTI_DIR}/sqanti3_qc.py" ] && [ ! -e "${SQANTI_DIR}/src/sqanti3_qc.py" ]; then
        echo "  WARNING: sqanti3_qc.py not found where expected after extraction - inspect ${SQANTI_DIR}/ manually before proceeding." >&2
    fi

    # bug fix: SQANTI3's own code looks for gtfToGenePred at utilities/gtfToGenePred (relative to sqanti3_qc.py),
    # but some releases actually do it one level deeper at src/utilities/gtfToGenePred
    # Symlink it into the location the code expects, if needed.
    if [ ! -e "${SQANTI_DIR}/utilities/gtfToGenePred" ] && [ -e "${SQANTI_DIR}/src/utilities/gtfToGenePred" ]; then
        echo "  gtfToGenePred found under src/utilities/ but not utilities/ - creating a symlink so SQANTI3's own lookup finds it."
        chmod +x "${SQANTI_DIR}/src/utilities/gtfToGenePred"
        mkdir -p "${SQANTI_DIR}/utilities"
        ln -sf "../src/utilities/gtfToGenePred" "${SQANTI_DIR}/utilities/gtfToGenePred"
    fi
fi

# Patch the downloaded env yaml so jq/curl are installed INTO the SQANTI3 env itself, not just somewhere on the system
SQANTI_YAML=$(find "$SQANTI_DIR" -iname "*.yml" -o -iname "*.yaml" 2>/dev/null | grep -i conda | head -1)
if [ -n "$SQANTI_YAML" ]; then
    for dep in jq curl; do
        if ! grep -qE "^\s*-\s*${dep}\b" "$SQANTI_YAML"; then
            echo "  Patching ${SQANTI_YAML}: adding '${dep}' as a dependency (fixes a known Bioconductor post-link script failure)"
            sed -i "/^dependencies:/a\\  - ${dep}" "$SQANTI_YAML"
        fi
    done
else
    echo "  WARNING: could not find SQANTI3's conda env yaml to patch - if env creation later fails with a 'jq: command not found' post-link error, add 'jq' and 'curl' to that yaml's dependencies manually."
fi

cd ..


# Singularity or Apptainer check - NB apptainer is the new name

echo "Checking for Singularity or Apptainer..."
if ! command -v singularity &> /dev/null && ! command -v apptainer &> /dev/null; then
    echo "ERROR: Singularity or Apptainer not found in PATH."
    echo "Please load it on your system (eg, 'module load singularity') or install it."
    echo "See: https://apptainer.org/docs/admin/main/installation.html"
    exit 1
else
    if command -v singularity &> /dev/null; then
        echo "Singularity found: $(singularity --version)"
        SING_CMD="singularity"
    else
        echo "Apptainer found: $(apptainer --version)"
        SING_CMD="apptainer"
    fi
fi

# Define temp/cache dirs
cd containers
mkdir -p tmp cache
export SINGULARITY_TMPDIR="$(pwd)/tmp"
export SINGULARITY_CACHEDIR="$(pwd)/cache"
echo "Using custom Singularity directories:"
echo "  SINGULARITY_TMPDIR   = $SINGULARITY_TMPDIR"
echo "  SINGULARITY_CACHEDIR = $SINGULARITY_CACHEDIR"


# Pull DeepVariant container
echo "Checking DeepVariant container..."
if [ ! -f deepvariant_1.9.0.sif ]; then
    echo "Pulling DeepVariant 1.9.0 SIF from DockerHub..."
    $SING_CMD pull deepvariant_1.9.0.sif docker://google/deepvariant:1.9.0
else
    echo "  DeepVariant SIF already exists - skipping pull."
fi


# Build DeepVariant sandbox
if [ ! -d deepvariant_1.9.0_sandbox ]; then
    echo "Building DeepVariant sandbox (first time setup only)..."
    $SING_CMD build --sandbox deepvariant_1.9.0_sandbox deepvariant_1.9.0.sif
    echo "DeepVariant sandbox built successfully."
else
    echo "  DeepVariant sandbox already exists - skipping rebuild."
fi

# Pull Clair3-RNA container
# Clair3-RNA is used for variant calling in ONT direct RNA workflows.
# The container includes pre-trained models for ONT R9/R10 and PacBio HiFi.
echo "Checking Clair3-RNA container..."
if [ ! -f clair3_rna_latest.sif ]; then
    echo "Pulling Clair3-RNA SIF from DockerHub..."
    $SING_CMD pull clair3_rna_latest.sif docker://hkubal/clair3-rna:latest
    echo "Clair3-RNA container pulled successfully."
else
    echo "  Clair3-RNA SIF already exists - skipping pull."
fi

# Finish

cd ..
echo ""
echo "=== Setup complete ==="
echo ""
echo "In your run config file, set the following under the process blocks:"
echo ""
echo "  DeepVariant (PacBio / ONT cDNA):"
echo "    container = \"\${projectDir}/containers/deepvariant_1.9.0_sandbox\""
echo ""
echo "  Clair3-RNA (ONT direct RNA):"
echo "    container = \"\${projectDir}/containers/clair3_rna_latest.sif\""
echo ""
echo "Available Clair3-RNA models (check inside container):"
echo "    singularity exec containers/clair3_rna_latest.sif ls /opt/models/"
echo "  Common models:"
echo "    ont_dorado_drna004        (ONT direct RNA R10, recommended)"
echo "    ont_r10_dorado_cdna       (ONT cDNA R10)"
echo "    ont_r9_guppy_cdna         (ONT cDNA R9)"
echo "    hifi_mas_minimap2         (PacBio MAS-Seq, minimap2 aligned)"
echo "    hifi_mas_pbmm2            (PacBio MAS-Seq, pbmm2 aligned)"
