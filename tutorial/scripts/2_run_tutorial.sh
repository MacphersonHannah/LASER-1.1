#!/bin/bash
#SBATCH -J LASER_tutorial
#SBATCH --cpus-per-task=16
#SBATCH --mem=54000M
#SBATCH --time=06:00:00

# LASER tutorial. Run from LASER-1.1/tutorial/scripts:
#     sbatch run_tutorial.sh
# Outputs go to LASER-1.1/run_tutorial/tutorial/

set -eo pipefail

# Load conda via its hook script by absolute path (NOT by adding miniforge/bin to PATH - breaks ORFANNOTATE)
source /home/hk627/rds/hpc-work/bin/miniforge3/etc/profile.d/conda.sh
conda activate base

module load singularity
export NXF_VER=25.04.8

echo "python:   $(command -v python) ($(python --version 2>&1))"
echo "nextflow: $(command -v nextflow)"
echo "conda:    $(command -v conda)"

bash ../../run_LASER.sh tutorial "$PWD/tutorial.config"
