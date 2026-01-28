#!/bin/bash

# Run CheckM2 to assess genome quality for NCBI submission
# This script:
# 1. Creates a conda environment for CheckM2 (if it doesn't exist)
# 2. Runs CheckM2 on the genomes in ncbi_submission/

set -e

source /home/bdemedeiros/miniconda3/etc/profile.d/conda.sh

NCBI_DIR="ncbi_submission"
CHECKM2_OUT="checkm2_out"
ENV_NAME="checkm2"

# Check if conda environment exists, create if not
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Creating conda environment for CheckM2..."
    conda create -n "$ENV_NAME" -c bioconda -c conda-forge checkm2 -y
else
    echo "Conda environment '$ENV_NAME' already exists."
fi

# Activate environment
conda activate "$ENV_NAME"

# Check if CheckM2 database is set up
if ! checkm2 database --current 2>/dev/null; then
    echo "CheckM2 database not found. Downloading..."
    checkm2 database --download --path /home/bdemedeiros/checkm2_db
fi

# Run CheckM2
echo "Running CheckM2 on genomes in $NCBI_DIR..."
checkm2 predict \
    --input "$NCBI_DIR" \
    --output-directory "$CHECKM2_OUT" \
    --extension .fasta \
    --threads 20

echo ""
echo "CheckM2 analysis complete. Results in $CHECKM2_OUT/"
echo "Quality report: $CHECKM2_OUT/quality_report.tsv"
