#!/bin/bash

# Prepare NCBI submission folder with classified genomes from GTDB-Tk
# This script:
# 1. Creates the ncbi_submission directory
# 2. Extracts genomes that were successfully classified by GTDB-Tk
# 3. Copies and renames FASTA files (replacing periods with underscores)

set -e

# Configuration
source /home/bdemedeiros/miniconda3/etc/profile.d/conda.sh
ENV_NAME="symbiont_phylo"
ENV_FILE="symbiont_phylo_env.yml"

# Create or activate conda environment
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Creating conda environment '$ENV_NAME' from $ENV_FILE..."
    conda env create -f "$ENV_FILE"
fi
conda activate "$ENV_NAME"

GTDBTK_OUT="gtdbtk_out"
META_ASSEMBLY="meta_assembly"
NCBI_DIR="ncbi_submission"

# Create output directory
mkdir -p "$NCBI_DIR"

# Function to check if a classification is valid (not Unclassified)
is_classified() {
    local classification="$1"
    # Return true if classification starts with d__ and is not just "Unclassified"
    if [[ "$classification" == d__* ]] && [[ "$classification" != "Unclassified" ]]; then
        return 0
    fi
    return 1
}

# Process bacterial summary
echo "Processing bacterial classifications..."
if [[ -f "$GTDBTK_OUT/classify/gtdbtk.bac120.summary.tsv" ]]; then
    while IFS=$'\t' read -r genome classification rest; do
        # Skip header
        [[ "$genome" == "user_genome" ]] && continue

        if is_classified "$classification"; then
            src_file="$META_ASSEMBLY/${genome}.fasta"
            # Replace periods with underscores for NCBI compatibility
            new_name=$(echo "$genome" | tr '.' '_')
            dst_file="$NCBI_DIR/${new_name}.fasta"

            if [[ -f "$src_file" ]]; then
                echo "  Copying $genome -> ${new_name}.fasta"
                cp "$src_file" "$dst_file"
            else
                echo "  WARNING: Source file not found: $src_file"
            fi
        fi
    done < "$GTDBTK_OUT/classify/gtdbtk.bac120.summary.tsv"
fi

# Process archaeal summary
echo "Processing archaeal classifications..."
if [[ -f "$GTDBTK_OUT/classify/gtdbtk.ar53.summary.tsv" ]]; then
    while IFS=$'\t' read -r genome classification rest; do
        # Skip header
        [[ "$genome" == "user_genome" ]] && continue

        if is_classified "$classification"; then
            src_file="$META_ASSEMBLY/${genome}.fasta"
            # Replace periods with underscores for NCBI compatibility
            new_name=$(echo "$genome" | tr '.' '_')
            dst_file="$NCBI_DIR/${new_name}.fasta"

            if [[ -f "$src_file" ]]; then
                echo "  Copying $genome -> ${new_name}.fasta"
                cp "$src_file" "$dst_file"
            else
                echo "  WARNING: Source file not found: $src_file"
            fi
        fi
    done < "$GTDBTK_OUT/classify/gtdbtk.ar53.summary.tsv"
fi

# Count results
count=$(ls -1 "$NCBI_DIR"/*.fasta 2>/dev/null | wc -l)
echo ""
echo "Done! Copied $count classified genomes to $NCBI_DIR/"
echo "Files have been renamed to replace periods with underscores."
