#!/bin/bash

# Configuration
source /home/bdemedeiros/miniconda3/etc/profile.d/conda.sh

# Input and output directories
FASTA_FILE="final_meta_assembly_palm_weevil.asm.p_ctg.fa"
OUTPUT_DIR="meta_assembly"
ENV_NAME="symbiont_phylo"
ENV_FILE="symbiont_phylo_env.yml"

# Create or activate conda environment
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Creating conda environment '$ENV_NAME' from $ENV_FILE..."
    conda env create -f "$ENV_FILE"
fi
conda activate "$ENV_NAME"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Function to process sequences
process_fasta() {
    local fasta_file=$1
    local output_dir=$2
    local seq_name=""
    local seq_data=""
    local min_length=100000

    while read -r line; do
        if [[ $line == ">"* ]]; then
            if [[ -n $seq_name && ${#seq_data} -ge $min_length ]]; then
                echo -e ">$seq_name\n$seq_data" > "$output_dir/$seq_name.fasta"
            fi
            seq_name=$(echo "$line" | sed 's/^>//')
            seq_data=""
        else
            seq_data+="$line"
        fi
    done < "$fasta_file"

    # Handle the last sequence in the file
    if [[ -n $seq_name && ${#seq_data} -ge $min_length ]]; then
        echo -e ">$seq_name\n$seq_data" > "$output_dir/$seq_name.fasta"
    fi
}

# Call the function
process_fasta "$FASTA_FILE" "$OUTPUT_DIR"

echo "Processing complete. Files are in the $OUTPUT_DIR directory."

