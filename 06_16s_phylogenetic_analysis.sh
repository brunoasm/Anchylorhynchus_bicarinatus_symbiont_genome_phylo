#!/bin/bash

# 16S rRNA Phylogenetic Analysis Pipeline
# This script:
# 1. Extracts 16S rRNA sequence from annotated GenBank file
# 2. BLASTs against NCBI nt database (filtered for 16S sequences)
# 3. Retrieves sequences and metadata from NCBI
# 4. Creates annotated metadata table (symbiont/host information)
# 5. Aligns sequences with MAFFT
# 6. Builds phylogenetic tree with IQ-TREE using UNREST model for compositional heterogeneity

set -e

# Configuration
source /home/bdemedeiros/miniconda3/etc/profile.d/conda.sh

GENBANK_FILE="ncbi_submission_quality_filtered/s7_ctg000008c/s7_ctg000008c.bgpipe.output_838931.gb"
OUTPUT_DIR="16s_phylogenetic_analysis"
NUM_BLAST_HITS=100
NCBI_EMAIL="bdemedeiros@fas.harvard.edu"  # Required for NCBI Entrez API
NCBI_API_KEY="${NCBI_API_KEY:-}"  # Optional: allows 10 requests/sec instead of 3

HELPER_SCRIPT="phylo_16s_helper.py"
ENV_NAME="phylo_16s"
ENV_FILE="phylo_16s_env.yml"

echo "=================================================="
echo "Pipeline Step 06: 16S rRNA Phylogenetic Analysis"
echo "=================================================="
echo ""

# Check required files exist
echo "Checking required input files..."
for file in "$GENBANK_FILE" "$HELPER_SCRIPT" "$ENV_FILE"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    fi
done
echo "All required files found."
echo ""

# Create or activate conda environment
echo "Setting up conda environment..."
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Creating conda environment '$ENV_NAME' from $ENV_FILE..."
    conda env create -f "$ENV_FILE"
else
    echo "Conda environment '$ENV_NAME' already exists."
fi

conda activate "$ENV_NAME"
echo "Environment activated: $ENV_NAME"
echo ""

# Create output directory
echo "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"
echo ""

# Step 1: Extract 16S rRNA from GenBank
echo "=================================================="
echo "Step 1: Extracting 16S rRNA sequence"
echo "=================================================="
QUERY_FASTA="$OUTPUT_DIR/query_16s.fasta"

python "$HELPER_SCRIPT" extract-16s \
    --genbank "$GENBANK_FILE" \
    --output "$QUERY_FASTA" \
    --name "query_16s"
echo ""

# Verify extraction
if [[ ! -s "$QUERY_FASTA" ]]; then
    echo "ERROR: Failed to extract 16S rRNA sequence"
    exit 1
fi

# Step 2: BLAST against NCBI using command-line BLAST+ (more reliable than web API)
echo "=================================================="
echo "Step 2: Running BLAST against NCBI nt database"
echo "=================================================="
BLAST_XML="$OUTPUT_DIR/blast_results.xml"

echo "Running remote BLAST against nt database with 16S filter..."
echo "This may take 10-30 minutes depending on NCBI server load..."

# Use BLAST+ command line with remote option
# The nt database with entrez_query filter is more reliable than the 16S_ribosomal_RNA database
blastn -query "$QUERY_FASTA" \
    -db nt \
    -remote \
    -max_target_seqs "$NUM_BLAST_HITS" \
    -outfmt 5 \
    -out "$BLAST_XML" \
    -task megablast \
    -entrez_query "all[filter] AND biomol_genomic[PROP] AND 16S[Title]"

echo "BLAST completed: $BLAST_XML"
echo ""

# Verify BLAST completed
if [[ ! -s "$BLAST_XML" ]]; then
    echo "ERROR: BLAST search failed"
    exit 1
fi

# Step 3: Retrieve sequences and metadata
echo "=================================================="
echo "Step 3: Retrieving sequences and metadata from NCBI"
echo "=================================================="
BLAST_HITS_FASTA="$OUTPUT_DIR/blast_hits.fasta"
RAW_METADATA="$OUTPUT_DIR/blast_hits_raw_metadata.tsv"

python "$HELPER_SCRIPT" retrieve-sequences \
    --blast-xml "$BLAST_XML" \
    --fasta-output "$BLAST_HITS_FASTA" \
    --metadata-output "$RAW_METADATA" \
    --email "$NCBI_EMAIL" \
    ${NCBI_API_KEY:+--api-key "$NCBI_API_KEY"}
echo ""

# Verify retrieval
if [[ ! -s "$BLAST_HITS_FASTA" ]]; then
    echo "ERROR: Failed to retrieve sequences"
    exit 1
fi

# Step 4: Parse and annotate metadata
echo "=================================================="
echo "Step 4: Parsing and annotating metadata"
echo "=================================================="
ANNOTATED_METADATA="$OUTPUT_DIR/sequence_metadata.tsv"

python "$HELPER_SCRIPT" parse-metadata \
    --raw-metadata "$RAW_METADATA" \
    --output "$ANNOTATED_METADATA" \
    --email "$NCBI_EMAIL" \
    ${NCBI_API_KEY:+--api-key "$NCBI_API_KEY"}
echo ""

# Step 5: Combine query and hits
echo "=================================================="
echo "Step 5: Combining query and hit sequences"
echo "=================================================="
ALL_SEQUENCES="$OUTPUT_DIR/all_sequences.fasta"

cat "$QUERY_FASTA" "$BLAST_HITS_FASTA" > "$ALL_SEQUENCES"
echo "Combined sequences written to: $ALL_SEQUENCES"
NUM_SEQS=$(grep -c "^>" "$ALL_SEQUENCES")
echo "Total sequences: $NUM_SEQS"
echo ""

# Step 6: Align with MAFFT
echo "=================================================="
echo "Step 6: Aligning sequences with MAFFT"
echo "=================================================="
ALIGNMENT="$OUTPUT_DIR/alignment.fasta"
ALIGNMENT_RAW="$OUTPUT_DIR/alignment_raw.fasta"

echo "Running MAFFT alignment..."
mafft --auto --thread -1 "$ALL_SEQUENCES" > "$ALIGNMENT_RAW" 2>/dev/null

# Clean MAFFT output (remove header comments that IQ-TREE can't parse)
sed -n '/^>/,$p' "$ALIGNMENT_RAW" > "$ALIGNMENT"
rm "$ALIGNMENT_RAW"

echo "Alignment written to: $ALIGNMENT"
echo ""

# Step 7: Find outgroup
echo "=================================================="
echo "Step 7: Finding most distant sequence for rooting"
echo "=================================================="

# The outgroup is selected as the sequence with highest average pairwise distance
# to all other sequences in the alignment
OUTGROUP=$(python "$HELPER_SCRIPT" find-outgroup --alignment "$ALIGNMENT" | tail -n 1)
echo "Selected outgroup: $OUTGROUP"
echo ""

# Step 8: Build tree with IQ-TREE
echo "=================================================="
echo "Step 8: Building phylogenetic tree with IQ-TREE"
echo "=================================================="
TREE_PREFIX="$OUTPUT_DIR/16s_tree"

echo "Running IQ-TREE with extended model selection..."
echo "  Testing GTR and UNREST (non-reversible) models"
echo "  UNREST is recommended for endosymbionts with heterogeneous base composition"
echo "  Outgroup: $OUTGROUP"
echo "  Bootstrap replicates: 1000 (ultrafast)"

# Use extended model selection including UNREST for compositional heterogeneity
# -T AUTO: automatically determine optimal number of threads
# -mset GTR,UNREST: test both standard GTR and non-reversible UNREST models
# -mfreq F: use empirical base frequencies
# -mrate E,I,G,R: test equal rates, invariant sites, gamma, and FreeRate
iqtree -s "$ALIGNMENT" \
    -m MFP \
    -mset GTR,UNREST \
    -mfreq F \
    -mrate E,I,G,R \
    -B 1000 \
    -T AUTO \
    -o "$OUTGROUP" \
    --prefix "$TREE_PREFIX" \
    -redo

echo ""
echo "Tree built successfully!"
echo "  Tree file: ${TREE_PREFIX}.treefile"
echo "  IQ-TREE report: ${TREE_PREFIX}.iqtree"
echo ""

# Step 9: Generate report
echo "=================================================="
echo "Step 9: Generating analysis report"
echo "=================================================="

python "$HELPER_SCRIPT" summary \
    --output-dir "$OUTPUT_DIR" \
    --query-name "query_16s"
echo ""

echo "=================================================="
echo "Pipeline Step 06 Complete"
echo "=================================================="
echo ""
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "Key output files:"
echo "  - query_16s.fasta              : Extracted 16S from query genome"
echo "  - blast_results.xml            : Raw NCBI BLAST output"
echo "  - blast_hits.fasta             : Retrieved sequences (original IDs)"
echo "  - sequence_metadata.tsv        : Annotated metadata table"
echo "  - alignment.fasta              : MAFFT alignment"
echo "  - 16s_tree.treefile            : Rooted Newick tree"
echo "  - 16s_tree.iqtree              : IQ-TREE analysis report"
echo "  - REPORT.md                    : Summary report"
echo ""
