#!/bin/bash

# 16S rRNA Phylogenetic Analysis Pipeline
# This script:
# 1. Extracts 16S rRNA sequence from annotated GenBank file
# 2. BLASTs against NCBI nt database (filtered for 16S sequences)
# 3. Retrieves nearest neighbours from SILVA SSU search results
# 4. Retrieves 16S sequences from GTDB WRAU01-classified genomes
# 5. Merges and deduplicates all sources (cd-hit-est at 99% identity)
# 6. Creates annotated metadata table (symbiont/host/microbial taxonomy)
# 7. Aligns sequences with MAFFT (with reverse-complement detection)
# 8. Builds phylogenetic tree with IQ-TREE using UNREST model for compositional heterogeneity

set -e

# Configuration
source /home/bdemedeiros/miniconda3/etc/profile.d/conda.sh

GENBANK_FILE="ncbi_submission_quality_filtered/s7_ctg000008c/s7_ctg000008c.bgpipe.output_838931.gb"
OUTPUT_DIR="16s_phylogenetic_analysis"
NUM_BLAST_HITS=100
NCBI_EMAIL="bdemedeiros@fas.harvard.edu"  # Required for NCBI Entrez API
NCBI_API_KEY="${NCBI_API_KEY:-}"  # Optional: allows 10 requests/sec instead of 3

# SILVA SSU search result (pre-downloaded from SILVA website)
SILVA_FASTA="${SILVA_FASTA:-}"  # Path to SILVA nearest-neighbour FASTA file

# GTDB MSA file (from GTDB-Tk classify step)
GTDB_MSA="${GTDB_MSA:-gtdbtk_out/classify/gtdbtk.bac120.msa.fasta.gz}"

HELPER_SCRIPT="phylo_16s_helper.py"
ENV_NAME="symbiont_phylo"
ENV_FILE="symbiont_phylo_env.yml"

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

# -------------------------------------------------------------------
# Step 1: Clean derived results (keep NCBI-fetched data to avoid re-downloading)
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 1: Cleaning derived results"
echo "=================================================="
# Files that depend on the combination of all sources are always regenerated.
# NCBI-fetched files (FASTA + metadata from Steps 3-8b) are preserved to avoid
# unnecessary re-downloading; individual retrieval steps skip if output exists.
for old_file in \
    "$OUTPUT_DIR/alignment.fasta" \
    "$OUTPUT_DIR/alignment_raw.fasta" \
    "$OUTPUT_DIR/all_sequences.fasta" \
    "$OUTPUT_DIR/all_sources_raw_metadata.tsv" \
    "$OUTPUT_DIR/deduplicated.fasta" \
    "$OUTPUT_DIR/deduplicated.fasta.clstr" \
    "$OUTPUT_DIR/sequence_metadata.tsv" \
    "$OUTPUT_DIR/16s_tree_rev.treefile" \
    "$OUTPUT_DIR/16s_tree_rev.iqtree" \
    "$OUTPUT_DIR/16s_tree_rev.log" \
    "$OUTPUT_DIR/16s_tree_rev.bionj" \
    "$OUTPUT_DIR/16s_tree_rev.ckp.gz" \
    "$OUTPUT_DIR/16s_tree_rev.contree" \
    "$OUTPUT_DIR/16s_tree_rev.mldist" \
    "$OUTPUT_DIR/16s_tree_rev.model.gz" \
    "$OUTPUT_DIR/16s_tree_rev.splits.nex" \
    "$OUTPUT_DIR/16s_tree_nonrev.treefile" \
    "$OUTPUT_DIR/16s_tree_nonrev.iqtree" \
    "$OUTPUT_DIR/16s_tree_nonrev.log" \
    "$OUTPUT_DIR/16s_tree_nonrev.bionj" \
    "$OUTPUT_DIR/16s_tree_nonrev.ckp.gz" \
    "$OUTPUT_DIR/16s_tree_nonrev.contree" \
    "$OUTPUT_DIR/16s_tree_nonrev.mldist" \
    "$OUTPUT_DIR/16s_tree_nonrev.model.gz" \
    "$OUTPUT_DIR/16s_tree_nonrev.rootstrap.nex" \
    "$OUTPUT_DIR/16s_tree_nonrev.splits.nex" \
    "$OUTPUT_DIR/16s_tree.treefile" \
    "$OUTPUT_DIR/16s_tree.iqtree" \
    "$OUTPUT_DIR/REPORT.md"; do
    if [[ -f "$old_file" ]]; then
        rm "$old_file"
    fi
done
echo "Cleaned derived output files (NCBI-fetched sequences preserved if present)."
echo ""

# -------------------------------------------------------------------
# Step 2: Extract 16S rRNA from GenBank
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 2: Extracting 16S rRNA sequence"
echo "=================================================="
QUERY_FASTA="$OUTPUT_DIR/query_16s.fasta"

if [[ -s "$QUERY_FASTA" ]]; then
    echo "Query FASTA already exists, skipping extraction."
else
    python "$HELPER_SCRIPT" extract-16s \
        --genbank "$GENBANK_FILE" \
        --output "$QUERY_FASTA" \
        --name "query_16s"
fi
echo ""

# Verify extraction
if [[ ! -s "$QUERY_FASTA" ]]; then
    echo "ERROR: Failed to extract 16S rRNA sequence"
    exit 1
fi

# -------------------------------------------------------------------
# Step 3: BLAST against NCBI (skip if blast_results.xml exists)
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 3: Running BLAST against NCBI nt database"
echo "=================================================="
BLAST_XML="$OUTPUT_DIR/blast_results.xml"

if [[ -s "$BLAST_XML" ]]; then
    echo "BLAST results already exist ($BLAST_XML), skipping BLAST search."
else
    echo "Running remote BLAST against nt database with 16S filter..."
    echo "This may take 10-30 minutes depending on NCBI server load..."

    blastn -query "$QUERY_FASTA" \
        -db nt \
        -remote \
        -max_target_seqs "$NUM_BLAST_HITS" \
        -outfmt 5 \
        -out "$BLAST_XML" \
        -task megablast \
        -entrez_query "all[filter] AND biomol_genomic[PROP] AND 16S[Title]"

    echo "BLAST completed: $BLAST_XML"
fi
echo ""

# Verify BLAST completed
if [[ ! -s "$BLAST_XML" ]]; then
    echo "ERROR: BLAST search failed"
    exit 1
fi

# -------------------------------------------------------------------
# Step 4: Retrieve BLAST sequences and metadata
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 4: Retrieving BLAST sequences and metadata from NCBI"
echo "=================================================="
BLAST_HITS_FASTA="$OUTPUT_DIR/blast_hits.fasta"
BLAST_RAW_METADATA="$OUTPUT_DIR/blast_hits_raw_metadata.tsv"

if [[ -s "$BLAST_HITS_FASTA" && -s "$BLAST_RAW_METADATA" ]]; then
    echo "BLAST sequences already retrieved, skipping."
else
    python "$HELPER_SCRIPT" retrieve-sequences \
        --blast-xml "$BLAST_XML" \
        --fasta-output "$BLAST_HITS_FASTA" \
        --metadata-output "$BLAST_RAW_METADATA" \
        --email "$NCBI_EMAIL" \
        ${NCBI_API_KEY:+--api-key "$NCBI_API_KEY"}
fi
echo ""

# Verify retrieval
if [[ ! -s "$BLAST_HITS_FASTA" ]]; then
    echo "ERROR: Failed to retrieve BLAST sequences"
    exit 1
fi

# -------------------------------------------------------------------
# Step 5: Parse SILVA file
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 5: Parsing SILVA nearest-neighbour results"
echo "=================================================="
SILVA_HITS_TSV="$OUTPUT_DIR/silva_hits.tsv"

if [[ -n "$SILVA_FASTA" && -f "$SILVA_FASTA" ]]; then
    python "$HELPER_SCRIPT" parse-silva \
        --silva-fasta "$SILVA_FASTA" \
        --output "$SILVA_HITS_TSV"
else
    echo "No SILVA FASTA provided (set SILVA_FASTA env var). Skipping SILVA source."
    # Create empty TSV with header
    echo -e "silva_id\taccession\tversion\tstart\tend\tsimilarity" > "$SILVA_HITS_TSV"
fi
echo ""

# -------------------------------------------------------------------
# Step 6: Retrieve SILVA sequences
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 6: Retrieving SILVA sequences from GenBank"
echo "=================================================="
SILVA_HITS_FASTA="$OUTPUT_DIR/silva_hits.fasta"
SILVA_RAW_METADATA="$OUTPUT_DIR/silva_hits_raw_metadata.tsv"

if [[ -s "$SILVA_HITS_FASTA" && -s "$SILVA_RAW_METADATA" ]]; then
    echo "SILVA sequences already retrieved, skipping."
else
    # Only retrieve if we have non-empty SILVA hits
    SILVA_COUNT=$(tail -n +2 "$SILVA_HITS_TSV" | wc -l)
    if [[ "$SILVA_COUNT" -gt 0 ]]; then
        python "$HELPER_SCRIPT" retrieve-silva-sequences \
            --silva-tsv "$SILVA_HITS_TSV" \
            --fasta-output "$SILVA_HITS_FASTA" \
            --metadata-output "$SILVA_RAW_METADATA" \
            --email "$NCBI_EMAIL" \
            ${NCBI_API_KEY:+--api-key "$NCBI_API_KEY"}
    else
        echo "No SILVA accessions to retrieve."
        touch "$SILVA_HITS_FASTA"
        echo -e "accession\torganism\thit_def\thost\tisolation_source\tstrain\tcountry\tidentity\talignment_length\tevalue\tbit_score\tsource" > "$SILVA_RAW_METADATA"
    fi
fi
echo ""

# -------------------------------------------------------------------
# Step 7: Parse GTDB WRAU01 genomes
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 7: Extracting WRAU01 genome accessions from GTDB"
echo "=================================================="
GTDB_ACCESSIONS="$OUTPUT_DIR/gtdb_wrau01_accessions.txt"

if [[ -f "$GTDB_MSA" ]]; then
    python "$HELPER_SCRIPT" parse-gtdb-wrau01 \
        --msa "$GTDB_MSA" \
        --output "$GTDB_ACCESSIONS"
else
    echo "GTDB MSA file not found ($GTDB_MSA). Skipping GTDB source."
    touch "$GTDB_ACCESSIONS"
fi
echo ""

# -------------------------------------------------------------------
# Step 8: Retrieve GTDB 16S sequences
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 8: Retrieving 16S from GTDB WRAU01 genomes"
echo "=================================================="
GTDB_HITS_FASTA="$OUTPUT_DIR/gtdb_hits.fasta"
GTDB_RAW_METADATA="$OUTPUT_DIR/gtdb_hits_raw_metadata.tsv"

GTDB_COUNT=$(wc -l < "$GTDB_ACCESSIONS" 2>/dev/null || echo 0)
if [[ "$GTDB_COUNT" -gt 0 ]]; then
    python "$HELPER_SCRIPT" retrieve-gtdb-16s \
        --accessions "$GTDB_ACCESSIONS" \
        --fasta-output "$GTDB_HITS_FASTA" \
        --metadata-output "$GTDB_RAW_METADATA" \
        --email "$NCBI_EMAIL" \
        --query-fasta "$QUERY_FASTA" \
        ${NCBI_API_KEY:+--api-key "$NCBI_API_KEY"}
else
    echo "No GTDB accessions to retrieve."
    touch "$GTDB_HITS_FASTA"
    echo -e "accession\torganism\thit_def\thost\tisolation_source\tstrain\tcountry\tidentity\talignment_length\tevalue\tbit_score\tsource" > "$GTDB_RAW_METADATA"
fi
echo ""

# -------------------------------------------------------------------
# Step 8b: Retrieve 16S from Castelli et al. 2025 genomes
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 8b: Retrieving 16S from Castelli et al. 2025 genomes"
echo "=================================================="
# Hepatincolaceae genomes from Castelli et al. 2025
# (https://pmc.ncbi.nlm.nih.gov/articles/PMC11724238/)
CASTELLI_ACCESSIONS="$OUTPUT_DIR/castelli_accessions.txt"
CASTELLI_HITS_FASTA="$OUTPUT_DIR/castelli_hits.fasta"
CASTELLI_RAW_METADATA="$OUTPUT_DIR/castelli_hits_raw_metadata.tsv"

cat > "$CASTELLI_ACCESSIONS" <<'EOF'
GCF_000688235.1
GCA_001510075.1
GCA_045504435.1
EOF

python "$HELPER_SCRIPT" retrieve-gtdb-16s \
    --accessions "$CASTELLI_ACCESSIONS" \
    --fasta-output "$CASTELLI_HITS_FASTA" \
    --metadata-output "$CASTELLI_RAW_METADATA" \
    --email "$NCBI_EMAIL" \
    --query-fasta "$QUERY_FASTA" \
    --source castelli \
    ${NCBI_API_KEY:+--api-key "$NCBI_API_KEY"}
echo ""

# -------------------------------------------------------------------
# Step 9: Merge all sources
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 9: Merging sequences from all sources"
echo "=================================================="
ALL_SEQUENCES="$OUTPUT_DIR/all_sequences.fasta"
ALL_RAW_METADATA="$OUTPUT_DIR/all_sources_raw_metadata.tsv"

python "$HELPER_SCRIPT" merge-sources \
    --fasta-inputs "$BLAST_HITS_FASTA" "$SILVA_HITS_FASTA" "$GTDB_HITS_FASTA" "$CASTELLI_HITS_FASTA" \
    --metadata-inputs "$BLAST_RAW_METADATA" "$SILVA_RAW_METADATA" "$GTDB_RAW_METADATA" "$CASTELLI_RAW_METADATA" \
    --query-fasta "$QUERY_FASTA" \
    --fasta-output "$ALL_SEQUENCES" \
    --metadata-output "$ALL_RAW_METADATA"

echo ""
NUM_SEQS=$(grep -c "^>" "$ALL_SEQUENCES")
echo "Total merged sequences: $NUM_SEQS"
echo ""

# -------------------------------------------------------------------
# Step 10: Deduplicate with cd-hit-est
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 10: Deduplicating sequences with cd-hit-est (99% identity)"
echo "=================================================="
DEDUPLICATED="$OUTPUT_DIR/deduplicated.fasta"

cd-hit-est -i "$ALL_SEQUENCES" -o "$DEDUPLICATED" \
    -c 0.99 -n 10 -aS 0.9 -T 0 -M 0

NUM_DEDUP=$(grep -c "^>" "$DEDUPLICATED")
echo "Before deduplication: $NUM_SEQS sequences"
echo "After deduplication:  $NUM_DEDUP sequences"
echo "Removed: $((NUM_SEQS - NUM_DEDUP)) sequences"
echo ""

# -------------------------------------------------------------------
# Step 11: Parse and annotate metadata
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 11: Parsing and annotating metadata"
echo "=================================================="
ANNOTATED_METADATA="$OUTPUT_DIR/sequence_metadata.tsv"

python "$HELPER_SCRIPT" parse-metadata \
    --raw-metadata "$ALL_RAW_METADATA" \
    --output "$ANNOTATED_METADATA" \
    --email "$NCBI_EMAIL" \
    --survivors-fasta "$DEDUPLICATED" \
    ${NCBI_API_KEY:+--api-key "$NCBI_API_KEY"}
echo ""

# -------------------------------------------------------------------
# Step 12: Align with MAFFT (with reverse-complement detection)
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 12: Aligning sequences with MAFFT"
echo "=================================================="
ALIGNMENT="$OUTPUT_DIR/alignment.fasta"
ALIGNMENT_RAW="$OUTPUT_DIR/alignment_raw.fasta"

echo "Running MAFFT alignment with --adjustdirection..."
mafft --auto --adjustdirection --thread -1 "$DEDUPLICATED" > "$ALIGNMENT_RAW" 2>/dev/null

# Log reverse-complemented sequences and clean _R_ prefix
RC_COUNT=0
while IFS= read -r line; do
    if [[ "$line" == ">_R_"* ]]; then
        RC_COUNT=$((RC_COUNT + 1))
        seq_name="${line#>_R_}"
        echo "  Reverse-complemented: $seq_name"
    fi
done < "$ALIGNMENT_RAW"

if [[ "$RC_COUNT" -gt 0 ]]; then
    echo "Total sequences reverse-complemented by MAFFT: $RC_COUNT"
fi

# Clean: remove _R_ prefix from sequence names and strip any header comments
sed -n '/^>/,$p' "$ALIGNMENT_RAW" | sed 's/^>_R_/>/' > "$ALIGNMENT"

echo "Alignment written to: $ALIGNMENT"
echo ""

# -------------------------------------------------------------------
# Step 13: Build tree with IQ-TREE
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 13: Building phylogenetic tree with IQ-TREE"
echo "=================================================="
TREE_PREFIX="$OUTPUT_DIR/16s_tree"

echo "Running two IQ-TREE analyses (reversible and non-reversible models)..."
echo "  IQ-TREE 3 requires separate runs for reversible and non-reversible models."
echo "  The best model will be selected by comparing BIC scores across both runs."
echo ""
echo "  Run 1: Reversible models (GTR variants)"
echo "  Run 2: Non-reversible models (UNREST + Lie Markov)"
echo "  Non-reversible models are recommended for endosymbionts with"
echo "  heterogeneous base composition and strand-asymmetric substitution patterns."
echo "  Bootstrap replicates: 1000 (ultrafast)"
echo ""

TREE_PREFIX_REV="${TREE_PREFIX}_rev"
TREE_PREFIX_NONREV="${TREE_PREFIX}_nonrev"

# Run 1: Reversible models (GTR)
echo "--- Run 1: Reversible models (GTR) ---"
iqtree -s "$ALIGNMENT" \
    -m MFP \
    -mset GTR \
    -mfreq F \
    -mrate E,I,G,R \
    -B 1000 \
    -T AUTO \
    --prefix "$TREE_PREFIX_REV" \
    -redo
echo ""

# Run 2: Non-reversible models (UNREST + Lie Markov)
echo "--- Run 2: Non-reversible models (UNREST + Lie Markov) ---"
iqtree -s "$ALIGNMENT" \
    -m MFP+LM \
    -mset UNREST \
    --nonrev-model \
    -mfreq F \
    -mrate E,I,G,R \
    -B 1000 \
    -T AUTO \
    --prefix "$TREE_PREFIX_NONREV" \
    -redo
echo ""

# Compare BIC scores and select the best model
echo "--- Comparing models across both runs ---"
BIC_REV=$(grep "^Best-fit model" "${TREE_PREFIX_REV}.iqtree" | head -1)
BIC_NONREV=$(grep "^Best-fit model" "${TREE_PREFIX_NONREV}.iqtree" | head -1)
echo "  Reversible:     $BIC_REV"
echo "  Non-reversible: $BIC_NONREV"

# Extract BIC values from the model tables (top line = best model by BIC)
BEST_BIC_REV=$(grep -A 3 "^List of models sorted by BIC" "${TREE_PREFIX_REV}.iqtree" | tail -1 | awk '{print $9}')
BEST_BIC_NONREV=$(grep -A 3 "^List of models sorted by BIC" "${TREE_PREFIX_NONREV}.iqtree" | tail -1 | awk '{print $9}')
echo "  Best reversible BIC:     $BEST_BIC_REV"
echo "  Best non-reversible BIC: $BEST_BIC_NONREV"

# Select the winner (lower BIC is better)
if python3 -c "import sys; sys.exit(0 if float('$BEST_BIC_NONREV') < float('$BEST_BIC_REV') else 1)"; then
    WINNER="nonrev"
    echo "  Winner: Non-reversible model (lower BIC)"
else
    WINNER="rev"
    echo "  Winner: Reversible model (lower BIC)"
fi

# Copy winning tree files to final names
WINNER_PREFIX="${TREE_PREFIX}_${WINNER}"
for ext in treefile iqtree log contree mldist model.gz splits.nex; do
    if [[ -f "${WINNER_PREFIX}.${ext}" ]]; then
        cp "${WINNER_PREFIX}.${ext}" "${TREE_PREFIX}.${ext}"
    fi
done
# Copy rootstrap.nex only from nonrev (reversible models don't produce it)
if [[ -f "${TREE_PREFIX_NONREV}.rootstrap.nex" ]]; then
    cp "${TREE_PREFIX_NONREV}.rootstrap.nex" "${TREE_PREFIX}.rootstrap.nex"
fi

echo ""
echo "Tree built successfully!"
echo "  Selected tree: ${TREE_PREFIX}.treefile"
echo "  IQ-TREE report: ${TREE_PREFIX}.iqtree"
echo "  Reversible report: ${TREE_PREFIX_REV}.iqtree"
echo "  Non-reversible report: ${TREE_PREFIX_NONREV}.iqtree"
echo ""

# -------------------------------------------------------------------
# Step 14: Generate report
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 14: Generating analysis report"
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
echo "  - blast_hits.fasta             : Retrieved BLAST sequences"
echo "  - silva_hits.fasta             : Retrieved SILVA sequences"
echo "  - gtdb_hits.fasta              : Retrieved GTDB 16S sequences"
echo "  - castelli_hits.fasta          : Retrieved Castelli et al. 2025 16S sequences"
echo "  - all_sequences.fasta          : Merged sequences (all sources + query)"
echo "  - deduplicated.fasta           : After cd-hit-est 99% clustering"
echo "  - sequence_metadata.tsv        : Annotated metadata (survivors only)"
echo "  - alignment.fasta              : MAFFT alignment"
echo "  - 16s_tree.treefile            : Best ML tree (Newick format)"
echo "  - 16s_tree.iqtree              : IQ-TREE report for best model"
echo "  - 16s_tree_rev.iqtree          : IQ-TREE report (reversible models)"
echo "  - 16s_tree_nonrev.iqtree       : IQ-TREE report (non-reversible models)"
echo "  - REPORT.md                    : Summary report"
echo ""
