#!/bin/bash

# Enrich Phylogenetic Analysis with GTDB WRAU01 Genomes
#
# This script:
# 1. Queries GTDB for all genomes in order WRAU01
# 2. Downloads protein FASTAs for new genomes (not already in our tree)
# 3. Runs eggNOG-mapper on each new genome
# 4. Maps proteins to Castelli OGs
# 5. Rebuilds the full phylogeny with all genomes (original 18 + new)
# 6. Updates metadata CSV with CheckM stats and host info
#
# Prerequisites: Run 06_genome_phylogenetic_analysis.sh first (needs its outputs).

set -e

# Configuration
source /home/bdemedeiros/miniconda3/etc/profile.d/conda.sh

PROTEIN_FILE="checkm2_out/protein_files/s7_ctg000008c.faa"
CASTELLI_OG_DIR="castelli_et_al/single_ogs"
OUTPUT_DIR="genome_phylogenetic_analysis"
HELPER_SCRIPT="phylo_genome_helper.py"
GTDB_HELPER="gtdb_enrichment_helper.py"
MERGE_SCRIPT="merge_checkm_stats.py"
ENV_NAME="symbiont_phylo"
EGGNOG_DB_DIR="${EGGNOG_DB_DIR:-eggnog_data}"
GTDB_DIR="gtdb_genomes"
MANIFEST="$GTDB_DIR/genome_manifest.json"
THREADS=20
COMP_BIAS_THRESHOLD=0.30

echo "=================================================="
echo "Pipeline Step 07: Enrich Phylogeny with GTDB WRAU01 Genomes"
echo "=================================================="
echo ""

# -------------------------------------------------------------------
# Step 0: Setup
# -------------------------------------------------------------------
echo "Checking required input files..."
for file in "$PROTEIN_FILE" "$HELPER_SCRIPT" "$GTDB_HELPER" "$MERGE_SCRIPT"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    fi
done

if [[ ! -d "$CASTELLI_OG_DIR" ]]; then
    echo "ERROR: Castelli OG directory not found: $CASTELLI_OG_DIR"
    exit 1
fi

# Need the OG mapping from step 06
OG_MAPPING="$OUTPUT_DIR/og_mapping.tsv"
if [[ ! -f "$OG_MAPPING" ]]; then
    echo "ERROR: OG mapping not found: $OG_MAPPING"
    echo "Run 06_genome_phylogenetic_analysis.sh first."
    exit 1
fi

echo "All required files found."
echo ""

conda activate "$ENV_NAME"
echo "Environment activated: $ENV_NAME"
echo ""

# -------------------------------------------------------------------
# Step 1: Query GTDB for WRAU01 genomes
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 1: Querying GTDB for WRAU01 genomes"
echo "=================================================="

if [[ -f "$MANIFEST" ]]; then
    echo "Manifest already exists: $MANIFEST"
    echo "Delete it to re-query GTDB."
else
    python "$GTDB_HELPER" query-gtdb --output-dir "$GTDB_DIR"
fi
echo ""

# Count genomes by status
N_NEW=$(python -c "import json; m=json.load(open('$MANIFEST')); print(sum(1 for g in m['genomes'] if g['status']=='new'))")
echo "New genomes to add: $N_NEW"
echo ""

# -------------------------------------------------------------------
# Step 2: Download protein FASTAs for new genomes
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 2: Downloading protein FASTAs"
echo "=================================================="

python "$GTDB_HELPER" download-proteins --manifest "$MANIFEST"
echo ""

# -------------------------------------------------------------------
# Step 3: Run eggNOG-mapper for new genomes
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 3: Running eggNOG-mapper for new genomes"
echo "=================================================="

python "$GTDB_HELPER" run-eggnog \
    --manifest "$MANIFEST" \
    --eggnog-db "$EGGNOG_DB_DIR" \
    --threads "$THREADS"
echo ""

# -------------------------------------------------------------------
# Step 4: Build OG mappings for new genomes
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 4: Building OG mappings for new genomes"
echo "=================================================="

python "$GTDB_HELPER" build-og-mappings \
    --manifest "$MANIFEST" \
    --castelli-og-dir "$CASTELLI_OG_DIR"
echo ""

# -------------------------------------------------------------------
# Step 5: Extract target taxa + all genomes, align
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 5: Extracting target taxa and aligning OGs (enriched)"
echo "=================================================="
ALIGNED_DIR="$OUTPUT_DIR/aligned_ogs"

# Remove previous alignments to rebuild with enriched set
if [[ -d "$ALIGNED_DIR" ]]; then
    echo "Removing previous alignments..."
    rm -rf "$ALIGNED_DIR"
fi

python "$HELPER_SCRIPT" extract-and-align \
    --og-mapping "$OG_MAPPING" \
    --castelli-og-dir "$CASTELLI_OG_DIR" \
    --protein-file "$PROTEIN_FILE" \
    --output-dir "$ALIGNED_DIR" \
    --threads "$THREADS" \
    --extra-genomes "$MANIFEST"
echo ""

# -------------------------------------------------------------------
# Step 6: Trim with BMGE
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 6: Trimming alignments with BMGE"
echo "=================================================="
TRIMMED_DIR="$OUTPUT_DIR/trimmed_ogs"

if [[ -d "$TRIMMED_DIR" ]]; then
    rm -rf "$TRIMMED_DIR"
fi

python "$HELPER_SCRIPT" trim-alignments \
    --input-dir "$ALIGNED_DIR" \
    --output-dir "$TRIMMED_DIR"
echo ""

# -------------------------------------------------------------------
# Step 7: Concatenate
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 7: Concatenating trimmed alignments"
echo "=================================================="
CONCATENATED="$OUTPUT_DIR/concatenated.fasta"
PARTITIONS="$OUTPUT_DIR/partitions.txt"

python "$HELPER_SCRIPT" concatenate \
    --input-dir "$TRIMMED_DIR" \
    --output "$CONCATENATED" \
    --partition-file "$PARTITIONS"
echo ""

# -------------------------------------------------------------------
# Step 8: Remove compositionally biased sites
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 8: Removing compositionally biased sites (top ${COMP_BIAS_THRESHOLD})"
echo "=================================================="
DEBIASED="$OUTPUT_DIR/concatenated_debiased.fasta"

python "$HELPER_SCRIPT" remove-comp-bias \
    --alignment "$CONCATENATED" \
    --threshold "$COMP_BIAS_THRESHOLD" \
    --output "$DEBIASED"
echo ""

# -------------------------------------------------------------------
# Step 9: Run IQ-TREE 3
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 9: Running IQ-TREE 3 with ModelFinder"
echo "=================================================="
echo "  Testing LG, WAG, JTT + mixture models LG+C20+F+R, LG+C60+F+R"
echo "  1000 ultrafast bootstraps + SH-aLRT"
echo ""

# Remove previous tree files to allow re-run
rm -f "$OUTPUT_DIR"/genome_tree.*

iqtree3 -s "$DEBIASED" \
    -m MFP \
    -mset LG,WAG,JTT \
    -madd LG+C20+F+R,LG+C60+F+R \
    -B 1000 \
    --alrt 1000 \
    -T AUTO \
    --prefix "$OUTPUT_DIR/genome_tree"

echo ""
echo "IQ-TREE complete."
echo "  Tree: $OUTPUT_DIR/genome_tree.treefile"
echo "  Report: $OUTPUT_DIR/genome_tree.iqtree"
echo ""

# -------------------------------------------------------------------
# Step 10: Download genomic FASTAs for CheckM
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 10: Downloading genomic FASTAs for CheckM"
echo "=================================================="

python "$GTDB_HELPER" download-genomic --manifest "$MANIFEST"
echo ""

# -------------------------------------------------------------------
# Step 11: Run CheckM v1 on new genomes
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 11: Running CheckM v1 on new genomes"
echo "=================================================="

# CheckM v1 is in a separate conda environment
CHECKM_BIN="$HOME/miniconda3/envs/checkm/bin/checkm"
if [[ ! -f "$CHECKM_BIN" ]]; then
    echo "ERROR: CheckM v1 not found at $CHECKM_BIN"
    exit 1
fi

python "$GTDB_HELPER" run-checkm \
    --manifest "$MANIFEST" \
    --threads "$THREADS"
echo ""

# -------------------------------------------------------------------
# Step 12: Update metadata CSV
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 12: Updating metadata CSV with GTDB genomes"
echo "=================================================="

python "$MERGE_SCRIPT" --gtdb-manifest "$MANIFEST"
echo ""

# -------------------------------------------------------------------
# Step 13: Summary report
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 13: Generating summary report"
echo "=================================================="

python "$HELPER_SCRIPT" summary \
    --output-dir "$OUTPUT_DIR"
echo ""

echo "=================================================="
echo "Pipeline Step 07 Complete"
echo "=================================================="
echo ""
echo "Output directory: $OUTPUT_DIR/"
echo "GTDB data directory: $GTDB_DIR/"
echo ""
echo "Key output files:"
echo "  - gtdb_genomes/genome_manifest.json  : GTDB genome inventory"
echo "  - aligned_ogs/                       : Per-OG alignments (enriched)"
echo "  - trimmed_ogs/                       : Per-OG trimmed alignments"
echo "  - concatenated.fasta                 : Supermatrix"
echo "  - concatenated_debiased.fasta        : After compositional bias removal"
echo "  - genome_tree.treefile               : ML tree (Newick format)"
echo "  - phylogeny_checkm_stats.csv         : Updated metadata"
echo "  - REPORT.md                          : Summary report"
echo ""
