#!/bin/bash

# Genome-based Phylogenetic Analysis Pipeline
# This script:
# 1. Runs eggNOG-mapper on s7_ctg000008c proteins
# 2. Maps OG assignments to Castelli et al. (2025) single-copy OG files
# 3. Extracts ~20 target taxa from each OG, adds the new MAG protein, aligns
# 4. Trims alignments with BMGE
# 5. Concatenates into a supermatrix
# 6. Removes compositionally biased sites (Munoz-Gomez et al. 2019)
# 7. Runs IQ-TREE 3 with ModelFinder (including mixture models)
#
# Goal: Place s7_ctg000008c within the Castelli et al. (2025) phylogenomic
# framework to verify its placement within/near Hepatincolaceae (WRAU01).

set -e

# Configuration
source /home/bdemedeiros/miniconda3/etc/profile.d/conda.sh

PROTEIN_FILE="checkm2_out/protein_files/s7_ctg000008c.faa"
CASTELLI_OG_DIR="castelli_et_al/single_ogs"
OUTPUT_DIR="genome_phylogenetic_analysis"
HELPER_SCRIPT="phylo_genome_helper.py"
ENV_NAME="symbiont_phylo"
ENV_FILE="symbiont_phylo_env.yml"
EGGNOG_DB_DIR="${EGGNOG_DB_DIR:-eggnog_data}"
THREADS=20
COMP_BIAS_THRESHOLD=0.30

echo "=================================================="
echo "Pipeline Step 06: Genome-based Phylogenetic Analysis"
echo "=================================================="
echo ""

# -------------------------------------------------------------------
# Step 0: Setup
# -------------------------------------------------------------------
echo "Checking required input files..."
for file in "$PROTEIN_FILE" "$HELPER_SCRIPT" "$ENV_FILE"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file"
        exit 1
    fi
done

if [[ ! -d "$CASTELLI_OG_DIR" ]]; then
    echo "ERROR: Castelli OG directory not found: $CASTELLI_OG_DIR"
    exit 1
fi
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
# Step 1: Download eggNOG database (if needed)
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 1: Checking eggNOG database"
echo "=================================================="

if [[ -f "$EGGNOG_DB_DIR/eggnog.db" && -f "$EGGNOG_DB_DIR/eggnog_proteins.dmnd" ]]; then
    echo "eggNOG database already present in $EGGNOG_DB_DIR"
else
    echo "Downloading eggNOG database to $EGGNOG_DB_DIR..."
    mkdir -p "$EGGNOG_DB_DIR"
    # download_eggnog_data.py has a known bug (issue #575) where the hostname
    # is wrong (eggnogdb.embl.de instead of eggnog5.embl.de), causing 404s.
    # Download manually with the corrected URL as a workaround.
    EGGNOG_BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"
    if [[ ! -f "$EGGNOG_DB_DIR/eggnog.db" ]]; then
        wget -c -O "$EGGNOG_DB_DIR/eggnog.db.gz" "$EGGNOG_BASE_URL/eggnog.db.gz"
        gunzip "$EGGNOG_DB_DIR/eggnog.db.gz"
    fi
    if [[ ! -f "$EGGNOG_DB_DIR/eggnog_proteins.dmnd" ]]; then
        wget -c -O "$EGGNOG_DB_DIR/eggnog_proteins.dmnd.gz" "$EGGNOG_BASE_URL/eggnog_proteins.dmnd.gz"
        gunzip "$EGGNOG_DB_DIR/eggnog_proteins.dmnd.gz"
    fi
    if [[ ! -f "$EGGNOG_DB_DIR/eggnog.taxa.db" ]]; then
        wget -c -O "$EGGNOG_DB_DIR/eggnog.taxa.tar.gz" "$EGGNOG_BASE_URL/eggnog.taxa.tar.gz"
        tar -xzf "$EGGNOG_DB_DIR/eggnog.taxa.tar.gz" -C "$EGGNOG_DB_DIR"
        rm -f "$EGGNOG_DB_DIR/eggnog.taxa.tar.gz"
    fi
    echo "eggNOG database download complete."
fi
echo ""

# -------------------------------------------------------------------
# Step 2: Run eggNOG-mapper
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 2: Running eggNOG-mapper"
echo "=================================================="
EGGNOG_OUTPUT="$OUTPUT_DIR/eggnog"

if [[ -s "${EGGNOG_OUTPUT}.emapper.annotations" ]]; then
    echo "eggNOG-mapper output already exists, skipping."
else
    echo "Running eggNOG-mapper on $PROTEIN_FILE..."
    # eggNOG-mapper's --iterate mode often hangs or causes bus errors with
    # large databases. Workaround: run diamond separately, then use emapper
    # with --no_search to annotate from the pre-computed hits.

    # Step 2a: Run diamond blastp directly (without --iterate)
    echo "  Running diamond search..."
    diamond blastp \
        -d "$EGGNOG_DB_DIR/eggnog_proteins.dmnd" \
        -q "$PROTEIN_FILE" \
        --threads "$THREADS" \
        -o "${EGGNOG_OUTPUT}.diamond.tmp" \
        -e 0.001 \
        --max-target-seqs 10 \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

    # Reformat to 4-column seed_orthologs format: query, hit, evalue, score
    echo "  Reformatting hits for eggNOG-mapper..."
    awk -F'\t' '{print $1"\t"$2"\t"$11"\t"$12}' "${EGGNOG_OUTPUT}.diamond.tmp" > "${EGGNOG_OUTPUT}.seed_orthologs"
    rm -f "${EGGNOG_OUTPUT}.diamond.tmp"

    # Step 2b: Run emapper annotation phase only
    echo "  Running eggNOG annotation..."
    emapper.py \
        -m no_search \
        --annotate_hits_table "${EGGNOG_OUTPUT}.seed_orthologs" \
        --no_file_comments \
        --output "$EGGNOG_OUTPUT" \
        --data_dir "$EGGNOG_DB_DIR" \
        --cpu "$THREADS"
    echo "eggNOG-mapper complete."
fi

NUM_ANNOTATED=$(grep -c -v "^#" "${EGGNOG_OUTPUT}.emapper.annotations" || true)
echo "Proteins annotated: $NUM_ANNOTATED"
echo ""

# -------------------------------------------------------------------
# Step 3: Map OG assignments to Castelli OG files
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 3: Mapping OG assignments to Castelli OG files"
echo "=================================================="
OG_MAPPING="$OUTPUT_DIR/og_mapping.tsv"

python "$HELPER_SCRIPT" map-ogs \
    --eggnog-annotations "${EGGNOG_OUTPUT}.emapper.annotations" \
    --castelli-og-dir "$CASTELLI_OG_DIR" \
    --output "$OG_MAPPING"

NUM_OGS=$(tail -n +2 "$OG_MAPPING" | wc -l)
echo "OGs mapped: $NUM_OGS / 179"
echo ""

# -------------------------------------------------------------------
# Step 4: Extract target taxa + new protein, re-align each OG
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 4: Extracting target taxa and aligning OGs"
echo "=================================================="
ALIGNED_DIR="$OUTPUT_DIR/aligned_ogs"

python "$HELPER_SCRIPT" extract-and-align \
    --og-mapping "$OG_MAPPING" \
    --castelli-og-dir "$CASTELLI_OG_DIR" \
    --protein-file "$PROTEIN_FILE" \
    --output-dir "$ALIGNED_DIR" \
    --threads "$THREADS"
echo ""

# -------------------------------------------------------------------
# Step 5: Trim with BMGE
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 5: Trimming alignments with BMGE"
echo "=================================================="
TRIMMED_DIR="$OUTPUT_DIR/trimmed_ogs"

python "$HELPER_SCRIPT" trim-alignments \
    --input-dir "$ALIGNED_DIR" \
    --output-dir "$TRIMMED_DIR"
echo ""

# -------------------------------------------------------------------
# Step 6: Concatenate
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 6: Concatenating trimmed alignments"
echo "=================================================="
CONCATENATED="$OUTPUT_DIR/concatenated.fasta"
PARTITIONS="$OUTPUT_DIR/partitions.txt"

python "$HELPER_SCRIPT" concatenate \
    --input-dir "$TRIMMED_DIR" \
    --output "$CONCATENATED" \
    --partition-file "$PARTITIONS"
echo ""

# -------------------------------------------------------------------
# Step 7: Remove compositionally biased sites
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 7: Removing compositionally biased sites (top ${COMP_BIAS_THRESHOLD})"
echo "=================================================="
DEBIASED="$OUTPUT_DIR/concatenated_debiased.fasta"

python "$HELPER_SCRIPT" remove-comp-bias \
    --alignment "$CONCATENATED" \
    --threshold "$COMP_BIAS_THRESHOLD" \
    --output "$DEBIASED"
echo ""

# -------------------------------------------------------------------
# Step 8: Run IQ-TREE 3
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 8: Running IQ-TREE 3 with ModelFinder"
echo "=================================================="
echo "  Testing LG, WAG, JTT + mixture models LG+C20+F+R, LG+C60+F+R"
echo "  1000 ultrafast bootstraps + SH-aLRT"
echo ""

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
# Step 9: Summary report
# -------------------------------------------------------------------
echo "=================================================="
echo "Step 9: Generating summary report"
echo "=================================================="

python "$HELPER_SCRIPT" summary \
    --output-dir "$OUTPUT_DIR"
echo ""

echo "=================================================="
echo "Pipeline Step 06 Complete"
echo "=================================================="
echo ""
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "Key output files:"
echo "  - eggnog.emapper.annotations   : eggNOG-mapper annotations"
echo "  - og_mapping.tsv               : OG-to-Castelli file mapping"
echo "  - aligned_ogs/                 : Per-OG alignments (MAFFT L-INS-i)"
echo "  - trimmed_ogs/                 : Per-OG trimmed alignments (BMGE)"
echo "  - concatenated.fasta           : Supermatrix (all OGs)"
echo "  - partitions.txt               : RAxML-style partition file"
echo "  - concatenated_debiased.fasta  : After compositional bias removal"
echo "  - genome_tree.treefile         : ML tree (Newick format)"
echo "  - genome_tree.iqtree           : IQ-TREE report"
echo "  - REPORT.md                    : Summary report"
echo ""
