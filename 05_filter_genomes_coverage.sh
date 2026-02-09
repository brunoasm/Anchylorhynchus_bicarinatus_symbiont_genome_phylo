#!/bin/bash

# Filter genomes by quality and calculate read coverage
# This script:
# 1. Filters genomes with >=90% completeness and <5% contamination
# 2. Copies qualifying genomes to ncbi_submission_quality_filtered/
# 3. Maps reads to filtered genomes using minimap2
# 4. Calculates coverage statistics using samtools
# 5. Generates README files with GTDB-tk, CheckM2, and coverage stats

set -e

# Configuration
source /home/bdemedeiros/miniconda3/etc/profile.d/conda.sh

NCBI_DIR="ncbi_submission"
OUTPUT_DIR="ncbi_submission_quality_filtered"
CHECKM2_RESULTS="checkm2_out/quality_report.tsv"
GTDBTK_RESULTS="gtdbtk_out/classify/gtdbtk.bac120.summary.tsv"
READS="reads/no_chordata_contaminant_reads.fq.gz"
HELPER_SCRIPT="filter_genomes_helper.py"
ENV_NAME="symbiont_phylo"
ENV_FILE="symbiont_phylo_env.yml"

MIN_COMPLETENESS=90.0
MAX_CONTAMINATION=5.0
THREADS=20

echo "=================================================="
echo "Pipeline Step 05: Filter Genomes and Calculate Coverage"
echo "=================================================="
echo ""

# Check required files exist
echo "Checking required input files..."
for file in "$CHECKM2_RESULTS" "$GTDBTK_RESULTS" "$READS" "$HELPER_SCRIPT" "$ENV_FILE"; do
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

# Filter genomes
echo "Filtering genomes (completeness >= ${MIN_COMPLETENESS}%, contamination < ${MAX_CONTAMINATION}%)..."
PASSING_GENOMES=$(python "$HELPER_SCRIPT" filter \
    --checkm2 "$CHECKM2_RESULTS" \
    --min-completeness "$MIN_COMPLETENESS" \
    --max-contamination "$MAX_CONTAMINATION")

if [[ -z "$PASSING_GENOMES" ]]; then
    echo "WARNING: No genomes passed the quality filter!"
    echo "Consider adjusting the thresholds."
    exit 0
fi

echo "Genomes passing quality filter:"
echo "$PASSING_GENOMES" | while read -r genome; do
    echo "  - $genome"
done
echo ""

# Create output directory
echo "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Process each passing genome
echo "$PASSING_GENOMES" | while read -r GENOME; do
    echo ""
    echo "=================================================="
    echo "Processing genome: $GENOME"
    echo "=================================================="

    GENOME_DIR="$OUTPUT_DIR/$GENOME"
    FASTA_SRC="$NCBI_DIR/${GENOME}.fasta"
    FASTA_DST="$GENOME_DIR/${GENOME}.fasta"
    BAM_FILE="$GENOME_DIR/${GENOME}.sorted.bam"
    COVERAGE_FILE="$GENOME_DIR/${GENOME}.coverage.tsv"
    README_FILE="$GENOME_DIR/README.md"

    # Create genome subdirectory
    mkdir -p "$GENOME_DIR"

    # Copy FASTA
    echo "  Copying FASTA file..."
    cp "$FASTA_SRC" "$FASTA_DST"

    # Map reads with minimap2
    echo "  Mapping reads with minimap2 (PacBio HiFi mode)..."
    minimap2 -ax map-hifi -t "$THREADS" "$FASTA_DST" "$READS" 2>/dev/null | \
        samtools view -bS -@ "$THREADS" - | \
        samtools sort -@ "$THREADS" -o "$BAM_FILE" -

    # Index BAM
    echo "  Indexing BAM file..."
    samtools index "$BAM_FILE"

    # Calculate coverage
    echo "  Calculating coverage statistics..."
    samtools coverage "$BAM_FILE" > "$COVERAGE_FILE"

    # Generate README with all statistics
    echo "  Generating README..."
    python "$HELPER_SCRIPT" readme \
        --genome "$GENOME" \
        --checkm2 "$CHECKM2_RESULTS" \
        --gtdbtk "$GTDBTK_RESULTS" \
        --coverage "$COVERAGE_FILE" \
        --output "$README_FILE"

    echo "  Done processing $GENOME"
done

echo ""
echo "=================================================="
echo "Pipeline Step 05 Complete"
echo "=================================================="
echo ""
echo "Output directory: $OUTPUT_DIR/"
echo ""
echo "Summary of filtered genomes:"
ls -la "$OUTPUT_DIR"/*/README.md 2>/dev/null | while read -r line; do
    echo "  $line"
done
echo ""
echo "Each genome folder contains:"
echo "  - <genome>.fasta       : Genome sequence"
echo "  - <genome>.sorted.bam  : Read alignments"
echo "  - <genome>.sorted.bam.bai : BAM index"
echo "  - <genome>.coverage.tsv : Coverage statistics"
echo "  - README.md            : Combined statistics"
