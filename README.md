# Metagenome-assembled genomes from *Anchylorhynchus bicarinatus*

This repository contains the analysis pipeline for bacterial genomes assembled from metagenomic data associated with the weevil *Anchylorhynchus bicarinatus*. Total DNA was extracted from the whole body of the insect, and non-insect reads were isolated using BlobTools. Bacterial genomes were then assembled from this metagenomic fraction using PacBio HiFi reads.

The pipeline identifies, quality-filters, and characterizes metagenome-assembled genomes (MAGs), and places them phylogenetically using genome-wide protein markers.

## Pipeline overview

The analysis proceeds through six numbered shell scripts, each corresponding to a discrete step. Scripts are designed to be run sequentially and are self-contained: each creates its own conda environment (if needed) and checks for required inputs before proceeding.

```
final_meta_assembly_palm_weevil.asm.p_ctg.fa
    |
    |-- [01] Split into individual contigs (>=100 kb)
    v
meta_assembly/ (75 contigs)
    |
    |-- [02] Taxonomic classification with GTDB-Tk
    v
gtdbtk_out/
    |
    |-- [03] Select classified genomes
    v
ncbi_submission/ (6 genomes)
    |
    |-- [04] Quality assessment with CheckM2
    v
checkm2_out/
    |
    |-- [05] Filter by quality; map reads and compute coverage
    v
ncbi_submission_quality_filtered/ (1 genome: s7_ctg000008c)
    |
    |-- [06] Genome-based phylogenetic analysis
    v
genome_phylogenetic_analysis/
```

## Input data

| File | Description | NCBI Accession |
|------|-------------|----------------|
| `final_meta_assembly_palm_weevil.asm.p_ctg.fa` | PacBio HiFi metagenome assembly (primary contigs) | [JBTXKG000000000](https://www.ncbi.nlm.nih.gov/nuccore/JBTXKG000000000) |
| `reads/no_chordata_contaminant_reads.fq.gz` | PacBio HiFi reads with chordate contaminant reads removed | [SRR36582439](https://www.ncbi.nlm.nih.gov/sra/SRR36582439) |

These files are not tracked in the repository due to their size. Download them from NCBI before running the pipeline.

## Step 01: Split assembly into individual contigs

**Script:** `01_preprocess_metagenome.sh`

Splits the monolithic FASTA assembly into individual contig files, retaining only contigs >= 100 kb.

**Environment:** `symbiont_phylo`

**Output:** `meta_assembly/` (75 FASTA files)

## Step 02: Taxonomic classification

**Script:** `02_run_gtdbtk.sh`

Runs GTDB-Tk `classify_wf` on all contigs using the bacterial (bac120) and archaeal (ar53) marker gene sets. Uses `--min_perc_aa 20` and `--full_tree` to maximize sensitivity for divergent genomes.

**Environment:** `gtdbtk-2.1.1` (pre-existing conda environment)

**Output:** `gtdbtk_out/` (classification summaries and intermediate files)

## Step 03: Select classified genomes

**Script:** `03_prepare_ncbi_submission.sh`

Extracts genomes that received a valid GTDB-Tk classification (i.e., were assigned to at least a domain) and copies them to a submission directory, renaming files to replace periods with underscores for NCBI compatibility.

**Environment:** `symbiont_phylo`

**Output:** `ncbi_submission/` (6 classified genomes)

## Step 04: Quality assessment

**Script:** `04_run_checkm2.sh`

Runs CheckM2 on the classified genomes to assess completeness and contamination.

**Environment:** `checkm2` (created automatically if absent)

**Output:** `checkm2_out/quality_report.tsv`

## Step 05: Quality filtering and coverage

**Script:** `05_filter_genomes_coverage.sh`
**Helper:** `filter_genomes_helper.py`
**Environment file:** `symbiont_phylo_env.yml`

Filters genomes by quality (>= 90% completeness, < 5% contamination), maps PacBio HiFi reads back to each passing genome with minimap2, and computes coverage statistics with samtools. Generates a README for each genome combining GTDB-Tk taxonomy, CheckM2 quality metrics, and coverage statistics.

**Dependencies:** Python >= 3.10, minimap2 >= 2.26, samtools >= 1.17, pandas >= 2.0

**Output:** `ncbi_submission_quality_filtered/`

A single genome passed the quality filter:

| Genome | Taxonomy (GTDB-Tk) | Completeness | Contamination | Size | Coverage |
|--------|---------------------|-------------|---------------|------|----------|
| s7_ctg000008c | Bacteria; Pseudomonadota; Alphaproteobacteria; o\_\_WRAU01 | 98.93% | 0.00% | 1,231,242 bp | 67.4x |

This genome represents a novel lineage within Alphaproteobacteria (order WRAU01, as classified by GTDB-Tk using RED placement). It consists of a single circular contig at 100% breadth of coverage.

## Step 06: Genome-based phylogenetic analysis

**Script:** `06_genome_phylogenetic_analysis.sh`
**Helper:** `phylo_genome_helper.py`
**Environment file:** `symbiont_phylo_env.yml`

Places the MAG within the phylogenomic framework of Castelli et al. (2025) to verify its relationship to Hepatincolaceae (WRAU01). Uses single-copy orthologous groups (OGs) from that study as a reference.

**Dependencies:** Python >= 3.10, BioPython >= 1.81, pandas >= 2.0, eggNOG-mapper >= 2.1, diamond >= 2.0, MAFFT >= 7.520, BMGE >= 1.12, IQ-TREE >= 2.2

### Workflow

1. **eggNOG-mapper** -- Annotates predicted proteins from the MAG against the eggNOG database to assign orthologous groups.

2. **OG mapping** -- Maps eggNOG OG assignments to the 179 single-copy OG alignments from Castelli et al. (2025).

3. **Target taxa extraction** -- For each mapped OG, extracts ~20 representative taxa from the Castelli alignment plus the MAG protein sequence.

4. **Alignment** -- Aligns each OG with MAFFT L-INS-i.

5. **Trimming** -- Trims poorly aligned regions with BMGE (BLOSUM30 matrix).

6. **Concatenation** -- Concatenates all trimmed alignments into a supermatrix with a RAxML-style partition file.

7. **Compositional bias removal** -- Removes the most compositionally heterogeneous sites (top 30%) following Munoz-Gomez et al. (2019), as endosymbiont genomes often exhibit strong compositional biases that can mislead phylogenetic inference.

8. **Tree inference** -- Runs IQ-TREE 3 with ModelFinder testing LG, WAG, JTT and mixture models (LG+C20+F+R, LG+C60+F+R). Branch support assessed with 1000 ultrafast bootstraps and SH-aLRT.

9. **Summary report** -- Generates `REPORT.md` with alignment statistics, model selection results, and tree summary.

### Input data

Step 06 requires:

- `checkm2_out/protein_files/s7_ctg000008c.faa` -- Predicted proteins from CheckM2 (generated automatically in step 04)
- `castelli_et_al/single_ogs/` -- Single-copy OG alignments from Castelli et al. (2025). Not tracked in repository; obtain from the original study.
- `eggnog_data/` -- eggNOG database files (downloaded automatically if absent)

### Output

All output is written to `genome_phylogenetic_analysis/`:

| File | Description |
|------|-------------|
| `eggnog.emapper.annotations` | eggNOG-mapper annotations |
| `eggnog.seed_orthologs` | Diamond hits in seed_orthologs format |
| `og_mapping.tsv` | Mapping of MAG proteins to Castelli OG files |
| `aligned_ogs/` | Per-OG alignments (MAFFT L-INS-i) |
| `trimmed_ogs/` | Per-OG trimmed alignments (BMGE) |
| `concatenated.fasta` | Supermatrix (all OGs) |
| `partitions.txt` | RAxML-style partition file |
| `concatenated_debiased.fasta` | After compositional bias removal |
| `genome_tree.treefile` | Maximum-likelihood tree (Newick format) |
| `genome_tree.iqtree` | IQ-TREE report |
| `REPORT.md` | Summary report |
| `phylogeny_checkm_stats.csv` | CheckM statistics for all phylogeny tips (see below) |

### CheckM statistics for phylogeny tips

**Script:** `merge_checkm_stats.py`

Compiles genome quality metrics (completeness, contamination, GC content, genome size) for all 18 taxa in the phylogenetic tree. Data are drawn from two sources:

1. **Castelli et al. (2025)** -- CheckM results from Supplementary Table 7 for 15 taxa
2. **This study** -- CheckM (v1.2.4) run on 3 additional genomes:
   - `s7_ctg000008c` (our MAG)
   - `Terasakiella_pusilla_DSM_6293` (GCA_000688235.1)
   - `Thalassospira_profundimaris` (GCA_000300275.1)

The two outgroup genomes were downloaded from NCBI and analyzed with CheckM `lineage_wf` to ensure comparability with the Castelli et al. data (which used CheckM, not CheckM2).

The script also adds a `new_name` column for display in figures, with format `"GCA_XXXXXXX.1 organism_name"` for genomes with NCBI accessions (organism names are fetched from the NCBI Assembly database via Entrez), `"tip_name (Castelli et al, 2025)"` for genomes without NCBI accessions, and `"A. bicarinatus symbiont"` for our MAG. Additionally, the script adds `host_species`, `host_family`, and `host_source` columns for symbiont genomes. Host data are sourced from NCBI BioSample records where available, supplemented with data from Dittmer et al. (2023) for *Hepatincola* symbionts, Castelli et al. (2025) for *Tardigradibacter* and *Haliotis* symbiont species identification, and this study for the *A. bicarinatus* symbiont. Non-symbiont genomes receive "NA" for all three host columns. **Running this script requires internet access** to query NCBI.

**Note on CheckM vs CheckM2:** CheckM uses lineage-specific marker genes with phylogenetic placement, while CheckM2 uses machine learning. Results are generally comparable but may differ for divergent lineages. This study uses CheckM2 for initial quality filtering (Step 04) but CheckM for the final comparative table to match Castelli et al.'s methodology.

## Conda environments

The pipeline uses three conda environments:

| Environment | Scripts | Environment file | Notes |
|-------------|---------|-----------------|-------|
| `gtdbtk-2.1.1` | 02 | (pre-existing) | Must be created manually before running step 02 |
| `checkm2` | 04 | (created by script) | Created automatically by `04_run_checkm2.sh` |
| `checkm` | merge_checkm_stats.py | (created manually) | For comparable CheckM stats; see above |
| `symbiont_phylo` | 01, 03, 05, 06 | `symbiont_phylo_env.yml` | Created automatically if absent |

## Reproducibility

All analyses are driven by shell scripts and Python helpers that can be re-run from scratch. Each step that requires specific software defines its dependencies in a conda environment YAML file, and the corresponding shell script creates the environment automatically if it does not already exist.

Step 06 accepts an optional environment variable:
- `EGGNOG_DB_DIR` -- Path to eggNOG database directory (defaults to `eggnog_data/`). The database is downloaded automatically if not present.

Scripts were developed with the assistance of Claude Code (Anthropic).

## Software versions

| Software | Version | Purpose |
|----------|---------|---------|
| GTDB-Tk | 2.1.1 | Taxonomic classification |
| CheckM | 1.2.4 | Genome quality assessment (for cross-study comparison) |
| CheckM2 | (via conda) | Genome quality assessment (for initial filtering) |
| minimap2 | >= 2.26 | Read mapping (PacBio HiFi) |
| samtools | >= 1.17 | BAM processing and coverage |
| eggNOG-mapper | >= 2.1 | Functional annotation and OG assignment |
| diamond | >= 2.0 | Protein sequence search |
| MAFFT | >= 7.520 | Multiple sequence alignment |
| BMGE | >= 1.12 | Alignment trimming |
| IQ-TREE | >= 2.2 | Phylogenetic inference and model selection |
| BioPython | >= 1.81 | Sequence parsing |
