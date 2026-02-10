# Metagenome-assembled genomes from *Anchylorhynchus bicarinatus*

This repository contains the analysis pipeline for bacterial genomes assembled from metagenomic data associated with the weevil *Anchylorhynchus bicarinatus*. Total DNA was extracted from the whole body of the insect, and non-insect reads were isolated using BlobTools. Bacterial genomes were then assembled from this metagenomic fraction using PacBio HiFi reads.

The pipeline identifies, quality-filters, and characterizes metagenome-assembled genomes (MAGs), and places them phylogenetically using genome-wide protein markers and 16S rRNA sequences.

## Pipeline overview

The analysis proceeds through eight numbered shell scripts, each corresponding to a discrete step. Scripts are designed to be run sequentially and are self-contained: each creates its own conda environment (if needed) and checks for required inputs before proceeding.

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
    |-- [06] Genome-based phylogenetic analysis (18 taxa)
    v
genome_phylogenetic_analysis/
    |
    |-- [07] Enrich phylogeny with GTDB WRAU01 genomes (24 taxa)
    v
genome_phylogenetic_analysis/ + gtdb_genomes/
    |
    |-- [08] 16S rRNA phylogenetic analysis
    v
16s_phylogenetic_analysis/
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

Places the MAG within the phylogenomic framework of Castelli et al. (2025) to verify its relationship to Hepatincolaceae (WRAU01). Uses single-copy orthologous groups (OGs) from that study as a reference. Produces an initial 18-taxon tree; step 07 enriches this with additional GTDB genomes.

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

Compiles genome quality metrics (completeness, contamination, GC content, genome size) for all 24 taxa in the phylogenetic tree. Data are drawn from three sources:

1. **Castelli et al. (2025)** -- CheckM results from Supplementary Table 7 for 15 taxa
2. **This study (original)** -- CheckM (v1.2.4) run on 3 genomes:
   - `s7_ctg000008c` (our MAG)
   - `Terasakiella_pusilla_DSM_6293` (GCA_000688235.1)
   - `Thalassospira_profundimaris` (GCA_000300275.1)
3. **This study (GTDB enrichment)** -- CheckM (v1.2.4) run on 6 additional GTDB genomes (see Step 07)

All CheckM analyses use `lineage_wf` with `--reduced_tree` to ensure comparability with the Castelli et al. data (which used CheckM, not CheckM2).

The script also adds a `new_name` column for display in figures, with format `"GCA_XXXXXXX.1 organism_name"` for genomes with NCBI accessions (organism names are fetched from the NCBI Assembly database via Entrez), `"tip_name (Castelli et al, 2025)"` for genomes without NCBI accessions, and `"A. bicarinatus symbiont"` for our MAG. Additionally, the script adds `host_species`, `host_family`, and `host_source` columns for symbiont genomes. Host data are sourced from NCBI BioSample records where available, supplemented with data from Dittmer et al. (2023) for *Hepatincola* symbionts, Castelli et al. (2025) for *Tardigradibacter* and *Haliotis* symbiont species identification, and this study for the *A. bicarinatus* symbiont. Non-symbiont genomes receive "NA" for all three host columns.

When called with `--gtdb-manifest`, the script also incorporates the GTDB enrichment genomes from step 07 into the metadata table.

**Running this script requires internet access** to query NCBI.

**Note on CheckM vs CheckM2:** CheckM uses lineage-specific marker genes with phylogenetic placement, while CheckM2 uses machine learning. Results are generally comparable but may differ for divergent lineages. This study uses CheckM2 for initial quality filtering (Step 04) but CheckM for the final comparative table to match Castelli et al.'s methodology.

## Step 07: Enrich phylogeny with GTDB WRAU01 genomes

**Script:** `07_enrich_with_gtdb.sh`
**Helper:** `gtdb_enrichment_helper.py`

Enriches the 18-taxon phylogeny from step 06 with all additional genomes classified in GTDB order WRAU01. Queries the GTDB API, downloads protein and genomic FASTAs from NCBI, runs eggNOG-mapper, rebuilds the phylogeny with the expanded taxon set, and runs CheckM v1 for comparable quality metrics.

**Dependencies:** Same as step 06, plus CheckM v1 (in the `checkm` conda environment)

### Workflow

1. **Query GTDB** -- Searches the GTDB API for all genomes in order WRAU01. Classifies each as "existing" (already in tree), "duplicate", or "new". Generates a manifest JSON.

2. **Download proteins** -- Downloads protein FASTAs from NCBI FTP for new genomes. Falls back to Prodigal gene prediction when protein files are unavailable (common for European genome assemblies from ENA).

3. **eggNOG-mapper** -- Annotates new genomes with eggNOG (diamond blastp + emapper.py).

4. **OG mapping** -- Maps eggNOG annotations to Castelli single-copy OGs.

5. **Extract and align** -- Rebuilds all OG alignments with the enriched taxon set (uses `--extra-genomes` flag in `phylo_genome_helper.py`).

6--8. **Trim, concatenate, debias** -- Same as step 06.

9. **IQ-TREE** -- Rebuilds the phylogeny with the enriched alignment.

10. **Download genomic FASTAs** -- Downloads nucleotide assemblies for CheckM.

11. **CheckM v1** -- Runs `lineage_wf` on new genomes for comparable quality assessment.

12. **Metadata** -- Updates `phylogeny_checkm_stats.csv` with CheckM stats and host metadata.

### GTDB WRAU01 genome inventory

Of 8 genomes in GTDB order WRAU01, 6 were new additions:

| Accession | Tip Name | Host |
|-----------|----------|------|
| GCA_031256515.1 | WRAU01_MAG_031256515 | *Cornitermes pugnax* (Termitidae) |
| GCA_031274895.1 | WRAU01_MAG_031274895 | *Jugositermes tuberculatus* (Termitidae) |
| GCA_945888065.1 | WRAU01_MAG_945888065 | — |
| GCA_947460385.1 | WRAU01_MAG_947460385 | — |
| GCA_964020055.1 | Symbiont_of_P_viduata | *Pipizella viduata* (Syrphidae) |
| GCA_964020105.1 | Symbiont_of_E_torrentis | *Ecdyonurus torrentis* (Heptageniidae) |

### Output

Output is written to `genome_phylogenetic_analysis/` (overwriting step 06 results; git preserves history) and `gtdb_genomes/`.

| File | Description |
|------|-------------|
| `gtdb_genomes/genome_manifest.json` | GTDB genome inventory with metadata |
| `gtdb_genomes/checkm_output/quality_report.tsv` | CheckM v1 results for new genomes |
| `genome_phylogenetic_analysis/phylogeny_checkm_stats.csv` | Updated metadata (24 taxa) |
| `genome_phylogenetic_analysis/genome_tree.treefile` | ML tree (24 taxa) |
| `genome_phylogenetic_analysis/REPORT.md` | Summary report |

## Step 08: 16S rRNA phylogenetic analysis

**Script:** `08_16s_phylogenetic_analysis.sh`
**Helper:** `phylo_16s_helper.py`

Builds a 16S rRNA phylogeny for the *A. bicarinatus* symbiont by combining sequences from multiple sources: NCBI BLAST hits, SILVA nearest-neighbour search, and 16S sequences extracted from GTDB WRAU01-classified genomes.

**Dependencies:** Python >= 3.10, BioPython >= 1.81, BLAST+ >= 2.14, MAFFT >= 7.520, CD-HIT >= 4.8, IQ-TREE >= 2.2

### Workflow

1. **16S extraction** -- Extracts the 16S rRNA sequence from the annotated GenBank file.
2. **BLAST search** -- Queries the NCBI nt database for nearest neighbours.
3. **SILVA search** -- Incorporates results from a pre-downloaded SILVA SSU search.
4. **GTDB sequences** -- Retrieves 16S sequences from GTDB WRAU01-classified genome assemblies.
5. **Deduplication** -- Merges all sources and deduplicates with CD-HIT-EST at 99% identity.
6. **Metadata** -- Creates an annotated metadata table for all sequences.
7. **Alignment** -- Aligns with MAFFT (with reverse-complement detection).
8. **Tree inference** -- Builds a phylogeny with IQ-TREE using the UNREST model for compositional heterogeneity.

### Output

All output is written to `16s_phylogenetic_analysis/`.

## Conda environments

The pipeline uses three conda environments:

| Environment | Scripts | Environment file | Notes |
|-------------|---------|-----------------|-------|
| `gtdbtk-2.1.1` | 02 | (pre-existing) | Must be created manually before running step 02 |
| `checkm2` | 04 | (created by script) | Created automatically by `04_run_checkm2.sh`; also provides Prodigal for step 07 |
| `checkm` | 07, merge_checkm_stats.py | (created manually) | CheckM v1 for comparable quality stats |
| `symbiont_phylo` | 01, 03, 05, 06, 07, 08 | `symbiont_phylo_env.yml` | Created automatically if absent |

## Reproducibility

All analyses are driven by shell scripts and Python helpers that can be re-run from scratch. Each step that requires specific software defines its dependencies in a conda environment YAML file, and the corresponding shell script creates the environment automatically if it does not already exist.

Steps 06 and 07 accept an optional environment variable:
- `EGGNOG_DB_DIR` -- Path to eggNOG database directory (defaults to `eggnog_data/`). The database is downloaded automatically if not present.

Step 07 requires internet access to query the GTDB API and download genomes from NCBI.

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
| IQ-TREE | 3.0.1 | Phylogenetic inference and model selection |
| BioPython | >= 1.81 | Sequence parsing |
| Prodigal | >= 2.6 | Gene prediction (fallback for genomes without protein files) |
| BLAST+ | >= 2.14 | 16S rRNA sequence search |
| CD-HIT | >= 4.8 | Sequence deduplication |
