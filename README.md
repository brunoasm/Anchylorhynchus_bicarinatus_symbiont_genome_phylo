# Metagenome-assembled genomes from *Anchylorhynchus bicarinatus*

This repository contains the analysis pipeline for bacterial genomes assembled from metagenomic data associated with the weevil *Anchylorhynchus bicarinatus*. Total DNA was extracted from the whole body of the insect, and non-insect reads were isolated using BlobTools. Bacterial genomes were then assembled from this metagenomic fraction using PacBio HiFi reads.

The pipeline identifies, quality-filters, and characterizes metagenome-assembled genomes (MAGs), and places them phylogenetically using 16S rRNA.

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
    |-- [06] 16S rRNA extraction, BLAST, alignment, phylogenetic inference
    v
16s_phylogenetic_analysis/
```

## Input data

| File | Description | NCBI Accession |
|------|-------------|----------------|
| `final_meta_assembly_palm_weevil.asm.p_ctg.fa` | PacBio HiFi metagenome assembly (primary contigs) | [JBTXKG000000000](https://www.ncbi.nlm.nih.gov/nuccore/JBTXKG000000000) |
| `reads/no_chordata_contaminant_reads.fq.gz` | PacBio HiFi reads with chordate contaminant reads removed | [SRR36582439](https://www.ncbi.nlm.nih.gov/sra/SRR36582439) |

These files are not tracked in the repository due to their size. Download them from NCBI before running the pipeline.

The GenBank annotation file (`ncbi_submission_quality_filtered/s7_ctg000008c/s7_ctg000008c.bgpipe.output_838931.gb`) is tracked in the repository for convenience. It can also be retrieved from genome accession [JBTXKG000000000](https://www.ncbi.nlm.nih.gov/nuccore/JBTXKG000000000). This file is required as input to step 06.

## Step 01: Split assembly into individual contigs

**Script:** `01_preprocess_metagenome.sh`

Splits the monolithic FASTA assembly into individual contig files, retaining only contigs >= 100 kb.

**Output:** `meta_assembly/` (75 FASTA files)

## Step 02: Taxonomic classification

**Script:** `02_run_gtdbtk.sh`

Runs GTDB-Tk `classify_wf` on all contigs using the bacterial (bac120) and archaeal (ar53) marker gene sets. Uses `--min_perc_aa 20` and `--full_tree` to maximize sensitivity for divergent genomes.

**Environment:** `gtdbtk-2.1.1` (pre-existing conda environment)

**Output:** `gtdbtk_out/` (classification summaries and intermediate files)

## Step 03: Select classified genomes

**Script:** `03_prepare_ncbi_submission.sh`

Extracts genomes that received a valid GTDB-Tk classification (i.e., were assigned to at least a domain) and copies them to a submission directory, renaming files to replace periods with underscores for NCBI compatibility.

**Output:** `ncbi_submission/` (6 classified genomes)

## Step 04: Quality assessment

**Script:** `04_run_checkm2.sh`

Runs CheckM2 on the classified genomes to assess completeness and contamination.

**Environment:** `checkm2` (created automatically if absent)

**Output:** `checkm2_out/quality_report.tsv`

## Step 05: Quality filtering and coverage

**Script:** `05_filter_genomes_coverage.sh`
**Helper:** `filter_genomes_helper.py`
**Environment file:** `genome_filter_coverage_env.yml`

Filters genomes by quality (>= 90% completeness, < 5% contamination), maps PacBio HiFi reads back to each passing genome with minimap2, and computes coverage statistics with samtools. Generates a README for each genome combining GTDB-Tk taxonomy, CheckM2 quality metrics, and coverage statistics.

**Dependencies:** Python >= 3.9, minimap2 >= 2.26, samtools >= 1.17, pandas >= 2.0

**Output:** `ncbi_submission_quality_filtered/`

A single genome passed the quality filter:

| Genome | Taxonomy (GTDB-Tk) | Completeness | Contamination | Size | Coverage |
|--------|---------------------|-------------|---------------|------|----------|
| s7_ctg000008c | Bacteria; Pseudomonadota; Alphaproteobacteria; o\_\_WRAU01 | 98.93% | 0.00% | 1,231,242 bp | 67.4x |

This genome represents a novel lineage within Alphaproteobacteria (order WRAU01, as classified by GTDB-Tk using RED placement). It consists of a single circular contig at 100% breadth of coverage. The genome was subsequently annotated externally (GenBank file: `s7_ctg000008c.bgpipe.output_838931.gb`).

## Step 06: 16S rRNA phylogenetic analysis

**Script:** `06_16s_phylogenetic_analysis.sh`
**Helper:** `phylo_16s_helper.py`
**Environment file:** `phylo_16s_env.yml`

Extracts the 16S rRNA gene from the annotated GenBank file, searches for related sequences in NCBI, and builds a phylogenetic tree to place the organism among its closest relatives.

**Dependencies:** Python >= 3.10, BioPython >= 1.81, pandas >= 2.0, BLAST+ >= 2.14, MAFFT >= 7.520, IQ-TREE >= 2.2

### Substeps

1. **16S extraction** -- Extracts the 16S rRNA gene from the GenBank annotation.

2. **BLAST search** -- Runs `blastn -remote` against the NCBI nt database with an Entrez query filter (`biomol_genomic[PROP] AND 16S[Title]`) to retrieve the 100 most similar 16S sequences. Command-line BLAST+ is used instead of the web API for reliability.

3. **Sequence and metadata retrieval** -- Fetches full GenBank records for each BLAST hit via the NCBI Entrez API. Extracts organism name, host, isolation source, strain, and country.

4. **Metadata annotation** -- Classifies each sequence as symbiont or non-symbiont based on the presence of a host field or symbiont-associated keywords. Parses host taxonomy (genus, family, order, phylum, and associated NCBI Taxonomy IDs) via Entrez. When the full host name does not match in NCBI Taxonomy (e.g., names with strain identifiers or author citations), the search falls back to genus + species, then genus alone.

5. **Alignment** -- Aligns all sequences (query + 100 hits) with MAFFT (`--auto`).

6. **Outgroup selection** -- Selects the sequence with the highest average pairwise distance to all others as the outgroup.

7. **Tree inference** -- Runs IQ-TREE with ModelFinder (`-m MFP`) testing both standard reversible (GTR) and non-reversible (UNREST) substitution models. UNREST models are included because endosymbiont genomes are known to exhibit compositional heterogeneity and strand-asymmetric substitution patterns, violating the assumptions of time-reversible models. The model search space includes empirical base frequencies (`-mfreq F`) and rate heterogeneity across sites (`-mrate E,I,G,R`). Branch support is assessed with 1000 ultrafast bootstrap replicates. Thread count is determined automatically (`-T AUTO`).

8. **Summary report** -- Generates `REPORT.md` with counts of symbionts, host taxonomy breakdown, alignment statistics, and the selected substitution model.

### Output

All output is written to `16s_phylogenetic_analysis/`:

| File | Description |
|------|-------------|
| `query_16s.fasta` | Extracted 16S rRNA from query genome |
| `blast_results.xml` | Raw BLAST output (XML format) |
| `blast_hits.fasta` | FASTA sequences for all BLAST hits |
| `blast_hits_raw_metadata.tsv` | Raw metadata from NCBI GenBank records |
| `sequence_metadata.tsv` | Annotated metadata with symbiont status and host taxonomy |
| `all_sequences.fasta` | Query + hits combined |
| `alignment.fasta` | MAFFT multiple sequence alignment |
| `16s_tree.treefile` | Maximum-likelihood tree (Newick format) |
| `16s_tree.iqtree` | Full IQ-TREE analysis report |
| `REPORT.md` | Summary report |

### Metadata columns

The annotated metadata table (`sequence_metadata.tsv`) contains:

| Column | Description |
|--------|-------------|
| `accession` | NCBI accession number |
| `taxon_name` | Cleaned organism name |
| `symbiont_status` | `symbiont` or `non-symbiont` |
| `host_name` | Host organism name |
| `host_taxon_id` | NCBI Taxonomy ID for the host |
| `host_genus` / `host_genus_id` | Host genus name and NCBI Taxonomy ID |
| `host_family` / `host_family_id` | Host family name and NCBI Taxonomy ID |
| `host_order` / `host_order_id` | Host order name and NCBI Taxonomy ID |
| `host_phylum` / `host_phylum_id` | Host phylum name and NCBI Taxonomy ID |
| `isolation_organ` | Tissue/organ of isolation (e.g., gut, bacteriome) |
| `annotation_source` | How host info was determined (`host_field`, `organism_name`, or `NA`) |

## Reproducibility

All analyses are driven by shell scripts and Python helpers that can be re-run from scratch. Each step that requires specific software defines its dependencies in a conda environment YAML file, and the corresponding shell script creates the environment automatically if it does not already exist.

An NCBI API key (environment variable `NCBI_API_KEY`) is recommended for steps that query NCBI Entrez, as it increases the rate limit from 3 to 10 requests per second.

Scripts were developed with the assistance of Claude Code (Anthropic).

## Software versions

| Software | Version | Purpose |
|----------|---------|---------|
| GTDB-Tk | 2.1.1 | Taxonomic classification |
| CheckM2 | (via conda) | Genome quality assessment |
| minimap2 | >= 2.26 | Read mapping (PacBio HiFi) |
| samtools | >= 1.17 | BAM processing and coverage |
| BLAST+ | >= 2.14 | Remote BLAST search |
| MAFFT | >= 7.520 | Multiple sequence alignment |
| IQ-TREE | >= 2.2 | Phylogenetic inference and model selection |
| BioPython | >= 1.81 | Sequence parsing and NCBI Entrez queries |
