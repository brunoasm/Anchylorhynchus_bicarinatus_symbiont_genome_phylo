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
    |-- [06] 16S rRNA extraction, multi-source search, alignment, phylogenetic inference
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

This genome represents a novel lineage within Alphaproteobacteria (order WRAU01, as classified by GTDB-Tk using RED placement). It consists of a single circular contig at 100% breadth of coverage. The genome was subsequently annotated externally (GenBank file: `s7_ctg000008c.bgpipe.output_838931.gb`).

## Step 06: 16S rRNA phylogenetic analysis

**Script:** `06_16s_phylogenetic_analysis.sh`
**Helper:** `phylo_16s_helper.py`
**Environment file:** `symbiont_phylo_env.yml`

Extracts the 16S rRNA gene from the annotated GenBank file, searches for related sequences from four independent sources, deduplicates, and builds a phylogenetic tree to place the organism among its closest relatives.

**Dependencies:** Python >= 3.10, BioPython >= 1.81, pandas >= 2.0, BLAST+ >= 2.14, MAFFT >= 7.520, IQ-TREE >= 2.2, cd-hit

### Data sources

The pipeline integrates 16S rRNA sequences from four independent sources:

1. **BLAST** -- Remote BLAST against the NCBI nt database with entrez query `biomol_genomic[PROP] AND 16S[Title]`, retrieving the 100 most similar sequences. Command-line BLAST+ is used instead of the web API for reliability.

2. **ARB SILVA** -- Nearest-neighbour sequences from a SILVA SSU search (minimum identity 0.7, 100 neighbours). The SILVA FASTA file is pre-downloaded from the [SILVA website](https://www.arb-silva.de/) and provided via the `SILVA_FASTA` environment variable. Accessions are parsed from SILVA `[nearest_slv=...]` headers. When the same genome has multiple 16S copies, only the longest region per accession is retained.

3. **GTDB WRAU01** -- 16S rRNA sequences from all genomes classified in the WRAU01 order within GTDB. Genome accessions are extracted from the GTDB-Tk MSA FASTA file by searching for WRAU01 in the taxonomic annotations. For each genome assembly, NCBI nucleotide records are searched and GenBank features are parsed for 16S rRNA annotations.

4. **Castelli et al. 2025** -- 16S rRNA sequences from three *Hepatincolaceae* genome assemblies reported by [Castelli et al. (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11724238/). That study showed that the phylogenetic grouping of Hepatincolaceae with Holosporales was artefactual due to compositional biases, and that the family instead represents an independent lineage within Rhodospirillales. The three genomes included here (GCF\_000688235.1, GCA\_001510075.1, GCA\_045504435.1) provide additional phylogenetic context for placing the WRAU01 lineage relative to described Hepatincolaceae members. For each assembly, 16S is recovered from GenBank feature annotations or, when annotations are absent, by local BLAST of the query 16S against the genome contigs.

### Deduplication

After merging sequences from all four sources, accession-level deduplication is performed: when the same accession appears in multiple sources, one copy of the sequence is kept and all contributing source labels are recorded (e.g., `blast,silva`). The merged set is then further clustered with **cd-hit-est** at 99% identity (`-c 0.99 -n 10 -aS 0.9`), keeping the longest representative per cluster.

### Reverse complement handling

MAFFT alignment is run with `--adjustdirection`, which detects and reverse-complements sequences that are in the wrong orientation. After alignment, the `_R_` prefix that MAFFT adds to reverse-complemented sequences is stripped, and the identities of flipped sequences are logged.

### Substeps

1. **Clean old results** -- Removes previous output files (alignment, tree, metadata, combined FASTAs) to avoid stale data. BLAST XML and query FASTA are preserved if they exist (expensive to regenerate).

2. **16S extraction** -- Extracts the 16S rRNA gene from the GenBank annotation.

3. **BLAST search** -- Runs `blastn -remote` against the NCBI nt database (skipped if `blast_results.xml` already exists).

4. **Retrieve BLAST sequences** -- Fetches full GenBank records for each BLAST hit via the NCBI Entrez API. Extracts organism name, host, isolation source, strain, and country.

5. **Parse SILVA file** -- Parses SILVA nearest-neighbour FASTA headers to extract accession coordinates and similarity scores.

6. **Retrieve SILVA sequences** -- Fetches 16S subsequences from GenBank for each SILVA hit, using `seq_start`/`seq_stop` parameters.

7. **Parse GTDB WRAU01** -- Searches the GTDB MSA FASTA for WRAU01-classified genomes and extracts NCBI assembly accessions.

8. **Retrieve GTDB 16S** -- For each WRAU01 genome, searches NCBI nucleotide records and parses GenBank features for 16S rRNA annotations.

8b. **Retrieve Castelli et al. 16S** -- Retrieves 16S rRNA sequences from three Hepatincolaceae genome assemblies (Castelli et al. 2025) using the same annotation scan and BLAST fallback approach as the GTDB step.

9. **Merge sources** -- Combines FASTA and metadata from all four sources plus the query, deduplicating by accession and tracking which sources contribute each sequence.

10. **Deduplicate with cd-hit-est** -- Clusters at 99% identity, keeping longest representatives.

11. **Metadata annotation** -- Classifies each sequence as symbiont or non-symbiont based on the presence of a host field or symbiont-associated keywords. Parses host taxonomy (genus, family, order, phylum, and associated NCBI Taxonomy IDs) via Entrez. Looks up microbial taxonomy for each organism and formats as a GTDB-style lineage string. Only sequences surviving deduplication are included.

12. **Alignment** -- Aligns all surviving sequences with MAFFT (`--auto --adjustdirection`). Reverse-complemented sequences are logged and the `_R_` prefix is cleaned.

13. **Tree inference** -- Runs IQ-TREE with ModelFinder in two separate passes (required by IQ-TREE 3, which does not allow mixing reversible and non-reversible models in a single run). The first pass tests the reversible GTR model (`-m MFP -mset GTR`). The second pass tests the non-reversible UNREST model plus all Lie Markov models (`-m MFP+LM -mset UNREST --nonrev-model`). Non-reversible and Lie Markov models are included because endosymbiont genomes are known to exhibit compositional heterogeneity and strand-asymmetric substitution patterns, violating the assumptions of time-reversible models. The Lie Markov models (Woodhams et al. 2015) provide a hierarchy of ~35 non-reversible models intermediate between GTR (6 rate parameters) and UNREST (12 rate parameters), allowing model selection to find the best trade-off between fit and complexity. No outgroup is specified for the non-reversible pass: these models infer the root position directly from substitution asymmetry in the data. The best model is selected by comparing BIC scores across both passes, and the corresponding tree files are copied to the final output names. Both model search spaces include empirical base frequencies (`-mfreq F`) and rate heterogeneity across sites (`-mrate E,I,G,R`). Branch support is assessed with 1000 ultrafast bootstrap replicates. Thread count is determined automatically (`-T AUTO`).

14. **Summary report** -- Generates `REPORT.md` with per-source sequence counts, deduplication statistics, reverse-complement counts, symbiont breakdown, host taxonomy, alignment statistics, and the selected substitution model.

### Output

All output is written to `16s_phylogenetic_analysis/`:

| File | Description |
|------|-------------|
| `query_16s.fasta` | Extracted 16S rRNA from query genome |
| `blast_results.xml` | Raw BLAST output (XML format) |
| `blast_hits.fasta` | FASTA sequences for BLAST hits |
| `blast_hits_raw_metadata.tsv` | Raw metadata from BLAST GenBank records |
| `silva_hits.tsv` | Parsed SILVA nearest-neighbour accessions |
| `silva_hits.fasta` | FASTA sequences for SILVA hits |
| `silva_hits_raw_metadata.tsv` | Raw metadata from SILVA GenBank records |
| `gtdb_wrau01_accessions.txt` | WRAU01 genome accessions extracted from GTDB |
| `gtdb_hits.fasta` | 16S sequences from GTDB WRAU01 genomes |
| `gtdb_hits_raw_metadata.tsv` | Raw metadata from GTDB GenBank records |
| `castelli_accessions.txt` | Castelli et al. 2025 genome accessions |
| `castelli_hits.fasta` | 16S sequences from Castelli et al. 2025 genomes |
| `castelli_hits_raw_metadata.tsv` | Raw metadata from Castelli GenBank records |
| `all_sequences.fasta` | Merged sequences from all sources + query |
| `all_sources_raw_metadata.tsv` | Merged raw metadata with source tracking |
| `deduplicated.fasta` | After cd-hit-est 99% clustering |
| `sequence_metadata.tsv` | Annotated metadata (deduplication survivors only) |
| `alignment_raw.fasta` | Raw MAFFT alignment (with `_R_` prefixes) |
| `alignment.fasta` | Cleaned MAFFT alignment |
| `16s_tree.treefile` | Maximum-likelihood tree from the best model (Newick format) |
| `16s_tree.iqtree` | IQ-TREE report for the best model |
| `16s_tree_rev.iqtree` | IQ-TREE report for reversible model search (GTR) |
| `16s_tree_nonrev.iqtree` | IQ-TREE report for non-reversible model search (UNREST + Lie Markov) |
| `REPORT.md` | Summary report |

### Metadata columns

The annotated metadata table (`sequence_metadata.tsv`) contains:

| Column | Description |
|--------|-------------|
| `accession` | NCBI accession number |
| `taxon_name` | Cleaned organism name |
| `source` | Data source(s): `blast`, `silva`, `gtdb`, `query`, or comma-separated combination |
| `symbiont_status` | `symbiont` or `non-symbiont` |
| `host_name` | Host organism name |
| `host_taxon_id` | NCBI Taxonomy ID for the host |
| `host_genus` / `host_genus_id` | Host genus name and NCBI Taxonomy ID |
| `host_family` / `host_family_id` | Host family name and NCBI Taxonomy ID |
| `host_order` / `host_order_id` | Host order name and NCBI Taxonomy ID |
| `host_phylum` / `host_phylum_id` | Host phylum name and NCBI Taxonomy ID |
| `isolation_organ` | Tissue/organ of isolation (e.g., gut, bacteriome) |
| `annotation_source` | How host info was determined (`host_field`, `organism_name`, or `NA`) |
| `microbial_lineage` | NCBI taxonomy of the microbe in GTDB-style format (`d__Name:ID;p__Name:ID;...`) |

## Conda environments

The pipeline uses three conda environments:

| Environment | Scripts | Environment file | Notes |
|-------------|---------|-----------------|-------|
| `gtdbtk-2.1.1` | 02 | (pre-existing) | Must be created manually before running step 02 |
| `checkm2` | 04 | (created by script) | Created automatically by `04_run_checkm2.sh` |
| `symbiont_phylo` | 01, 03, 05, 06 | `symbiont_phylo_env.yml` | Created automatically if absent |

## Reproducibility

All analyses are driven by shell scripts and Python helpers that can be re-run from scratch. Each step that requires specific software defines its dependencies in a conda environment YAML file, and the corresponding shell script creates the environment automatically if it does not already exist.

An NCBI API key (environment variable `NCBI_API_KEY`) is recommended for steps that query NCBI Entrez, as it increases the rate limit from 3 to 10 requests per second.

Step 06 accepts two optional environment variables:
- `SILVA_FASTA` -- Path to a SILVA SSU search result FASTA file (nearest neighbours, min identity 0.7, 100 neighbours). If not set, the SILVA source is skipped.
- `GTDB_MSA` -- Path to the GTDB-Tk MSA FASTA file (defaults to `gtdbtk_out/classify/gtdbtk.bac120.msa.fasta.gz`). If not found, the GTDB source is skipped.

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
| cd-hit | (via conda) | Sequence clustering and deduplication |
| BioPython | >= 1.81 | Sequence parsing and NCBI Entrez queries |
