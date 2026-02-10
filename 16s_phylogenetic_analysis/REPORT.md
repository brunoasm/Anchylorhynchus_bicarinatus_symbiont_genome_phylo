# 16S rRNA Phylogenetic Analysis Report

## Overview

- **Total reference sequences (merged)**: 148
  - from blast: 100
  - from castelli: 2
  - from gtdb: 5
  - from silva: 47

- **Before cd-hit deduplication**: 149 sequences
- **After cd-hit deduplication (99% identity)**: 54 sequences
- **Sequences removed**: 95

- **Symbionts identified**: 19
- **Non-symbionts**: 35
- **Sequences with known host**: 54

## Source Breakdown

| Source | Count |
|--------|-------|
| blast | 35 |
| silva | 9 |
| blast,silva | 4 |
| gtdb | 3 |
| castelli | 2 |
| query | 1 |

## Host Orders

| Order | Count |
|-------|-------|
| Fabales | 4 |
| Hymenoptera | 3 |
| Blattodea | 2 |
| Diplostraca | 1 |
| Priapulimorphida | 1 |
| Hemiptera | 1 |
| Ectocarpales | 1 |
| Asparagales | 1 |
| Scorpiones | 1 |
| Brassicales | 1 |
| Sapindales | 1 |
| Isopoda | 1 |
| Parachela | 1 |

- **Sequences reverse-complemented by MAFFT**: 2

## Alignment

- **Total sequences**: 54 (query + references)
- **Alignment length**: 6883 bp

## Phylogenetic Tree

- **Tree file**: `16s_tree.treefile`
- **Bootstrap replicates**: 1000 (ultrafast)

## Output Files

| File | Description |
|------|-------------|
| `query_16s.fasta` | Extracted 16S from query genome |
| `blast_results.xml` | Raw NCBI BLAST output |
| `blast_hits.fasta` | Retrieved BLAST sequences |
| `blast_hits_raw_metadata.tsv` | Raw BLAST metadata from NCBI |
| `silva_hits.tsv` | Parsed SILVA nearest-neighbour accessions |
| `silva_hits.fasta` | Retrieved SILVA sequences |
| `silva_hits_raw_metadata.tsv` | Raw SILVA metadata from NCBI |
| `gtdb_wrau01_accessions.txt` | WRAU01 genome accessions from GTDB |
| `gtdb_hits.fasta` | Retrieved GTDB 16S sequences |
| `gtdb_hits_raw_metadata.tsv` | Raw GTDB metadata from NCBI |
| `castelli_accessions.txt` | Castelli et al. 2025 genome accessions |
| `castelli_hits.fasta` | Retrieved Castelli et al. 2025 16S sequences |
| `castelli_hits_raw_metadata.tsv` | Raw Castelli metadata from NCBI |
| `all_sequences.fasta` | All sources merged + query |
| `all_sources_raw_metadata.tsv` | Merged raw metadata |
| `deduplicated.fasta` | After cd-hit-est 99% clustering |
| `sequence_metadata.tsv` | Annotated metadata (survivors only) |
| `alignment.fasta` | MAFFT alignment |
| `16s_tree.treefile` | Best ML tree (Newick format) |
| `16s_tree.iqtree` | IQ-TREE report for best model |
| `16s_tree_rev.iqtree` | IQ-TREE report (reversible models) |
| `16s_tree_nonrev.iqtree` | IQ-TREE report (non-reversible models) |
