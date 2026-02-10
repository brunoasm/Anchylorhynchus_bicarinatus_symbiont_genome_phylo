# Genome-based Phylogenetic Analysis Report

## Overview

This phylogeny was built using the framework of Castelli et al. (2025), enriched
with all GTDB genomes classified in order WRAU01.

- **Original taxa (step 06)**: 18 (Castelli targets + outgroups + *A. bicarinatus* symbiont)
- **GTDB WRAU01 genomes found**: 8 (1 existing, 1 duplicate, 6 new)
- **Final taxa**: 24

## Ortholog Groups

- **OGs mapped from eggNOG to Castelli**: 169 / 179 (primary MAG)
- **Additional OGs from GTDB genomes**: 6
- **OGs aligned**: 175
- **OGs after trimming**: 175

## Concatenated Alignment

- **Taxa**: 24
- **Total alignment length**: 64,496 AA sites

## After Compositional Bias Removal

- **Alignment length**: 45,147 AA sites
- **Sites retained**: 70.0%

## Phylogenetic Tree

- **Tree file**: `genome_tree.treefile`
- **Model**: LG+I+R5 (selected by ModelFinder, BIC)
- **Bootstrap replicates**: 1000 (ultrafast) + SH-aLRT
- **Best log-likelihood**: -1,399,360.213

## CheckM v1 Statistics

All genome quality estimates use CheckM v1 (`lineage_wf`) for comparability
with Castelli et al. (2025). Sources:

- **Castelli et al. (2025)**: 15 taxa (from Suppl. Table 7)
- **This study**: 9 taxa (3 original + 6 GTDB enrichment genomes)

## GTDB Enrichment Genomes

| Accession | Tip Name | Comp% | Cont% | GC | Size (bp) | Contigs | Host |
|-----------|----------|-------|-------|----|-----------|---------|------|
| GCA_031256515.1 | WRAU01_MAG_031256515 | 98.1 | 0.0 | 0.294 | 1,487,858 | 104 | *Cornitermes pugnax* (Termitidae) |
| GCA_031274895.1 | WRAU01_MAG_031274895 | 97.3 | 0.0 | 0.360 | 1,499,021 | 137 | *Jugositermes tuberculatus* (Termitidae) |
| GCA_945888065.1 | WRAU01_MAG_945888065 | 86.8 | 0.0 | 0.254 | 867,558 | 26 | — |
| GCA_947460385.1 | WRAU01_MAG_947460385 | 88.7 | 1.6 | 0.296 | 1,093,185 | 66 | — |
| GCA_964020055.1 | Symbiont_of_P_viduata | 97.8 | 0.0 | 0.268 | 1,105,045 | 1 | *Pipizella viduata* (Syrphidae) |
| GCA_964020105.1 | Symbiont_of_E_torrentis | 98.9 | 0.0 | 0.312 | 1,383,406 | 1 | *Ecdyonurus torrentis* (Heptageniidae) |

## Taxa in Final Tree

| Taxon | Group |
|-------|-------|
| Hepatincola_Av | Hepatincolaceae |
| Hepatincola_Pdp | Hepatincolaceae |
| Hepatincola_Pp | Hepatincolaceae |
| Marine_MAG_001510075 | Thalassospiraceae |
| Marine_MAG_001830425 | Thalassospiraceae |
| Marine_MAG_002327565 | Thalassospiraceae |
| Marine_MAG_002687515 | Thalassospiraceae |
| Marine_MAG_009694195 | Thalassospiraceae |
| Marine_MAG_013204045 | Thalassospiraceae |
| Marine_MAG_014859895 | Thalassospiraceae |
| Marine_MAG_018662225 | Thalassospiraceae |
| Outgroup_009649675 | Outgroup |
| Symbiont_of_E_torrentis | Thalassospiraceae |
| Symbiont_of_Haliotis | Thalassospiraceae |
| Symbiont_of_L_labralis | Hepatincolaceae |
| Symbiont_of_P_viduata | Thalassospiraceae |
| Tardigradibacter_bertolanii | Hepatincolaceae |
| Terasakiella_pusilla_DSM_6293 | Thalassospiraceae |
| Thalassospira_profundimaris | Outgroup |
| WRAU01_MAG_031256515 | Thalassospiraceae |
| WRAU01_MAG_031274895 | Thalassospiraceae |
| WRAU01_MAG_945888065 | Thalassospiraceae |
| WRAU01_MAG_947460385 | Thalassospiraceae |
| s7_ctg000008c | Hepatincolaceae |
