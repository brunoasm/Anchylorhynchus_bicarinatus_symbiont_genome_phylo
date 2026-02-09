#!/usr/bin/env python3
"""
Merge CheckM statistics from multiple sources for phylogeny tip taxa.

This script combines:
1. CheckM results from Castelli et al. (2025) supplementary data
2. CheckM results generated in this study for missing taxa

Input files:
- castelli_et_al/Suppl_tab_7_newer_checkM_results.xlsx: CheckM data from Castelli et al.
- genome_phylogenetic_analysis/checkm_output/quality_report.tsv: Our CheckM results

Output:
- genome_phylogenetic_analysis/phylogeny_checkm_stats.csv: Merged stats for all phylogeny tips

Usage:
    python merge_checkm_stats.py
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


# Mapping from phylogeny tip names to search patterns in Castelli CheckM data
# Based on TARGET_TAXA in phylo_genome_helper.py
TIP_TO_CASTELLI_PATTERN = {
    "Hepatincola_Av": "Hepatincola_Av",
    "Hepatincola_Pdp": "Hepatincola_Pdp",
    "Hepatincola_Pp": "Hepatincola_Pp",
    "Tardigradibacter_bertolanii": "Tardigradibacter_bertolanii",
    "Symbiont_of_Haliotis": "GCA_002632265",
    "Marine_MAG_001510075": "GCA_001510075",
    "Marine_MAG_001830425": "GCA_001830425",
    "Marine_MAG_002327565": "GCA_002327565",
    "Marine_MAG_002687515": "GCA_002687515",
    "Marine_MAG_009694195": "GCA_009694195",
    "Marine_MAG_013204045": "GCA_013204045",
    "Marine_MAG_014859895": "GCA_014859895",
    "Marine_MAG_018662225": "GCA_018662225",
    "Outgroup_009649675": "GCA_009649675",
    "Symbiont_of_L_labralis": "GCA_009780035",
}

# Mapping from CheckM bin IDs to phylogeny tip names (for our CheckM results)
BIN_ID_TO_TIP = {
    "s7_ctg000008c": "s7_ctg000008c",
    "Terasakiella_pusilla_DSM_6293": "Terasakiella_pusilla_DSM_6293",
    "Thalassospira_profundimaris": "Thalassospira_profundimaris",
}


def load_castelli_checkm(xlsx_path: Path) -> pd.DataFrame:
    """Load CheckM data from Castelli et al. (2025) supplementary table."""
    df = pd.read_excel(xlsx_path)
    print(f"Loaded {len(df)} organisms from Castelli et al. CheckM data")
    return df


def load_our_checkm(tsv_path: Path) -> pd.DataFrame:
    """Load CheckM results from our analysis."""
    df = pd.read_csv(tsv_path, sep='\t')
    print(f"Loaded {len(df)} organisms from our CheckM results")
    return df


def match_castelli_taxa(castelli_df: pd.DataFrame) -> list[dict]:
    """Match phylogeny tips to Castelli et al. CheckM data."""
    results = []

    for tip_name, pattern in TIP_TO_CASTELLI_PATTERN.items():
        matches = castelli_df[
            castelli_df['Organism'].str.contains(pattern, case=False, na=False)
        ]

        if len(matches) >= 1:
            row = matches.iloc[0]
            results.append({
                'tip_name': tip_name,
                'completeness': row['Completeness'],
                'contamination': row['Contamination'],
                'gc_content': row['GC'],
                'genome_size': row['Genome size'],
                'n_contigs': row['# contigs'],
                'data_source': 'Castelli_et_al_2025_CheckM'
            })
            print(f"  {tip_name} -> {row['Organism']} (Castelli)")
        else:
            print(f"  WARNING: {tip_name} not found in Castelli data (pattern: {pattern})")

    return results


def add_our_checkm_taxa(our_df: pd.DataFrame) -> list[dict]:
    """Extract CheckM results for taxa we analyzed."""
    results = []

    for _, row in our_df.iterrows():
        bin_id = row['Bin Id']

        if bin_id not in BIN_ID_TO_TIP:
            continue

        tip_name = BIN_ID_TO_TIP[bin_id]

        results.append({
            'tip_name': tip_name,
            'completeness': row['Completeness'],
            'contamination': row['Contamination'],
            'gc_content': row['GC'] / 100,  # Convert % to fraction (like Castelli data)
            'genome_size': row['Genome size (bp)'],
            'n_contigs': row['# contigs'],
            'data_source': 'this_study_CheckM'
        })
        print(f"  {tip_name} (this study)")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Merge CheckM statistics for phylogeny tip taxa"
    )
    parser.add_argument(
        "--castelli-checkm",
        type=Path,
        default=Path("castelli_et_al/Suppl_tab_7_newer_checkM_results.xlsx"),
        help="Path to Castelli et al. CheckM xlsx file"
    )
    parser.add_argument(
        "--our-checkm",
        type=Path,
        default=Path("genome_phylogenetic_analysis/checkm_output/quality_report.tsv"),
        help="Path to our CheckM TSV output"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("genome_phylogenetic_analysis/phylogeny_checkm_stats.csv"),
        help="Output CSV path"
    )

    args = parser.parse_args()

    # Validate input files
    if not args.castelli_checkm.exists():
        print(f"ERROR: Castelli CheckM file not found: {args.castelli_checkm}")
        sys.exit(1)
    if not args.our_checkm.exists():
        print(f"ERROR: Our CheckM file not found: {args.our_checkm}")
        sys.exit(1)

    print("=" * 60)
    print("Merging CheckM statistics for phylogeny tips")
    print("=" * 60)
    print()

    # Load data
    castelli_df = load_castelli_checkm(args.castelli_checkm)
    our_df = load_our_checkm(args.our_checkm)
    print()

    # Match taxa
    print("Matching Castelli et al. data:")
    castelli_results = match_castelli_taxa(castelli_df)
    print()

    print("Adding this study's CheckM results:")
    our_results = add_our_checkm_taxa(our_df)
    print()

    # Combine and save
    all_results = castelli_results + our_results
    df = pd.DataFrame(all_results)
    df = df.sort_values('tip_name')

    print("=" * 60)
    print(f"Final results: {len(df)} taxa")
    print("=" * 60)
    print(df.to_string(index=False))
    print()

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(args.output, index=False)
    print(f"Saved to: {args.output}")


if __name__ == "__main__":
    main()
