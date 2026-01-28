#!/usr/bin/env python3
"""
Helper script for filtering genomes and generating README files.
Used by 05_filter_genomes_coverage.sh
"""

import argparse
import pandas as pd
from pathlib import Path
import sys


def load_checkm2_results(checkm2_path: str) -> pd.DataFrame:
    """Load CheckM2 quality report TSV."""
    return pd.read_csv(checkm2_path, sep='\t')


def load_gtdbtk_results(gtdbtk_path: str) -> pd.DataFrame:
    """Load GTDB-tk summary TSV."""
    return pd.read_csv(gtdbtk_path, sep='\t')


def convert_name_for_gtdbtk(checkm2_name: str) -> str:
    """
    Convert CheckM2 naming (s7_ctg000008c) to GTDB-tk naming (s7.ctg000008c).
    GTDB-tk uses dots instead of underscores between sample and contig IDs.
    """
    parts = checkm2_name.split('_')
    if len(parts) >= 2:
        return f"{parts[0]}.{'_'.join(parts[1:])}"
    return checkm2_name


def filter_genomes(checkm2_df: pd.DataFrame,
                   min_completeness: float = 90.0,
                   max_contamination: float = 5.0) -> list:
    """
    Filter genomes based on quality criteria.
    Returns list of genome names that pass the filter.
    """
    filtered = checkm2_df[
        (checkm2_df['Completeness'] >= min_completeness) &
        (checkm2_df['Contamination'] < max_contamination)
    ]
    return filtered['Name'].tolist()


def get_checkm2_stats(checkm2_df: pd.DataFrame, genome_name: str) -> dict:
    """Get CheckM2 statistics for a genome."""
    row = checkm2_df[checkm2_df['Name'] == genome_name].iloc[0]
    return {
        'completeness': row['Completeness'],
        'contamination': row['Contamination'],
        'genome_size': row['Genome_Size'],
        'n50': row['Contig_N50'],
        'gc_content': row['GC_Content'],
        'coding_density': row['Coding_Density'],
        'total_contigs': row['Total_Contigs'],
        'total_cds': row['Total_Coding_Sequences'],
    }


def get_gtdbtk_stats(gtdbtk_df: pd.DataFrame, genome_name: str) -> dict:
    """Get GTDB-tk statistics for a genome."""
    gtdbtk_name = convert_name_for_gtdbtk(genome_name)
    row = gtdbtk_df[gtdbtk_df['user_genome'] == gtdbtk_name]

    if row.empty:
        return {
            'classification': 'Not found in GTDB-tk results',
            'reference': 'N/A',
            'ani': 'N/A',
            'red_value': 'N/A',
            'method': 'N/A',
            'warnings': 'N/A',
        }

    row = row.iloc[0]
    classification = row['classification']
    reference = row['closest_genome_reference'] if pd.notna(row['closest_genome_reference']) else 'N/A'
    ani = row['closest_genome_ani'] if pd.notna(row['closest_genome_ani']) else 'N/A'
    red_value = row['red_value'] if pd.notna(row['red_value']) else 'N/A'
    method = row['classification_method'] if pd.notna(row['classification_method']) else 'N/A'
    warnings = row['warnings'] if pd.notna(row['warnings']) else 'None'

    return {
        'classification': classification,
        'reference': reference,
        'ani': ani,
        'red_value': red_value,
        'method': method,
        'warnings': warnings,
    }


def generate_readme(genome_name: str, checkm2_stats: dict, gtdbtk_stats: dict,
                    coverage_stats: dict = None) -> str:
    """Generate README content for a genome."""

    # Format classification for readability
    classification = gtdbtk_stats['classification']
    if classification and classification != 'N/A':
        tax_levels = classification.split(';')
        tax_formatted = '\n'.join([f"  - {level}" for level in tax_levels])
    else:
        tax_formatted = f"  - {classification}"

    readme = f"""# Genome: {genome_name}

## Taxonomy (GTDB-tk)

**Full classification:**
{tax_formatted}

- **Classification method**: {gtdbtk_stats['method']}
- **Closest genome reference**: {gtdbtk_stats['reference']}
- **ANI to reference**: {gtdbtk_stats['ani']}
- **RED value**: {gtdbtk_stats['red_value']}
- **Warnings**: {gtdbtk_stats['warnings']}

## Quality Metrics (CheckM2)

| Metric | Value |
|--------|-------|
| Completeness | {checkm2_stats['completeness']:.2f}% |
| Contamination | {checkm2_stats['contamination']:.2f}% |
| Genome size | {checkm2_stats['genome_size']:,} bp |
| N50 | {checkm2_stats['n50']:,} bp |
| GC content | {checkm2_stats['gc_content'] * 100:.2f}% |
| Coding density | {checkm2_stats['coding_density']:.3f} |
| Total contigs | {checkm2_stats['total_contigs']} |
| Total CDS | {checkm2_stats['total_cds']} |
"""

    if coverage_stats:
        readme += f"""
## Coverage Statistics

| Metric | Value |
|--------|-------|
| Mean depth | {coverage_stats['meandepth']:.2f}x |
| Breadth of coverage | {coverage_stats['coverage']:.2f}% |
| Mapped reads | {coverage_stats['numreads']:,} |
| Total bases covered | {coverage_stats['covbases']:,} bp |
"""

    return readme


def main():
    parser = argparse.ArgumentParser(description='Filter genomes and generate README files')
    subparsers = parser.add_subparsers(dest='command', required=True)

    # Filter command
    filter_parser = subparsers.add_parser('filter', help='Filter genomes by quality')
    filter_parser.add_argument('--checkm2', required=True, help='Path to CheckM2 quality_report.tsv')
    filter_parser.add_argument('--min-completeness', type=float, default=90.0)
    filter_parser.add_argument('--max-contamination', type=float, default=5.0)

    # Readme command
    readme_parser = subparsers.add_parser('readme', help='Generate README for a genome')
    readme_parser.add_argument('--genome', required=True, help='Genome name')
    readme_parser.add_argument('--checkm2', required=True, help='Path to CheckM2 quality_report.tsv')
    readme_parser.add_argument('--gtdbtk', required=True, help='Path to GTDB-tk summary TSV')
    readme_parser.add_argument('--coverage', help='Path to samtools coverage output')
    readme_parser.add_argument('--output', required=True, help='Output README path')

    args = parser.parse_args()

    if args.command == 'filter':
        checkm2_df = load_checkm2_results(args.checkm2)
        passing_genomes = filter_genomes(
            checkm2_df,
            min_completeness=args.min_completeness,
            max_contamination=args.max_contamination
        )
        for genome in passing_genomes:
            print(genome)

    elif args.command == 'readme':
        checkm2_df = load_checkm2_results(args.checkm2)
        gtdbtk_df = load_gtdbtk_results(args.gtdbtk)

        checkm2_stats = get_checkm2_stats(checkm2_df, args.genome)
        gtdbtk_stats = get_gtdbtk_stats(gtdbtk_df, args.genome)

        coverage_stats = None
        if args.coverage and Path(args.coverage).exists():
            cov_df = pd.read_csv(args.coverage, sep='\t')
            if not cov_df.empty:
                coverage_stats = {
                    'meandepth': cov_df['meandepth'].iloc[0],
                    'coverage': cov_df['coverage'].iloc[0],
                    'numreads': int(cov_df['numreads'].iloc[0]),
                    'covbases': int(cov_df['covbases'].iloc[0]),
                }

        readme_content = generate_readme(args.genome, checkm2_stats, gtdbtk_stats, coverage_stats)

        with open(args.output, 'w') as f:
            f.write(readme_content)

        print(f"README generated: {args.output}")


if __name__ == '__main__':
    main()
