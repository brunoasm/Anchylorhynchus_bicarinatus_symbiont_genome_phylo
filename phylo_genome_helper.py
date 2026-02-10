#!/usr/bin/env python3
"""
Helper script for genome-based phylogenetic analysis.
Used by 07_genome_phylogenetic_analysis.sh

Subcommands:
  map-ogs            Parse eggNOG-mapper output, match OG IDs to Castelli filenames
  extract-and-align  Filter taxa from OG files, add new protein, run MAFFT L-INS-i
  trim-alignments    Run BMGE on each OG alignment
  concatenate        Build supermatrix from trimmed OGs + partition file
  remove-comp-bias   Remove compositionally biased sites (chi-squared method)
  summary            Generate markdown report
"""

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats


# ---------------------------------------------------------------------------
# Target taxa definitions
# ---------------------------------------------------------------------------

# Map from a lookup key to a clean display name.
# Keys are either exact header matches or GCA 9-digit accession numbers.
TARGET_TAXA = {
    # Hepatincolaceae (WRAU01)
    "Tardigrade.Symbiont": "Tardigradibacter_bertolanii",
    "Hepatincola_Pp_genome.fasta": "Hepatincola_Pp",
    "Hepatincola_Av_genome.fasta": "Hepatincola_Av",
    "Hepatincola_Pdp_scaffold.fasta": "Hepatincola_Pdp",
    "009780035": "Symbiont_of_L_labralis",
    # Thalassospiraceae
    "000688235": "Terasakiella_pusilla_DSM_6293",
    "002632265": "Symbiont_of_Haliotis",
    "001510075": "Marine_MAG_001510075",
    "001830425": "Marine_MAG_001830425",
    "002327565": "Marine_MAG_002327565",
    "002687515": "Marine_MAG_002687515",
    "009694195": "Marine_MAG_009694195",
    "013204045": "Marine_MAG_013204045",
    "014859895": "Marine_MAG_014859895",
    "018662225": "Marine_MAG_018662225",
    # Outgroup
    "009649675": "Outgroup_009649675",
    "000300275": "Thalassospira_profundimaris",
}

# Regex to extract 9-digit GCA accession from any header format
GCA_PATTERN = re.compile(r"GCA[._](\d{9})")


def _match_header(header: str) -> str | None:
    """Return the display name if *header* matches a target taxon, else None.

    Matching rules:
    1. Exact match on the full header string.
    2. Extract GCA 9-digit accession via regex and look up.
    """
    header = header.strip()
    if header in TARGET_TAXA:
        return TARGET_TAXA[header]
    m = GCA_PATTERN.search(header)
    if m and m.group(1) in TARGET_TAXA:
        return TARGET_TAXA[m.group(1)]
    return None


# ---------------------------------------------------------------------------
# map-ogs
# ---------------------------------------------------------------------------

def map_ogs(eggnog_annotations: str, castelli_og_dir: str, output: str) -> None:
    """Parse eggNOG-mapper annotations and match OG IDs to Castelli OG files.

    For each protein in the eggNOG output, extract its OG assignments from
    the ``eggNOG_OGs`` column.  Match those OG IDs against the filename
    prefixes of the Castelli single-copy OG files (the part before the first
    ``_aggiunta``).

    When a protein maps to multiple OGs that all correspond to the same
    Castelli file, keep the one with the best eggNOG score.  When multiple
    proteins map to the same OG (paralogs), keep the one with the best score.
    """
    # Build lookup: OG prefix -> Castelli filename
    og_dir = Path(castelli_og_dir)
    og_files = sorted(og_dir.glob("*.faa"))
    prefix_to_file = {}
    for f in og_files:
        prefix = f.name.split("_aggiunta")[0]
        prefix_to_file[prefix] = f.name

    print(f"Loaded {len(prefix_to_file)} Castelli OG file prefixes")

    # Parse eggNOG annotations
    mappings = []  # list of (protein_id, og_prefix, castelli_file, score)
    with open(eggnog_annotations) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            protein_id = fields[0]
            score = float(fields[3]) if fields[3] else 0.0
            eggnog_ogs = fields[4] if len(fields) > 4 else ""

            # Parse OG IDs: format like "COG0001@1|root,2TQKA@1224|Proteobacteria,..."
            for og_entry in eggnog_ogs.split(","):
                og_entry = og_entry.strip()
                if not og_entry:
                    continue
                og_id = og_entry.split("@")[0]
                if og_id in prefix_to_file:
                    mappings.append((protein_id, og_id, prefix_to_file[og_id], score))

    print(f"Found {len(mappings)} raw protein-to-OG mappings")

    # Deduplicate: keep best-scoring protein per OG
    best_per_og = {}
    for protein_id, og_id, castelli_file, score in mappings:
        if og_id not in best_per_og or score > best_per_og[og_id][3]:
            best_per_og[og_id] = (protein_id, og_id, castelli_file, score)

    print(f"After paralog filtering: {len(best_per_og)} OGs with a mapped protein")

    # Write output
    with open(output, "w") as fh:
        fh.write("protein_id\tog_id\tcastelli_file\tscore\n")
        for og_id in sorted(best_per_og):
            protein_id, _, castelli_file, score = best_per_og[og_id]
            fh.write(f"{protein_id}\t{og_id}\t{castelli_file}\t{score}\n")

    print(f"OG mapping written to: {output}")


# ---------------------------------------------------------------------------
# extract-and-align
# ---------------------------------------------------------------------------

def _read_fasta_dict(path: str) -> dict[str, str]:
    """Read a FASTA file into a dict {header: sequence}."""
    seqs = {}
    current_header = None
    current_seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    seqs[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_header is not None:
        seqs[current_header] = "".join(current_seq)
    return seqs


def _degap(seq: str) -> str:
    """Remove gap characters from a sequence."""
    return seq.replace("-", "").replace(".", "")


def _write_fasta(records: dict[str, str], path: str) -> None:
    """Write a dict {name: sequence} to a FASTA file."""
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def _load_extra_genomes(manifest_path: str) -> list[dict]:
    """Load extra genome data from a GTDB manifest JSON.

    Returns a list of dicts with keys: tip_name, proteins (dict), og_entries (list).
    """
    with open(manifest_path) as fh:
        manifest = json.load(fh)

    extras = []
    for genome in manifest["genomes"]:
        if genome["status"] != "new":
            continue
        tip_name = genome.get("tip_name")
        protein_path = genome.get("protein_fasta")
        og_mapping_path = genome.get("og_mapping")

        if not tip_name or not protein_path or not og_mapping_path:
            print(f"  Skipping {genome['accession']}: missing data")
            continue

        if not os.path.exists(protein_path) or not os.path.exists(og_mapping_path):
            print(f"  Skipping {genome['accession']}: files not found")
            continue

        proteins = _read_fasta_dict(protein_path)

        og_entries = []
        with open(og_mapping_path) as fh:
            fh.readline()  # skip header
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    og_entries.append({
                        "protein_id": fields[0],
                        "og_id": fields[1],
                        "castelli_file": fields[2],
                    })

        extras.append({
            "tip_name": tip_name,
            "proteins": proteins,
            "og_entries": og_entries,
            "og_lookup": {e["og_id"]: e for e in og_entries},
        })
        print(f"  Loaded extra genome: {tip_name} ({len(proteins)} proteins, {len(og_entries)} OGs)")

    return extras


def extract_and_align(og_mapping: str, castelli_og_dir: str, protein_file: str,
                      output_dir: str, threads: int = 4,
                      extra_genomes_manifest: str | None = None) -> None:
    """For each mapped OG: extract target taxa, add new MAG protein, align.

    If extra_genomes_manifest is provided, also adds proteins from GTDB genomes
    listed in the manifest.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Read MAG proteins
    mag_proteins = _read_fasta_dict(protein_file)
    print(f"Loaded {len(mag_proteins)} proteins from MAG")

    # Load extra genomes if provided
    extra_genomes = []
    if extra_genomes_manifest:
        print(f"\nLoading extra genomes from manifest:")
        extra_genomes = _load_extra_genomes(extra_genomes_manifest)
        print(f"  Total extra genomes: {len(extra_genomes)}\n")

    # Read primary OG mapping
    og_entries = []
    with open(og_mapping) as fh:
        header = fh.readline()  # skip header
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                og_entries.append({
                    "protein_id": fields[0],
                    "og_id": fields[1],
                    "castelli_file": fields[2],
                })

    # Build union of all OG IDs (primary + extra genomes)
    all_og_ids = {e["og_id"]: e["castelli_file"] for e in og_entries}
    primary_og_lookup = {e["og_id"]: e for e in og_entries}

    for extra in extra_genomes:
        for e in extra["og_entries"]:
            if e["og_id"] not in all_og_ids:
                all_og_ids[e["og_id"]] = e["castelli_file"]

    # Build unified OG list sorted by ID
    unified_entries = []
    for og_id in sorted(all_og_ids.keys()):
        unified_entries.append({
            "og_id": og_id,
            "castelli_file": all_og_ids[og_id],
        })

    print(f"Processing {len(unified_entries)} OGs ({len(og_entries)} from primary MAG, "
          f"{len(unified_entries) - len(og_entries)} additional from extra genomes)...")

    # Catalog all unique headers across OG files for S. maritima detection
    all_headers_count = defaultdict(int)

    n_aligned = 0
    n_skipped_few_taxa = 0
    n_skipped_no_proteins = 0

    for entry in unified_entries:
        og_id = entry["og_id"]
        castelli_file = entry["castelli_file"]
        og_path = Path(castelli_og_dir) / castelli_file

        if not og_path.exists():
            print(f"  [{og_id}] WARNING: Castelli file not found: {castelli_file}")
            continue

        # Read Castelli OG sequences
        og_seqs = _read_fasta_dict(str(og_path))

        # Track headers for S. maritima detection
        for hdr in og_seqs:
            all_headers_count[hdr] += 1

        # Filter to target taxa
        filtered = {}
        for hdr, seq in og_seqs.items():
            display_name = _match_header(hdr)
            if display_name is not None:
                filtered[display_name] = _degap(seq)

        # Add primary MAG protein if mapped to this OG
        has_any_new = False
        primary_entry = primary_og_lookup.get(og_id)
        if primary_entry:
            protein_id = primary_entry["protein_id"]
            if protein_id in mag_proteins:
                filtered["s7_ctg000008c"] = _degap(mag_proteins[protein_id])
                has_any_new = True

        # Add extra genome proteins if mapped to this OG
        for extra in extra_genomes:
            extra_entry = extra["og_lookup"].get(og_id)
            if extra_entry:
                protein_id = extra_entry["protein_id"]
                if protein_id in extra["proteins"]:
                    filtered[extra["tip_name"]] = _degap(extra["proteins"][protein_id])
                    has_any_new = True

        # Must have at least one new genome protein for this OG
        if not has_any_new:
            n_skipped_no_proteins += 1
            continue

        # Check minimum ingroup taxa (excluding outgroup)
        ingroup_names = [n for n in filtered if n not in
                         ("Outgroup_009649675", "Thalassospira_profundimaris")]
        if len(ingroup_names) < 4:
            n_skipped_few_taxa += 1
            continue

        # Write unaligned FASTA
        unaligned_path = Path(output_dir) / f"{og_id}_unaligned.fasta"
        _write_fasta(filtered, str(unaligned_path))

        # Align with MAFFT L-INS-i
        aligned_path = Path(output_dir) / f"{og_id}.fasta"
        try:
            result = subprocess.run(
                ["mafft", "--localpair", "--maxiterate", "1000",
                 "--thread", str(threads), "--quiet", str(unaligned_path)],
                capture_output=True, text=True, check=True,
            )
            with open(aligned_path, "w") as fh:
                fh.write(result.stdout)
            n_aligned += 1
            n_taxa = len(filtered)
            aln_len = len(result.stdout.split("\n")[1]) if "\n" in result.stdout else 0
            print(f"  [{og_id}] {n_taxa} taxa, aligned")
        except subprocess.CalledProcessError as e:
            print(f"  [{og_id}] MAFFT failed: {e.stderr[:200]}")
            continue

        # Clean up unaligned file
        unaligned_path.unlink(missing_ok=True)

    print(f"\nAlignment summary:")
    print(f"  OGs aligned: {n_aligned}")
    print(f"  Skipped (< 4 ingroup taxa): {n_skipped_few_taxa}")
    print(f"  Skipped (no new protein found): {n_skipped_no_proteins}")

    # S. maritima symbiont detection: find headers not matching any target
    # that appear in most OG files
    matched_headers = set()
    for hdr in all_headers_count:
        if _match_header(hdr) is not None:
            matched_headers.add(hdr)

    unmatched_frequent = []
    threshold = len(unified_entries) * 0.5  # present in >50% of OGs
    for hdr, count in sorted(all_headers_count.items(), key=lambda x: -x[1]):
        if hdr not in matched_headers and count > threshold:
            unmatched_frequent.append((hdr, count))

    if unmatched_frequent:
        print(f"\nPotential missed taxa (unmatched headers in >50% of OGs):")
        for hdr, count in unmatched_frequent[:10]:
            print(f"  {hdr} ({count}/{len(unified_entries)} OGs)")
        print("  Consider adding these to TARGET_TAXA if they are relevant.")


# ---------------------------------------------------------------------------
# trim-alignments
# ---------------------------------------------------------------------------

def trim_alignments(input_dir: str, output_dir: str) -> None:
    """Run BMGE on each aligned OG to trim poorly-aligned regions."""
    os.makedirs(output_dir, exist_ok=True)

    input_path = Path(input_dir)
    alignments = sorted(input_path.glob("*.fasta"))
    # Exclude unaligned files
    alignments = [a for a in alignments if "_unaligned" not in a.name]

    n_trimmed = 0
    n_empty = 0

    for aln_file in alignments:
        og_id = aln_file.stem
        out_file = Path(output_dir) / f"{og_id}.fasta"

        try:
            result = subprocess.run(
                ["bmge", "-i", str(aln_file), "-t", "AA", "-m", "BLOSUM30",
                 "-of", str(out_file)],
                capture_output=True, text=True, check=True,
            )

            # Check if output is non-empty
            if out_file.exists() and out_file.stat().st_size > 0:
                # Verify at least one sequence has non-gap characters
                seqs = _read_fasta_dict(str(out_file))
                has_content = any(_degap(s) for s in seqs.values())
                if has_content:
                    n_trimmed += 1
                    print(f"  [{og_id}] Trimmed OK")
                else:
                    out_file.unlink()
                    n_empty += 1
                    print(f"  [{og_id}] Empty after trimming, skipped")
            else:
                n_empty += 1
                print(f"  [{og_id}] Empty after trimming, skipped")

        except subprocess.CalledProcessError as e:
            print(f"  [{og_id}] BMGE failed: {e.stderr[:200]}")
            n_empty += 1

    print(f"\nTrimming summary:")
    print(f"  OGs trimmed: {n_trimmed}")
    print(f"  OGs empty after trimming: {n_empty}")


# ---------------------------------------------------------------------------
# concatenate
# ---------------------------------------------------------------------------

def concatenate(input_dir: str, output: str, partition_file: str) -> None:
    """Build a supermatrix from trimmed OG alignments + RAxML partition file.

    Missing taxa in individual OGs are filled with gap characters.
    """
    input_path = Path(input_dir)
    trimmed_files = sorted(input_path.glob("*.fasta"))

    if not trimmed_files:
        print("ERROR: No trimmed alignment files found", file=sys.stderr)
        sys.exit(1)

    # Collect all taxa across all OGs
    all_taxa = set()
    og_data = []  # list of (og_id, {taxon: aligned_seq}, aln_length)

    for f in trimmed_files:
        og_id = f.stem
        seqs = _read_fasta_dict(str(f))
        if not seqs:
            continue
        aln_len = len(next(iter(seqs.values())))
        all_taxa.update(seqs.keys())
        og_data.append((og_id, seqs, aln_len))

    all_taxa = sorted(all_taxa)
    print(f"Taxa: {len(all_taxa)}")
    print(f"OGs:  {len(og_data)}")

    # Build concatenated alignment
    concatenated = {taxon: "" for taxon in all_taxa}
    partitions = []
    current_pos = 1

    for og_id, seqs, aln_len in og_data:
        for taxon in all_taxa:
            if taxon in seqs:
                concatenated[taxon] += seqs[taxon]
            else:
                concatenated[taxon] += "-" * aln_len
        end_pos = current_pos + aln_len - 1
        partitions.append(f"AUTO, {og_id} = {current_pos}-{end_pos}")
        current_pos = end_pos + 1

    total_len = len(next(iter(concatenated.values())))
    print(f"Total alignment length: {total_len} sites")

    # Write concatenated FASTA
    _write_fasta(concatenated, output)
    print(f"Concatenated alignment: {output}")

    # Write partition file
    with open(partition_file, "w") as fh:
        for p in partitions:
            fh.write(p + "\n")
    print(f"Partition file: {partition_file}")

    # Per-taxon presence stats
    print(f"\nPer-taxon OG coverage:")
    for taxon in all_taxa:
        n_present = sum(1 for _, seqs, _ in og_data if taxon in seqs)
        print(f"  {taxon}: {n_present}/{len(og_data)} OGs")


# ---------------------------------------------------------------------------
# remove-comp-bias
# ---------------------------------------------------------------------------

# Standard amino acid alphabet (20 canonical AAs)
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def remove_comp_bias(alignment: str, threshold: float, output: str) -> None:
    """Remove the most compositionally biased sites from a concatenated alignment.

    Uses a per-site chi-squared test (Munoz-Gomez et al. 2019):
    1. For each column, count amino acid frequencies (ignoring gaps).
    2. Compute global expected frequencies across all non-gap residues.
    3. For each column: chi2 = sum_aa [(obs - exp)^2 / exp]
       where exp = global_freq_aa * n_valid_taxa_at_column.
    4. Rank columns by chi2 (descending), remove top *threshold* fraction.
    """
    seqs = _read_fasta_dict(alignment)
    if not seqs:
        print("ERROR: Empty alignment", file=sys.stderr)
        sys.exit(1)

    taxa = sorted(seqs.keys())
    aln_len = len(next(iter(seqs.values())))
    n_taxa = len(taxa)

    print(f"Input: {n_taxa} taxa, {aln_len} sites")
    print(f"Threshold: remove top {threshold * 100:.0f}% most biased sites")

    # Build sequence matrix
    aa_to_idx = {aa: i for i, aa in enumerate(AA_ALPHABET)}
    n_aa = len(AA_ALPHABET)

    # Count global AA frequencies
    global_counts = np.zeros(n_aa, dtype=np.float64)
    # Per-site counts: shape (aln_len, n_aa)
    site_counts = np.zeros((aln_len, n_aa), dtype=np.float64)
    site_valid = np.zeros(aln_len, dtype=np.float64)  # non-gap count per site

    for taxon in taxa:
        seq = seqs[taxon].upper()
        for j, aa in enumerate(seq):
            if aa in aa_to_idx:
                idx = aa_to_idx[aa]
                site_counts[j, idx] += 1
                global_counts[idx] += 1
                site_valid[j] += 1

    total_residues = global_counts.sum()
    if total_residues == 0:
        print("ERROR: No valid amino acid residues found", file=sys.stderr)
        sys.exit(1)

    global_freq = global_counts / total_residues

    # Chi-squared per site
    chi2_values = np.zeros(aln_len, dtype=np.float64)
    for j in range(aln_len):
        n_valid = site_valid[j]
        if n_valid < 2:
            chi2_values[j] = 0.0
            continue
        expected = global_freq * n_valid
        # Only compute for AAs with non-zero expected
        mask = expected > 0
        if mask.sum() == 0:
            chi2_values[j] = 0.0
            continue
        obs = site_counts[j, mask]
        exp = expected[mask]
        chi2_values[j] = np.sum((obs - exp) ** 2 / exp)

    # Rank and remove
    n_remove = int(np.ceil(aln_len * threshold))
    n_keep = aln_len - n_remove

    # Get indices of sites to keep (lowest chi2)
    sorted_indices = np.argsort(chi2_values)
    keep_indices = sorted(sorted_indices[:n_keep])

    print(f"Sites to remove: {n_remove}")
    print(f"Sites to keep: {n_keep}")
    print(f"Chi2 range: {chi2_values.min():.2f} - {chi2_values.max():.2f}")
    print(f"Chi2 threshold: {chi2_values[sorted_indices[n_keep - 1]]:.2f}")

    # Build filtered alignment
    filtered = {}
    for taxon in taxa:
        seq = seqs[taxon]
        filtered[taxon] = "".join(seq[j] for j in keep_indices)

    _write_fasta(filtered, output)
    print(f"Debiased alignment: {output}")


# ---------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------

def generate_summary(output_dir: str) -> None:
    """Generate a markdown summary report for the genome phylogenetic analysis."""
    out = Path(output_dir)

    report = [
        "# Genome-based Phylogenetic Analysis Report",
        "",
        "## Overview",
        "",
    ]

    # OG mapping stats
    og_mapping = out / "og_mapping.tsv"
    if og_mapping.exists():
        n_mapped = sum(1 for i, line in enumerate(open(og_mapping)) if i > 0 and line.strip())
        report.append(f"- **OGs mapped from eggNOG to Castelli**: {n_mapped} / 179")
        report.append("")

    # Aligned OGs
    aligned_dir = out / "aligned_ogs"
    if aligned_dir.exists():
        aligned_files = list(aligned_dir.glob("*.fasta"))
        aligned_files = [f for f in aligned_files if "_unaligned" not in f.name]
        report.append(f"- **OGs aligned**: {len(aligned_files)}")

    # Trimmed OGs
    trimmed_dir = out / "trimmed_ogs"
    if trimmed_dir.exists():
        trimmed_files = list(trimmed_dir.glob("*.fasta"))
        report.append(f"- **OGs after trimming**: {len(trimmed_files)}")

    report.append("")

    # Concatenated alignment stats
    concat_file = out / "concatenated.fasta"
    if concat_file.exists():
        seqs = _read_fasta_dict(str(concat_file))
        n_taxa = len(seqs)
        aln_len = len(next(iter(seqs.values()))) if seqs else 0
        report.extend([
            "## Concatenated Alignment",
            "",
            f"- **Taxa**: {n_taxa}",
            f"- **Total alignment length**: {aln_len} AA sites",
            "",
        ])

    # Debiased alignment stats
    debiased_file = out / "concatenated_debiased.fasta"
    if debiased_file.exists():
        seqs = _read_fasta_dict(str(debiased_file))
        aln_len = len(next(iter(seqs.values()))) if seqs else 0
        report.extend([
            "## After Compositional Bias Removal",
            "",
            f"- **Alignment length**: {aln_len} AA sites",
            "",
        ])

        if concat_file.exists():
            orig_seqs = _read_fasta_dict(str(concat_file))
            orig_len = len(next(iter(orig_seqs.values()))) if orig_seqs else 0
            if orig_len > 0:
                pct_kept = aln_len / orig_len * 100
                report.append(f"- **Sites retained**: {pct_kept:.1f}%")
                report.append("")

    # IQ-TREE results
    tree_file = out / "genome_tree.treefile"
    iqtree_file = out / "genome_tree.iqtree"
    if tree_file.exists():
        report.extend([
            "## Phylogenetic Tree",
            "",
            "- **Tree file**: `genome_tree.treefile`",
        ])

        if iqtree_file.exists():
            with open(iqtree_file) as fh:
                for line in fh:
                    if "Best-fit model:" in line:
                        model = line.split(":")[-1].strip()
                        report.append(f"- **Best-fit model**: {model}")
                        break

        report.extend([
            "- **Bootstrap replicates**: 1000 (ultrafast) + SH-aLRT",
            "",
        ])

    # Taxa list
    if debiased_file.exists() or concat_file.exists():
        src = debiased_file if debiased_file.exists() else concat_file
        seqs = _read_fasta_dict(str(src))
        if seqs:
            report.extend([
                "## Taxa in Final Tree",
                "",
                "| Taxon | Group |",
                "|-------|-------|",
            ])
            hepatincola_keys = {
                "Tardigradibacter_bertolanii", "Hepatincola_Pp",
                "Hepatincola_Av", "Hepatincola_Pdp",
                "Symbiont_of_L_labralis", "s7_ctg000008c",
            }
            outgroup_keys = {"Outgroup_009649675", "Thalassospira_profundimaris"}

            for taxon in sorted(seqs.keys()):
                if taxon in hepatincola_keys:
                    group = "Hepatincolaceae"
                elif taxon in outgroup_keys:
                    group = "Outgroup"
                else:
                    group = "Thalassospiraceae"
                report.append(f"| {taxon} | {group} |")
            report.append("")

    # Write report
    report_content = "\n".join(report)
    report_path = out / "REPORT.md"
    with open(report_path, "w") as fh:
        fh.write(report_content)

    print(f"Report generated: {report_path}")
    print()
    print(report_content)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Helper script for genome-based phylogenetic analysis"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # map-ogs
    p = subparsers.add_parser("map-ogs", help="Map eggNOG OGs to Castelli OG files")
    p.add_argument("--eggnog-annotations", required=True, help="eggNOG-mapper .annotations file")
    p.add_argument("--castelli-og-dir", required=True, help="Directory with Castelli single-copy OG files")
    p.add_argument("--output", required=True, help="Output TSV mapping file")

    # extract-and-align
    p = subparsers.add_parser("extract-and-align", help="Extract target taxa, add MAG protein, align")
    p.add_argument("--og-mapping", required=True, help="OG mapping TSV from map-ogs")
    p.add_argument("--castelli-og-dir", required=True, help="Directory with Castelli single-copy OG files")
    p.add_argument("--protein-file", required=True, help="MAG protein FASTA (.faa)")
    p.add_argument("--output-dir", required=True, help="Output directory for aligned OGs")
    p.add_argument("--threads", type=int, default=4, help="Threads for MAFFT")
    p.add_argument("--extra-genomes", default=None,
                   help="Path to GTDB genome_manifest.json for adding extra genomes")

    # trim-alignments
    p = subparsers.add_parser("trim-alignments", help="Trim alignments with BMGE")
    p.add_argument("--input-dir", required=True, help="Directory with aligned OG FASTAs")
    p.add_argument("--output-dir", required=True, help="Output directory for trimmed OGs")

    # concatenate
    p = subparsers.add_parser("concatenate", help="Build supermatrix from trimmed OGs")
    p.add_argument("--input-dir", required=True, help="Directory with trimmed OG FASTAs")
    p.add_argument("--output", required=True, help="Output concatenated FASTA")
    p.add_argument("--partition-file", required=True, help="Output RAxML-style partition file")

    # remove-comp-bias
    p = subparsers.add_parser("remove-comp-bias", help="Remove compositionally biased sites")
    p.add_argument("--alignment", required=True, help="Concatenated alignment FASTA")
    p.add_argument("--threshold", type=float, default=0.30, help="Fraction of sites to remove (default: 0.30)")
    p.add_argument("--output", required=True, help="Output debiased alignment FASTA")

    # summary
    p = subparsers.add_parser("summary", help="Generate analysis report")
    p.add_argument("--output-dir", required=True, help="Output directory")

    args = parser.parse_args()

    if args.command == "map-ogs":
        map_ogs(args.eggnog_annotations, args.castelli_og_dir, args.output)
    elif args.command == "extract-and-align":
        extract_and_align(args.og_mapping, args.castelli_og_dir, args.protein_file,
                          args.output_dir, args.threads, args.extra_genomes)
    elif args.command == "trim-alignments":
        trim_alignments(args.input_dir, args.output_dir)
    elif args.command == "concatenate":
        concatenate(args.input_dir, args.output, args.partition_file)
    elif args.command == "remove-comp-bias":
        remove_comp_bias(args.alignment, args.threshold, args.output)
    elif args.command == "summary":
        generate_summary(args.output_dir)


if __name__ == "__main__":
    main()
