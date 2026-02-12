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
import json
import subprocess
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
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
    "Symbiont_of_S_maritima": "Symbiont_of_Strigamia_maritima",
}

# Mapping from CheckM bin IDs to phylogeny tip names (for our CheckM results)
BIN_ID_TO_TIP = {
    "s7_ctg000008c": "s7_ctg000008c",
    "Terasakiella_pusilla_DSM_6293": "Terasakiella_pusilla_DSM_6293",
    "Thalassospira_profundimaris": "Thalassospira_profundimaris",
}

# Mapping from tip names to NCBI GenBank accessions (with version)
TIP_TO_ACCESSION = {
    "Symbiont_of_Haliotis": "GCA_002632265.1",
    "Marine_MAG_001510075": "GCA_001510075.1",
    "Marine_MAG_001830425": "GCA_001830425.1",
    "Marine_MAG_002327565": "GCA_002327565.1",
    "Marine_MAG_002687515": "GCA_002687515.1",
    "Marine_MAG_009694195": "GCA_009694195.1",
    "Marine_MAG_013204045": "GCA_013204045.1",
    "Marine_MAG_014859895": "GCA_014859895.1",
    "Marine_MAG_018662225": "GCA_018662225.1",
    "Outgroup_009649675": "GCA_009649675.1",
    "Symbiont_of_L_labralis": "GCA_009780035.1",
    "Terasakiella_pusilla_DSM_6293": "GCA_000688235.1",
    "Thalassospira_profundimaris": "GCA_000300275.1",
}

# Host info for tips where NCBI BioSample doesn't have host data
# Values are (host_species, host_family, source)
HARDCODED_HOSTS = {
    "s7_ctg000008c": ("Anchylorhynchus bicarinatus", "Curculionidae", "this_study"),
    "Hepatincola_Av": ("Armadillidium vulgare", "Armadillidiidae", "Dittmer_et_al_2023"),
    "Hepatincola_Pp": ("Porcellionides pruinosus", "Porcellionidae", "Dittmer_et_al_2023"),
    "Hepatincola_Pdp": ("Porcellio dilatatus petiti", "Porcellionidae", "Dittmer_et_al_2023"),
    "Tardigradibacter_bertolanii": ("Richtersius cf. coronifer", "Richtersiidae", "Castelli_et_al_2025"),
    "Symbiont_of_S_maritima": ("Strigamia maritima", "Linotaeniidae", "Castelli_et_al_2025"),
}

# Map NCBI BioSample host names to taxonomic families
HOST_TO_FAMILY = {
    "Labiotermes labralis": "Termitidae",
    "abalone": "Haliotidae",
    "Haliotis discus hannai": "Haliotidae",
    "Cornitermes pugnax": "Termitidae",
    "Jugositermes tuberculatus": "Termitidae",
}

# Special enrichment: BioSample says "abalone", Castelli et al. provide species ID
BIOSAMPLE_HOST_ENRICHMENT = {
    "abalone": ("Haliotis discus hannai", "NCBI_BioSample+Castelli_et_al_2025"),
}


def _fetch_url(url: str, retries: int = 3) -> str:
    """Fetch URL content using system curl (works around conda SSL issues)."""
    for attempt in range(retries):
        result = subprocess.run(
            ["/usr/bin/curl", "-s", "-L", url],
            capture_output=True, text=True, check=True,
        )
        text = result.stdout
        # Retry on empty or HTML error responses (NCBI rate limiting)
        if text.strip() and not text.strip().startswith("<!"):
            return text
        if attempt < retries - 1:
            time.sleep(1 + attempt * 2)
    return result.stdout


def lookup_organism_names(accessions: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    """Query NCBI Entrez for organism names and BioSample accessions.

    Uses the NCBI Assembly database to look up organism names and BioSample IDs.
    Returns (organism_names, biosample_accns) dicts, both mapping accession -> value.
    """
    organism_names = {}
    biosample_accns = {}
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    for accession in accessions:
        try:
            # Search Assembly database for this accession
            term = urllib.parse.quote(f"{accession}[Assembly Accession]")
            search_url = (
                f"{base_url}/esearch.fcgi?"
                f"db=assembly&term={term}&retmode=xml"
            )
            xml_text = _fetch_url(search_url)
            try:
                tree = ET.fromstring(xml_text)
            except ET.ParseError:
                print(f"  WARNING: Bad XML response for {accession} (search), skipping")
                time.sleep(2)
                continue
            id_list = tree.findall(".//Id")
            if not id_list:
                print(f"  WARNING: No assembly found for {accession}")
                continue

            assembly_id = id_list[0].text

            # Fetch assembly summary to get organism name and BioSample accession
            summary_url = (
                f"{base_url}/esummary.fcgi?"
                f"db=assembly&id={assembly_id}&retmode=xml"
            )
            xml_text = _fetch_url(summary_url)
            try:
                tree = ET.fromstring(xml_text)
            except ET.ParseError:
                print(f"  WARNING: Bad XML response for {accession} (summary), skipping")
                time.sleep(2)
                continue

            # Look for Organism element in the DocumentSummary
            organism_elem = tree.find(".//Organism")
            organism = organism_elem.text if organism_elem is not None else None

            if organism:
                organism_names[accession] = organism
                print(f"  {accession} -> {organism}")
            else:
                print(f"  WARNING: No organism name found for {accession}")

            # Look for BioSampleAccn element
            biosample_elem = tree.find(".//BioSampleAccn")
            if biosample_elem is not None and biosample_elem.text:
                biosample_accns[accession] = biosample_elem.text
                print(f"    BioSample: {biosample_elem.text}")

        except Exception as e:
            print(f"  WARNING: Error looking up {accession}: {e}")

        # Be polite to NCBI: rate limit
        time.sleep(0.5)

    return organism_names, biosample_accns


def lookup_biosample_hosts(biosample_accns: dict[str, str]) -> dict[str, str]:
    """Fetch host attribute from NCBI BioSample for each accession.

    Takes dict of assembly_accession -> biosample_accession.
    Returns dict of assembly_accession -> host_name.
    """
    results = {}
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    for assembly_acc, biosample_acc in biosample_accns.items():
        try:
            # Fetch BioSample XML
            fetch_url = (
                f"{base_url}/efetch.fcgi?"
                f"db=biosample&id={biosample_acc}&retmode=xml"
            )
            xml_text = _fetch_url(fetch_url)
            try:
                tree = ET.fromstring(xml_text)
            except ET.ParseError:
                print(f"  {assembly_acc} ({biosample_acc}) -> bad XML response, skipping")
                time.sleep(2)
                continue

            # Look for host attribute
            for attr in tree.findall(".//Attribute"):
                if attr.get("attribute_name") == "host" or attr.get("harmonized_name") == "host":
                    host = attr.text
                    if host and host.lower() not in ("not applicable", "missing", "not collected"):
                        results[assembly_acc] = host
                        print(f"  {assembly_acc} ({biosample_acc}) -> host: {host}")
                    break
            else:
                print(f"  {assembly_acc} ({biosample_acc}) -> no host attribute")
        except Exception as e:
            print(f"  {assembly_acc} ({biosample_acc}) -> error: {e}")

        # Be polite to NCBI: rate limit
        time.sleep(0.5)

    return results


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


def load_gtdb_manifest(manifest_path: Path) -> list[dict]:
    """Load CheckM stats and metadata from a GTDB genome manifest.

    Returns a list of dicts with the same fields as match_castelli_taxa output,
    plus extra metadata for host resolution.
    """
    with open(manifest_path) as fh:
        manifest = json.load(fh)

    results = []
    for genome in manifest["genomes"]:
        if genome["status"] != "new":
            continue

        tip_name = genome.get("tip_name")
        if not tip_name:
            continue

        results.append({
            'tip_name': tip_name,
            'completeness': genome.get("checkm_completeness"),
            'contamination': genome.get("checkm_contamination"),
            'gc_content': genome.get("gc_content"),
            'genome_size': genome.get("genome_size"),
            'n_contigs': genome.get("n_contigs"),
            'data_source': 'this_study_CheckM',
        })
        print(f"  {tip_name} ({genome['accession']})")

    return results, manifest


def _parse_host_from_organism(organism: str) -> str | None:
    """Extract host species from organism name like 'endosymbiont of Pipizella viduata'."""
    if not organism:
        return None
    lower = organism.lower()
    if "endosymbiont of" in lower:
        host = organism.split("endosymbiont of")[-1].strip()
        # Clean up asterisks, quotes, etc.
        host = host.strip("*\"' ")
        if host and host.lower() not in ("", "not applicable", "missing"):
            return host
    return None


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
    parser.add_argument(
        "--gtdb-manifest",
        type=Path,
        default=None,
        help="Path to GTDB genome_manifest.json for enrichment genomes"
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

    # Combine results
    all_results = castelli_results + our_results

    # Load GTDB manifest if provided
    gtdb_manifest = None
    gtdb_tip_accessions = {}  # tip_name -> accession for GTDB genomes
    gtdb_tip_organisms = {}   # tip_name -> ncbi_organism
    gtdb_tip_biosamples = {}  # tip_name -> biosample_accn
    if args.gtdb_manifest and args.gtdb_manifest.exists():
        print("Adding GTDB enrichment genomes:")
        gtdb_results, gtdb_manifest = load_gtdb_manifest(args.gtdb_manifest)
        # Dedup: skip GTDB tips that already appeared in Castelli/our results
        existing_tips = {r['tip_name'] for r in all_results}
        gtdb_results = [r for r in gtdb_results if r['tip_name'] not in existing_tips]
        all_results.extend(gtdb_results)

        # Build lookup dicts for GTDB genomes
        for genome in gtdb_manifest["genomes"]:
            if genome["status"] == "new" and genome.get("tip_name"):
                tip = genome["tip_name"]
                gtdb_tip_accessions[tip] = genome["accession"]
                gtdb_tip_organisms[tip] = genome.get("ncbi_organism", "")
                if genome.get("biosample_accn"):
                    gtdb_tip_biosamples[tip] = genome["biosample_accn"]
        print()

    # Merge TIP_TO_ACCESSION with GTDB genomes for NCBI lookups
    all_tip_to_accession = dict(TIP_TO_ACCESSION)
    all_tip_to_accession.update(gtdb_tip_accessions)

    df = pd.DataFrame(all_results)
    df = df.sort_values('tip_name')

    # Build new_name column using NCBI organism names
    all_accessions = list(all_tip_to_accession.values())
    print("Looking up organism names from NCBI:")
    organism_names, biosample_accns = lookup_organism_names(all_accessions)
    print()

    def make_new_name(tip_name):
        if tip_name == "s7_ctg000008c":
            return "A. bicarinatus symbiont"
        if tip_name == "Symbiont_of_S_maritima":
            return "JBEXBW000000000 endosymbiont of Strigamia maritima"
        accession = all_tip_to_accession.get(tip_name)
        if accession and accession in organism_names:
            return f"{accession} {organism_names[accession]}"
        # Genomes without NCBI accessions (Castelli-only)
        return f"{tip_name} (Castelli et al, 2025)"

    df['new_name'] = df['tip_name'].apply(make_new_name)

    # Build host metadata columns
    # Merge GTDB BioSample accessions into the lookup
    for tip, bs in gtdb_tip_biosamples.items():
        acc = gtdb_tip_accessions.get(tip)
        if acc and acc not in biosample_accns:
            biosample_accns[acc] = bs

    print("Looking up host info from NCBI BioSample:")
    biosample_hosts = lookup_biosample_hosts(biosample_accns)
    print()

    # Host family mappings for GTDB symbiont genomes
    gtdb_host_families = {
        "Pipizella viduata": "Syrphidae",
        "Ecdyonurus torrentis": "Heptageniidae",
    }

    def resolve_host(tip_name):
        """Resolve host_species, host_family, host_source for a tip."""
        # Check hardcoded hosts first (no NCBI data available)
        if tip_name in HARDCODED_HOSTS:
            return HARDCODED_HOSTS[tip_name]

        accession = all_tip_to_accession.get(tip_name)

        # Check NCBI BioSample
        if accession and accession in biosample_hosts:
            host = biosample_hosts[accession]
            # Apply enrichment if available (e.g., "abalone" -> full species)
            if host in BIOSAMPLE_HOST_ENRICHMENT:
                enriched_species, source = BIOSAMPLE_HOST_ENRICHMENT[host]
                family = HOST_TO_FAMILY.get(host, "NA")
                return (enriched_species, family, source)
            family = HOST_TO_FAMILY.get(host, gtdb_host_families.get(host, "NA"))
            return (host, family, "NCBI_BioSample")

        # For GTDB genomes: try parsing host from organism name
        if tip_name in gtdb_tip_organisms:
            organism = gtdb_tip_organisms[tip_name]
            host = _parse_host_from_organism(organism)
            if host:
                family = gtdb_host_families.get(host, HOST_TO_FAMILY.get(host, "NA"))
                return (host, family, "GTDB_organism_name")

        # Not a symbiont (or no host data found)
        return ("NA", "NA", "NA")

    host_data = df['tip_name'].apply(resolve_host)
    df['host_species'] = host_data.apply(lambda x: x[0])
    df['host_family'] = host_data.apply(lambda x: x[1])
    df['host_source'] = host_data.apply(lambda x: x[2])

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
