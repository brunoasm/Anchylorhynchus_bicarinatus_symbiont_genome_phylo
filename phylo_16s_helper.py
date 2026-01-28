#!/usr/bin/env python3
"""
Helper script for 16S rRNA phylogenetic analysis.
Used by 06_16s_phylogenetic_analysis.sh

Subcommands:
  extract-16s      Extract 16S rRNA sequence from GenBank file
  blast-ncbi       Run BLAST against NCBI 16S database
  retrieve-sequences  Fetch sequences and metadata from BLAST results
  parse-metadata   Parse and annotate sequence metadata
  find-outgroup    Find most distant sequence for tree rooting
  summary          Generate analysis report
"""

import argparse
import socket
import sys
import time
from pathlib import Path
from io import StringIO

import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Keywords indicating symbiont status
SYMBIONT_KEYWORDS = [
    'symbiont', 'endosymbiont', 'primary symbiont', 'secondary symbiont',
    'Candidatus', 'intracellular', 'obligate symbiont', 'facultative symbiont',
    'P-symbiont', 'S-symbiont', 'bacteriocyte', 'mycetocyte',
]

# Keywords for organ/tissue parsing
ORGAN_KEYWORDS = {
    'gut': ['gut', 'intestine', 'midgut', 'hindgut', 'foregut', 'digestive tract'],
    'bacteriome': ['bacteriome', 'mycetome', 'bacteriocyte', 'mycetocyte'],
    'hemolymph': ['hemolymph', 'haemolymph', 'blood'],
    'body': ['whole body', 'whole insect', 'body', 'abdomen', 'thorax'],
    'ovary': ['ovary', 'ovaries', 'oocyte', 'egg'],
    'fat_body': ['fat body', 'fat tissue'],
    'salivary_gland': ['salivary gland', 'salivary'],
}


def extract_16s(genbank_path: str, output_path: str, query_name: str = "query_16s") -> None:
    """Extract 16S rRNA sequence from a GenBank file."""
    record = SeqIO.read(genbank_path, "genbank")

    rrna_features = []
    for feature in record.features:
        if feature.type == "rRNA":
            product = feature.qualifiers.get("product", [""])[0]
            if "16S" in product:
                rrna_features.append(feature)

    if not rrna_features:
        print("ERROR: No 16S rRNA feature found in GenBank file", file=sys.stderr)
        sys.exit(1)

    # Take the first 16S rRNA feature
    feature = rrna_features[0]
    seq = feature.extract(record.seq)

    # Create output record
    output_record = SeqRecord(
        seq,
        id=query_name,
        description=f"16S rRNA from {record.id}"
    )

    SeqIO.write(output_record, output_path, "fasta")
    print(f"Extracted 16S rRNA sequence: {len(seq)} bp")
    print(f"Location: {feature.location}")
    print(f"Output: {output_path}")


def blast_ncbi(query_path: str, output_path: str, num_hits: int = 100,
               email: str = "user@example.com", api_key: str = None) -> None:
    """Run BLAST against NCBI 16S ribosomal RNA database using REST API."""
    import urllib.request
    import urllib.parse
    import re

    # Read query sequence
    query_record = SeqIO.read(query_path, "fasta")
    query_seq = str(query_record.seq)

    print(f"Running BLAST against NCBI 16S_ribosomal_RNA database...")
    print(f"Query: {query_record.id} ({len(query_seq)} bp)")
    print(f"Requesting {num_hits} hits...")

    socket.setdefaulttimeout(120)

    blast_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    # Step 1: Submit BLAST job
    submit_params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": "16S_ribosomal_RNA",
        "QUERY": query_seq,
        "MEGABLAST": "on",
        "HITLIST_SIZE": str(num_hits),
        "FORMAT_TYPE": "XML",
    }
    if api_key:
        submit_params["API_KEY"] = api_key

    print("Submitting BLAST job...")
    max_retries = 5
    for attempt in range(max_retries):
        try:
            data = urllib.parse.urlencode(submit_params).encode("utf-8")
            request = urllib.request.Request(blast_url, data=data)
            request.add_header("User-Agent", f"BiopythonClient ({email})")
            response = urllib.request.urlopen(request)
            response_text = response.read().decode("utf-8")
            break
        except Exception as e:
            if attempt < max_retries - 1:
                wait_time = 15 * (2 ** attempt)
                print(f"  Submit attempt {attempt + 1} failed: {e}")
                print(f"  Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                raise

    # Extract RID
    rid_match = re.search(r"RID = (\S+)", response_text)
    rtoe_match = re.search(r"RTOE = (\d+)", response_text)

    if not rid_match:
        print("ERROR: Failed to get BLAST RID from response", file=sys.stderr)
        sys.exit(1)

    rid = rid_match.group(1)
    rtoe = int(rtoe_match.group(1)) if rtoe_match else 60

    print(f"  RID: {rid}")
    print(f"  Estimated time: {rtoe} seconds")
    print(f"  Waiting {rtoe} seconds before first check...")
    time.sleep(rtoe)

    # Step 2: Poll for results
    poll_params = {
        "CMD": "Get",
        "FORMAT_OBJECT": "SearchInfo",
        "RID": rid,
    }
    if api_key:
        poll_params["API_KEY"] = api_key

    while True:
        try:
            poll_url = blast_url + "?" + urllib.parse.urlencode(poll_params)
            request = urllib.request.Request(poll_url)
            request.add_header("User-Agent", f"BiopythonClient ({email})")
            response = urllib.request.urlopen(request)
            status_text = response.read().decode("utf-8")
        except Exception as e:
            print(f"  Poll error: {e}, retrying in 30 seconds...")
            time.sleep(30)
            continue

        if "Status=WAITING" in status_text:
            print("  Still running, checking again in 30 seconds...")
            time.sleep(30)
        elif "Status=FAILED" in status_text:
            print("ERROR: BLAST search failed on NCBI server", file=sys.stderr)
            sys.exit(1)
        elif "Status=UNKNOWN" in status_text:
            print("ERROR: BLAST RID expired or unknown", file=sys.stderr)
            sys.exit(1)
        elif "Status=READY" in status_text:
            if "ThereAreHits=yes" in status_text:
                print("  Results ready! Downloading...")
                break
            else:
                print("WARNING: No hits found")
                break
        else:
            print("  Unknown status, checking again in 30 seconds...")
            time.sleep(30)

    # Step 3: Retrieve results
    get_params = {
        "CMD": "Get",
        "FORMAT_TYPE": "XML",
        "RID": rid,
    }
    if api_key:
        get_params["API_KEY"] = api_key

    get_url = blast_url + "?" + urllib.parse.urlencode(get_params)
    request = urllib.request.Request(get_url)
    request.add_header("User-Agent", f"BiopythonClient ({email})")
    response = urllib.request.urlopen(request)
    result_xml = response.read().decode("utf-8")

    with open(output_path, "w") as f:
        f.write(result_xml)

    print(f"BLAST results saved to: {output_path}")


def retrieve_sequences(blast_xml_path: str, fasta_output: str, metadata_output: str,
                       email: str = "user@example.com", api_key: str = None) -> None:
    """Retrieve sequences and raw metadata from BLAST results."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Rate limit: 10/sec with API key, 3/sec without
    rate_delay = 0.1 if api_key else 0.4

    # Parse BLAST results
    with open(blast_xml_path) as f:
        blast_records = NCBIXML.parse(f)
        blast_record = next(blast_records)

    print(f"Found {len(blast_record.alignments)} BLAST hits")

    sequences = []
    metadata_rows = []

    for i, alignment in enumerate(blast_record.alignments):
        # Extract accession from hit_id or hit_accession
        accession = alignment.accession
        if not accession:
            # Try to extract from hit_id (format: gi|xxx|ref|NR_xxx.x|)
            parts = alignment.hit_id.split('|')
            for j, part in enumerate(parts):
                if part in ('ref', 'gb', 'emb', 'dbj'):
                    accession = parts[j + 1].split('.')[0]
                    break

        if not accession:
            accession = alignment.hit_id

        # Get the best HSP
        hsp = alignment.hsps[0]

        print(f"  [{i+1}/{len(blast_record.alignments)}] Fetching {accession}...")

        # Rate limiting
        time.sleep(rate_delay)

        try:
            # Fetch GenBank record for full metadata
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            gb_text = handle.read()
            handle.close()

            # Parse GenBank record
            gb_record = SeqIO.read(StringIO(gb_text), "genbank")

            # Extract metadata from source features
            organism = gb_record.annotations.get("organism", "")
            source_features = [f for f in gb_record.features if f.type == "source"]

            host = ""
            isolation_source = ""
            strain = ""
            country = ""

            if source_features:
                source = source_features[0]
                host = source.qualifiers.get("host", [""])[0]
                isolation_source = source.qualifiers.get("isolation_source", [""])[0]
                strain = source.qualifiers.get("strain", [""])[0]
                country = source.qualifiers.get("country", [""])[0]

            # Get sequence from BLAST hit (use the subject sequence)
            # We'll fetch the actual sequence from GenBank
            seq_record = SeqRecord(
                gb_record.seq,
                id=accession,
                description=alignment.hit_def
            )
            sequences.append(seq_record)

            metadata_rows.append({
                'accession': accession,
                'organism': organism,
                'hit_def': alignment.hit_def,
                'host': host,
                'isolation_source': isolation_source,
                'strain': strain,
                'country': country,
                'identity': hsp.identities / hsp.align_length * 100,
                'alignment_length': hsp.align_length,
                'evalue': hsp.expect,
                'bit_score': hsp.bits,
            })

        except Exception as e:
            print(f"    Warning: Failed to fetch {accession}: {e}")
            continue

    # Write sequences
    SeqIO.write(sequences, fasta_output, "fasta")
    print(f"Wrote {len(sequences)} sequences to: {fasta_output}")

    # Write metadata
    metadata_df = pd.DataFrame(metadata_rows)
    metadata_df.to_csv(metadata_output, sep='\t', index=False)
    print(f"Wrote metadata to: {metadata_output}")


def _str_or_empty(val) -> str:
    """Convert a value to string, treating NaN/None as empty string."""
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return ""
    return str(val).strip()


def is_symbiont(organism: str, host, isolation_source) -> bool:
    """Determine if an organism is likely a symbiont."""
    organism = _str_or_empty(organism)
    host = _str_or_empty(host)
    isolation_source = _str_or_empty(isolation_source)
    text_to_check = (organism + " " + isolation_source).lower()

    # Check for explicit host field
    if host:
        return True

    # Check for symbiont keywords
    for keyword in SYMBIONT_KEYWORDS:
        if keyword.lower() in text_to_check:
            return True

    return False


def parse_host_from_organism(organism: str) -> str:
    """Try to extract host name from organism name patterns."""
    organism_lower = organism.lower()

    # Pattern: "symbiont of Genus species"
    patterns = [
        'symbiont of ', 'endosymbiont of ', 'primary symbiont of ',
        'secondary symbiont of ', 'intracellular symbiont of ',
    ]

    for pattern in patterns:
        if pattern in organism_lower:
            idx = organism_lower.find(pattern) + len(pattern)
            host_part = organism[idx:].strip()
            # Clean up - take first two words (genus species)
            words = host_part.split()
            if len(words) >= 2:
                return f"{words[0]}_{words[1]}"
            elif words:
                return words[0]

    return ""


def parse_organ(isolation_source: str) -> str:
    """Parse organ/tissue from isolation source."""
    if not isolation_source:
        return "NA"

    source_lower = isolation_source.lower()

    for organ, keywords in ORGAN_KEYWORDS.items():
        for keyword in keywords:
            if keyword in source_lower:
                return organ

    return "NA"


def get_host_taxonomy_entrez(host_name: str, email: str, api_key: str = None) -> dict:
    """
    Get taxonomic information for a host using NCBI Taxonomy via Entrez.
    Returns dict with taxon_id, genus, family, order, phylum (names and IDs).

    Falls back through progressively shorter names:
    full name -> genus + species -> genus only.
    """
    result = {
        "taxon_id": "NA",
        "genus": "NA", "genus_id": "NA",
        "family": "NA", "family_id": "NA",
        "order": "NA", "order_id": "NA",
        "phylum": "NA", "phylum_id": "NA",
    }

    if not host_name or host_name == "NA":
        return result

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Rate limit: 10/sec with API key, 3/sec without
    rate_delay = 0.1 if api_key else 0.4

    # Clean host name for search
    search_term = host_name.replace('_', ' ')
    words = search_term.split()

    # Build list of search terms: full name -> genus species -> genus
    # Use dict.fromkeys to deduplicate while preserving order
    candidates = [search_term]
    if len(words) >= 2:
        candidates.append(f"{words[0]} {words[1]}")
    if words:
        candidates.append(words[0])
    search_terms = list(dict.fromkeys(candidates))

    tax_id = None
    for term in search_terms:
        try:
            time.sleep(rate_delay)
            handle = Entrez.esearch(db="taxonomy", term=term, retmax=1)
            record = Entrez.read(handle)
            handle.close()

            if record["IdList"]:
                tax_id = record["IdList"][0]
                break
        except Exception as e:
            print(f"    Warning: Search failed for '{term}': {e}", file=sys.stderr)
            continue

    if not tax_id:
        return result

    result["taxon_id"] = tax_id

    try:
        # Fetch the taxonomy record
        time.sleep(rate_delay)
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        if not records:
            return result

        # Look through lineage for genus, family, order, phylum
        lineage = records[0].get("LineageEx", [])
        for taxon in lineage:
            rank = taxon.get("Rank")
            name = taxon.get("ScientificName", "NA")
            taxon_id = taxon.get("TaxId", "NA")
            if rank == "genus":
                result["genus"] = name
                result["genus_id"] = taxon_id
            elif rank == "family":
                result["family"] = name
                result["family_id"] = taxon_id
            elif rank == "order":
                result["order"] = name
                result["order_id"] = taxon_id
            elif rank == "phylum":
                result["phylum"] = name
                result["phylum_id"] = taxon_id

        return result

    except Exception as e:
        print(f"    Warning: Failed to fetch taxonomy for {host_name}: {e}", file=sys.stderr)
        return result


def clean_taxon_name(organism: str) -> str:
    """Clean organism name for use as taxon name."""
    # Remove common prefixes
    prefixes_to_remove = [
        'Candidatus ', 'uncultured ', 'unidentified ',
    ]

    name = organism
    for prefix in prefixes_to_remove:
        if name.startswith(prefix):
            name = name[len(prefix):]

    # Replace spaces with underscores
    name = name.replace(' ', '_')

    # Remove special characters
    name = ''.join(c for c in name if c.isalnum() or c == '_')

    return name


def parse_metadata(raw_metadata_path: str, output_path: str, email: str, api_key: str = None) -> None:
    """Parse and annotate sequence metadata."""
    df = pd.read_csv(raw_metadata_path, sep='\t')

    annotated_rows = []

    # Cache for host taxonomy lookups to avoid redundant API calls
    host_taxonomy_cache = {}

    for idx, row in df.iterrows():
        organism = _str_or_empty(row['organism'])
        host = _str_or_empty(row['host'])
        isolation_source = _str_or_empty(row['isolation_source'])

        # Determine symbiont status
        symbiont_status = "symbiont" if is_symbiont(organism, host, isolation_source) else "non-symbiont"

        # Get host name
        if host:
            host_name = host.replace(' ', '_')
            annotation_source = "host_field"
        else:
            host_name = parse_host_from_organism(organism)
            annotation_source = "organism_name" if host_name else "NA"

        if not host_name:
            host_name = "NA"

        # Get organ
        isolation_organ = parse_organ(isolation_source)

        # Get host taxonomy via Entrez with caching
        if host_name != "NA":
            if host_name not in host_taxonomy_cache:
                print(f"  [{idx+1}/{len(df)}] Looking up taxonomy for: {host_name}")
                host_taxonomy_cache[host_name] = get_host_taxonomy_entrez(host_name, email, api_key)
            host_taxonomy = host_taxonomy_cache[host_name]
        else:
            host_taxonomy = {
                "taxon_id": "NA",
                "genus": "NA", "genus_id": "NA",
                "family": "NA", "family_id": "NA",
                "order": "NA", "order_id": "NA",
                "phylum": "NA", "phylum_id": "NA",
            }

        # Clean taxon name
        taxon_name = clean_taxon_name(organism)

        annotated_rows.append({
            'accession': row['accession'],
            'taxon_name': taxon_name,
            'symbiont_status': symbiont_status,
            'host_name': host_name,
            'host_taxon_id': host_taxonomy["taxon_id"],
            'host_genus': host_taxonomy["genus"],
            'host_genus_id': host_taxonomy["genus_id"],
            'host_family': host_taxonomy["family"],
            'host_family_id': host_taxonomy["family_id"],
            'host_order': host_taxonomy["order"],
            'host_order_id': host_taxonomy["order_id"],
            'host_phylum': host_taxonomy["phylum"],
            'host_phylum_id': host_taxonomy["phylum_id"],
            'isolation_organ': isolation_organ,
            'annotation_source': annotation_source,
        })

    annotated_df = pd.DataFrame(annotated_rows)
    annotated_df.to_csv(output_path, sep='\t', index=False)

    # Print summary
    symbiont_count = (annotated_df['symbiont_status'] == 'symbiont').sum()
    host_known_count = (annotated_df['host_name'] != 'NA').sum()
    host_order_known = (annotated_df['host_order'] != 'NA').sum()
    host_phylum_known = (annotated_df['host_phylum'] != 'NA').sum()

    print(f"Metadata annotation complete:")
    print(f"  Total sequences: {len(annotated_df)}")
    print(f"  Symbionts: {symbiont_count}")
    print(f"  Non-symbionts: {len(annotated_df) - symbiont_count}")
    print(f"  With known host: {host_known_count}")
    print(f"  With host order resolved: {host_order_known}")
    print(f"  With host phylum resolved: {host_phylum_known}")
    print(f"  Unique hosts looked up: {len(host_taxonomy_cache)}")
    print(f"Output: {output_path}")


def find_outgroup(alignment_path: str) -> str:
    """Find the most distant sequence to use as outgroup."""
    from Bio import AlignIO

    alignment = AlignIO.read(alignment_path, "fasta")
    n_seqs = len(alignment)
    seq_length = alignment.get_alignment_length()

    print(f"Analyzing alignment: {n_seqs} sequences, {seq_length} positions")

    # Calculate pairwise distances (simple identity-based)
    distances = {}
    for i, rec1 in enumerate(alignment):
        total_dist = 0
        for j, rec2 in enumerate(alignment):
            if i != j:
                # Calculate identity
                matches = sum(a == b for a, b in zip(str(rec1.seq), str(rec2.seq)))
                identity = matches / seq_length
                total_dist += (1 - identity)
        distances[rec1.id] = total_dist / (n_seqs - 1)

    # Find sequence with highest average distance
    outgroup = max(distances, key=distances.get)
    avg_dist = distances[outgroup]

    print(f"Most distant sequence: {outgroup}")
    print(f"Average distance: {avg_dist:.4f}")

    # Print just the outgroup ID for shell script capture
    return outgroup


def generate_summary(output_dir: str, query_name: str = "query_16s") -> None:
    """Generate markdown summary report."""
    output_dir = Path(output_dir)

    report_lines = [
        "# 16S rRNA Phylogenetic Analysis Report",
        "",
        "## Overview",
        "",
    ]

    # Check for blast hits
    blast_hits_path = output_dir / "blast_hits.fasta"
    if blast_hits_path.exists():
        sequences = list(SeqIO.parse(blast_hits_path, "fasta"))
        report_lines.append(f"- **BLAST hits retrieved**: {len(sequences)}")

    # Check for metadata
    metadata_path = output_dir / "sequence_metadata.tsv"
    if metadata_path.exists():
        df = pd.read_csv(metadata_path, sep='\t')
        symbiont_count = (df['symbiont_status'] == 'symbiont').sum()
        host_known = (df['host_name'] != 'NA').sum()

        report_lines.extend([
            f"- **Symbionts identified**: {symbiont_count}",
            f"- **Non-symbionts**: {len(df) - symbiont_count}",
            f"- **Sequences with known host**: {host_known}",
            "",
        ])

        # Host order breakdown
        order_counts = df['host_order'].value_counts()
        if len(order_counts) > 1 or (len(order_counts) == 1 and order_counts.index[0] != 'NA'):
            report_lines.extend([
                "## Host Orders",
                "",
                "| Order | Count |",
                "|-------|-------|",
            ])
            for order, count in order_counts.items():
                report_lines.append(f"| {order} | {count} |")
            report_lines.append("")

    # Check for alignment
    alignment_path = output_dir / "alignment.fasta"
    if alignment_path.exists():
        from Bio import AlignIO
        alignment = AlignIO.read(str(alignment_path), "fasta")
        report_lines.extend([
            "## Alignment",
            "",
            f"- **Total sequences**: {len(alignment)} (query + hits)",
            f"- **Alignment length**: {alignment.get_alignment_length()} bp",
            "",
        ])

    # Check for tree
    tree_path = output_dir / "16s_tree.treefile"
    iqtree_report_path = output_dir / "16s_tree.iqtree"
    if tree_path.exists():
        report_lines.extend([
            "## Phylogenetic Tree",
            "",
            f"- **Tree file**: `16s_tree.treefile`",
        ])

        # Try to extract model from IQ-TREE report
        if iqtree_report_path.exists():
            with open(iqtree_report_path) as f:
                for line in f:
                    if "Best-fit model:" in line:
                        model = line.split(":")[-1].strip()
                        report_lines.append(f"- **Best-fit model**: {model}")
                        break

        report_lines.extend([
            "- **Bootstrap replicates**: 1000 (ultrafast)",
            "",
        ])

    # Output files section
    report_lines.extend([
        "## Output Files",
        "",
        "| File | Description |",
        "|------|-------------|",
        "| `query_16s.fasta` | Extracted 16S from query genome |",
        "| `blast_results.xml` | Raw NCBI BLAST output |",
        "| `blast_hits.fasta` | Retrieved sequences (original IDs) |",
        "| `blast_hits_raw_metadata.tsv` | Raw metadata from NCBI |",
        "| `sequence_metadata.tsv` | Annotated metadata table |",
        "| `all_sequences.fasta` | Query + hits combined |",
        "| `alignment.fasta` | MAFFT alignment |",
        "| `16s_tree.treefile` | Rooted Newick tree |",
        "| `16s_tree.iqtree` | IQ-TREE report |",
        "",
    ])

    report_content = '\n'.join(report_lines)
    report_path = output_dir / "REPORT.md"

    with open(report_path, 'w') as f:
        f.write(report_content)

    print(f"Report generated: {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Helper script for 16S rRNA phylogenetic analysis'
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # extract-16s command
    extract_parser = subparsers.add_parser('extract-16s', help='Extract 16S rRNA from GenBank')
    extract_parser.add_argument('--genbank', required=True, help='Path to GenBank file')
    extract_parser.add_argument('--output', required=True, help='Output FASTA path')
    extract_parser.add_argument('--name', default='query_16s', help='Sequence name')

    # blast-ncbi command
    blast_parser = subparsers.add_parser('blast-ncbi', help='Run BLAST against NCBI 16S database')
    blast_parser.add_argument('--query', required=True, help='Query FASTA file')
    blast_parser.add_argument('--output', required=True, help='Output XML path')
    blast_parser.add_argument('--num-hits', type=int, default=100, help='Number of hits')
    blast_parser.add_argument('--email', required=True, help='Email for NCBI Entrez')
    blast_parser.add_argument('--api-key', help='NCBI API key (optional, increases rate limit)')

    # retrieve-sequences command
    retrieve_parser = subparsers.add_parser('retrieve-sequences', help='Retrieve sequences from BLAST results')
    retrieve_parser.add_argument('--blast-xml', required=True, help='BLAST XML file')
    retrieve_parser.add_argument('--fasta-output', required=True, help='Output FASTA')
    retrieve_parser.add_argument('--metadata-output', required=True, help='Output metadata TSV')
    retrieve_parser.add_argument('--email', required=True, help='Email for NCBI Entrez')
    retrieve_parser.add_argument('--api-key', help='NCBI API key (optional, increases rate limit)')

    # parse-metadata command
    metadata_parser = subparsers.add_parser('parse-metadata', help='Parse and annotate metadata')
    metadata_parser.add_argument('--raw-metadata', required=True, help='Raw metadata TSV')
    metadata_parser.add_argument('--output', required=True, help='Annotated metadata output')
    metadata_parser.add_argument('--email', required=True, help='Email for NCBI Entrez')
    metadata_parser.add_argument('--api-key', help='NCBI API key (optional, increases rate limit)')

    # find-outgroup command
    outgroup_parser = subparsers.add_parser('find-outgroup', help='Find most distant sequence')
    outgroup_parser.add_argument('--alignment', required=True, help='Alignment FASTA file')

    # summary command
    summary_parser = subparsers.add_parser('summary', help='Generate summary report')
    summary_parser.add_argument('--output-dir', required=True, help='Output directory')
    summary_parser.add_argument('--query-name', default='query_16s', help='Query sequence name')

    args = parser.parse_args()

    if args.command == 'extract-16s':
        extract_16s(args.genbank, args.output, args.name)

    elif args.command == 'blast-ncbi':
        blast_ncbi(args.query, args.output, args.num_hits, args.email, args.api_key)

    elif args.command == 'retrieve-sequences':
        retrieve_sequences(args.blast_xml, args.fasta_output, args.metadata_output, args.email, args.api_key)

    elif args.command == 'parse-metadata':
        parse_metadata(args.raw_metadata, args.output, args.email, args.api_key)

    elif args.command == 'find-outgroup':
        outgroup = find_outgroup(args.alignment)
        # Print just the outgroup ID for shell script capture
        print(outgroup)

    elif args.command == 'summary':
        generate_summary(args.output_dir, args.query_name)


if __name__ == '__main__':
    main()
