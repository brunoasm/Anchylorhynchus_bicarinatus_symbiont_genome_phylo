#!/usr/bin/env python3
"""
Helper script for 16S rRNA phylogenetic analysis.
Used by 06_16s_phylogenetic_analysis.sh

Subcommands:
  extract-16s             Extract 16S rRNA sequence from GenBank file
  blast-ncbi              Run BLAST against NCBI 16S database
  retrieve-sequences      Fetch sequences and metadata from BLAST results
  parse-silva             Parse SILVA nearest-neighbour FASTA headers
  retrieve-silva-sequences  Fetch sequences for SILVA hits from GenBank
  parse-gtdb-wrau01       Extract WRAU01 genome accessions from GTDB MSA
  retrieve-gtdb-16s       Fetch 16S sequences from GTDB WRAU01 genomes
  merge-sources           Merge FASTA and metadata from all sources
  parse-metadata          Parse and annotate sequence metadata
  find-outgroup           Find most distant sequence for tree rooting
  summary                 Generate analysis report
"""

import argparse
import gzip
import os
import re
import socket
import sys
import time
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML, NCBIWWW
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

# Microbial taxonomy ranks and prefixes
MICROBIAL_RANKS = [
    ('domain', 'd__'),
    ('phylum', 'p__'),
    ('class', 'c__'),
    ('order', 'o__'),
    ('family', 'f__'),
    ('genus', 'g__'),
    ('species', 's__'),
]


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
    import urllib.parse
    import urllib.request

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
                'source': 'blast',
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


# ---------------------------------------------------------------------------
# SILVA subcommands
# ---------------------------------------------------------------------------

def parse_silva(silva_fasta_path: str, output_tsv: str) -> None:
    """Parse SILVA nearest-neighbour FASTA to extract accessions and similarity.

    SILVA headers contain a ``[nearest_slv=...]`` block with entries like::

        AC001339.1055.2587~96.8;JN392909.1.1490~96.2;...

    Each entry is ``accession.version.start.end~similarity``.  When the same
    NCBI accession appears more than once (multiple 16S copies on the same
    genome), we keep only the entry with the longest region (end - start).
    """
    entries = []

    with open(silva_fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue

            # Extract the nearest_slv block
            match = re.search(r'\[nearest_slv=([^\]]+)\]', line)
            if not match:
                continue

            raw = match.group(1)
            # SILVA uses spaces between entries (not semicolons)
            for token in raw.split():
                token = token.strip()
                if not token:
                    continue

                # Format: accession.version.start.end~similarity
                if '~' not in token:
                    continue

                id_part, sim_str = token.rsplit('~', 1)
                try:
                    similarity = float(sim_str)
                except ValueError:
                    continue

                # Split id_part on '.' -- expect accession.version.start.end
                parts = id_part.split('.')
                if len(parts) < 4:
                    # Fallback: store as-is
                    entries.append({
                        'silva_id': id_part,
                        'accession': parts[0] if parts else id_part,
                        'version': parts[1] if len(parts) > 1 else '',
                        'start': 0,
                        'end': 0,
                        'similarity': similarity,
                    })
                    continue

                accession = parts[0]
                version = parts[1]
                try:
                    start = int(parts[2])
                    end = int(parts[3])
                except ValueError:
                    start = 0
                    end = 0

                entries.append({
                    'silva_id': id_part,
                    'accession': accession,
                    'version': version,
                    'start': start,
                    'end': end,
                    'similarity': similarity,
                })

    if not entries:
        print("WARNING: No nearest_slv entries found in SILVA FASTA", file=sys.stderr)
        df = pd.DataFrame(columns=['silva_id', 'accession', 'version', 'start', 'end', 'similarity'])
        df.to_csv(output_tsv, sep='\t', index=False)
        return

    df = pd.DataFrame(entries)
    n_total = len(df)

    # Deduplicate by accession: keep longest region per accession
    df['region_length'] = (df['end'] - df['start']).abs()
    df = df.sort_values('region_length', ascending=False).drop_duplicates(subset='accession', keep='first')
    df = df.drop(columns='region_length').sort_values('similarity', ascending=False).reset_index(drop=True)

    n_dedup = len(df)
    print(f"Parsed {n_total} SILVA entries, {n_dedup} unique accessions (deduplicated by longest region)")
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Output: {output_tsv}")


def retrieve_silva_sequences(silva_tsv_path: str, fasta_output: str, metadata_output: str,
                             email: str = "user@example.com", api_key: str = None) -> None:
    """Fetch 16S subsequences from GenBank for SILVA hits."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    rate_delay = 0.1 if api_key else 0.4

    df = pd.read_csv(silva_tsv_path, sep='\t')
    print(f"Retrieving sequences for {len(df)} SILVA accessions...")

    sequences = []
    metadata_rows = []

    for idx, row in df.iterrows():
        acc = row['accession']
        ver = row.get('version', '')
        start = int(row.get('start', 0))
        end = int(row.get('end', 0))
        accession_full = f"{acc}.{ver}" if ver else acc

        print(f"  [{idx+1}/{len(df)}] Fetching {accession_full} ({start}-{end})...")
        time.sleep(rate_delay)

        try:
            # Fetch with subsequence coordinates if available
            fetch_kwargs = dict(db="nucleotide", id=accession_full, rettype="gb", retmode="text")
            if start > 0 and end > 0:
                fetch_kwargs['seq_start'] = str(start)
                fetch_kwargs['seq_stop'] = str(end)

            handle = Entrez.efetch(**fetch_kwargs)
            gb_text = handle.read()
            handle.close()

            gb_record = SeqIO.read(StringIO(gb_text), "genbank")

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

            seq_record = SeqRecord(
                gb_record.seq,
                id=acc,
                description=f"SILVA hit {accession_full} region {start}-{end}"
            )
            sequences.append(seq_record)

            metadata_rows.append({
                'accession': acc,
                'organism': organism,
                'hit_def': gb_record.description,
                'host': host,
                'isolation_source': isolation_source,
                'strain': strain,
                'country': country,
                'identity': row.get('similarity', ''),
                'alignment_length': '',
                'evalue': '',
                'bit_score': '',
                'source': 'silva',
            })

        except Exception as e:
            print(f"    Warning: Failed to fetch {accession_full}: {e}")
            continue

    SeqIO.write(sequences, fasta_output, "fasta")
    print(f"Wrote {len(sequences)} SILVA sequences to: {fasta_output}")

    metadata_df = pd.DataFrame(metadata_rows)
    metadata_df.to_csv(metadata_output, sep='\t', index=False)
    print(f"Wrote SILVA metadata to: {metadata_output}")


# ---------------------------------------------------------------------------
# GTDB WRAU01 subcommands
# ---------------------------------------------------------------------------

def parse_gtdb_wrau01(msa_path: str, output_path: str) -> None:
    """Extract WRAU01-classified genome accessions from a GTDB MSA FASTA file.

    The GTDB MSA file (``gtdbtk.bac120.msa.fasta.gz``) contains headers like::

        >GB_GCA_009780035.1  d__Bacteria;...;o__WRAU01;...

    We search for headers containing 'WRAU01' and extract the NCBI genome
    accession (stripping the ``GB_`` or ``RS_`` prefix).
    """
    open_fn = gzip.open if str(msa_path).endswith('.gz') else open

    accessions = []
    with open_fn(msa_path, 'rt') as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            if 'WRAU01' not in line:
                continue

            # Header format: >GB_GCA_009780035.1 or >RS_GCF_...
            header_id = line[1:].split()[0]
            # Strip GB_ or RS_ prefix
            if header_id.startswith(('GB_', 'RS_')):
                ncbi_acc = header_id[3:]
            else:
                ncbi_acc = header_id

            accessions.append(ncbi_acc)

    accessions = sorted(set(accessions))

    with open(output_path, 'w') as fh:
        for acc in accessions:
            fh.write(acc + '\n')

    print(f"Found {len(accessions)} WRAU01 genome accessions in GTDB MSA")
    for acc in accessions:
        print(f"  {acc}")
    print(f"Output: {output_path}")


def _get_assembly_metadata(assembly_acc: str, rate_delay: float) -> dict:
    """Fetch metadata from NCBI Assembly and BioSample for a genome assembly.

    Returns a dict with:
        - organism: organism name from assembly
        - host: host from biosample (if available)
        - isolation_source: isolation source from biosample (if available)
        - strain: strain from biosample (if available)
    """
    result = {
        'organism': '',
        'host': '',
        'isolation_source': '',
        'strain': '',
    }

    try:
        # Search for assembly
        time.sleep(rate_delay)
        handle = Entrez.esearch(db="assembly", term=assembly_acc, retmax=1)
        search_result = Entrez.read(handle)
        handle.close()

        assembly_ids = search_result.get("IdList", [])
        if not assembly_ids:
            return result

        assembly_uid = assembly_ids[0]

        # Fetch assembly summary
        time.sleep(rate_delay)
        handle = Entrez.esummary(db="assembly", id=assembly_uid)
        summary = Entrez.read(handle, validate=False)
        handle.close()

        doc = summary['DocumentSummarySet']['DocumentSummary'][0]
        result['organism'] = doc.get('Organism', '')
        biosample_acc = doc.get('BioSampleAccn', '')

        # Fetch BioSample if available
        if biosample_acc:
            try:
                time.sleep(rate_delay)
                handle = Entrez.efetch(db="biosample", id=biosample_acc, retmode="xml")
                biosample_xml = handle.read()
                handle.close()

                # Parse BioSample XML for attributes
                import xml.etree.ElementTree as ET
                root = ET.fromstring(biosample_xml)

                for attr in root.findall(".//Attribute"):
                    attr_name = attr.get('attribute_name', '').lower()
                    attr_value = attr.text or ''

                    if attr_name == 'host':
                        result['host'] = attr_value
                    elif attr_name in ('isolation_source', 'isolation source'):
                        result['isolation_source'] = attr_value
                    elif attr_name == 'strain':
                        result['strain'] = attr_value
            except Exception as e:
                # BioSample fetch failed, but we still have organism from assembly
                pass

    except Exception as e:
        # Return empty result if assembly fetch fails
        pass

    return result


def _get_wgs_contig_ids(assembly_acc: str, rate_delay: float) -> list:
    """Resolve an assembly accession to its WGS contig IDs via NCBI.

    Strategy:
    1. Search the assembly DB for the accession.
    2. Use elink to get the WGS master nucleotide record.
    3. Fetch the master record to extract the WGS contig range.
    4. Return the list of contig accessions.
    """
    # Step 1: find assembly UID
    time.sleep(rate_delay)
    handle = Entrez.esearch(db="assembly", term=assembly_acc, retmax=1)
    result = Entrez.read(handle)
    handle.close()

    assembly_ids = result.get("IdList", [])
    if not assembly_ids:
        return []

    assembly_uid = assembly_ids[0]

    # Step 2: elink to nucleotide (gets WGS master)
    time.sleep(rate_delay)
    handle = Entrez.elink(dbfrom="assembly", db="nuccore", id=assembly_uid)
    links = Entrez.read(handle)
    handle.close()

    nuc_ids = []
    for linkset in links:
        for db_link in linkset.get("LinkSetDb", []):
            for link in db_link["Link"]:
                nuc_ids.append(link["Id"])

    if not nuc_ids:
        return []

    # Step 3: fetch the WGS master to get contig range
    time.sleep(rate_delay)
    handle = Entrez.efetch(db="nucleotide", id=nuc_ids[0], rettype="gb", retmode="text")
    text = handle.read()
    handle.close()

    # Parse WGS line, e.g. "WGS         WRAU01000001-WRAU01000017"
    wgs_match = re.search(r'WGS\s+(\w+)-(\w+)', text)
    if not wgs_match:
        # Not a WGS record -- the nuc_ids themselves are the contigs
        return nuc_ids

    first = wgs_match.group(1)
    last = wgs_match.group(2)

    # Extract prefix and numeric range
    prefix_match = re.match(r'([A-Z]+\d+)(\d{6,})', first)
    if not prefix_match:
        return nuc_ids

    prefix = prefix_match.group(1)
    start_num = int(prefix_match.group(2))
    width = len(prefix_match.group(2))

    last_match = re.match(r'([A-Z]+\d+)(\d{6,})', last)
    end_num = int(last_match.group(2)) if last_match else start_num

    contig_ids = [f"{prefix}{i:0{width}d}" for i in range(start_num, end_num + 1)]
    return contig_ids


def _blast_16s_from_contigs(contig_ids: list, query_fasta: str,
                            assembly_acc: str, rate_delay: float,
                            source: str = "gtdb",
                            assembly_metadata: dict = None) -> tuple:
    """BLAST query 16S against genome contigs to find unannotated 16S.

    Downloads contigs as FASTA, builds a local BLAST database, and runs
    blastn.  Returns ``(SeqRecord, organism, metadata_dict)`` on success,
    or ``None`` if no hit passes filters (>500 bp, >70% identity).

    If assembly_metadata is provided, it will be used for organism/host info.
    """
    if assembly_metadata is None:
        assembly_metadata = {
            'organism': '',
            'host': '',
            'isolation_source': '',
            'strain': '',
        }
    import subprocess
    import tempfile

    print(f"    Fallback: BLASTing query 16S against {assembly_acc} contigs...")

    tmpdir = tempfile.mkdtemp(prefix="gtdb_blast_")
    contigs_fasta = Path(tmpdir) / "contigs.fasta"
    blast_db = str(contigs_fasta)

    try:
        # Fetch contigs as FASTA (in batches to avoid timeouts)
        with open(contigs_fasta, "w") as fh:
            batch_size = 50
            for i in range(0, len(contig_ids), batch_size):
                batch = contig_ids[i:i + batch_size]
                time.sleep(rate_delay)
                handle = Entrez.efetch(db="nucleotide", id=",".join(batch),
                                       rettype="fasta", retmode="text")
                fh.write(handle.read())
                handle.close()

        # Build BLAST database
        subprocess.run(
            ["makeblastdb", "-in", str(contigs_fasta), "-dbtype", "nucl"],
            capture_output=True, check=True,
        )

        # Run BLAST
        result = subprocess.run(
            ["blastn", "-query", query_fasta, "-db", blast_db,
             "-outfmt", "6 sseqid sstart send slen length pident evalue bitscore",
             "-evalue", "1e-10", "-max_target_seqs", "1"],
            capture_output=True, text=True, check=True,
        )

        if not result.stdout.strip():
            print(f"    No BLAST hit found in {assembly_acc} contigs")
            return None

        # Parse best hit
        fields = result.stdout.strip().split("\n")[0].split("\t")
        sseqid = fields[0]
        sstart = int(fields[1])
        send = int(fields[2])
        aln_length = int(fields[4])
        pident = float(fields[5])
        evalue = fields[6]
        bitscore = fields[7]

        if aln_length < 500 or pident < 70.0:
            print(f"    BLAST hit too short or low identity: {aln_length} bp, {pident:.1f}%")
            return None

        # Extract the hit region from the contig FASTA
        for record in SeqIO.parse(str(contigs_fasta), "fasta"):
            if record.id == sseqid:
                start = min(sstart, send) - 1  # 0-based
                end = max(sstart, send)
                seq_16s = record.seq[start:end]
                if sstart > send:
                    seq_16s = seq_16s.reverse_complement()

                # Get organism from assembly metadata first, fallback to contig FASTA header
                organism = assembly_metadata.get('organism', '')
                if not organism:
                    for rec in SeqIO.parse(str(contigs_fasta), "fasta"):
                        desc = rec.description
                        # NCBI FASTA headers: "accession description [organism]"
                        org_match = re.search(r'\[([^\]]+)\]', desc)
                        if org_match:
                            organism = org_match.group(1)
                        break

                print(f"    BLAST fallback found 16S: {len(seq_16s)} bp, "
                      f"{pident:.1f}% identity from {sseqid}")

                seq_record = SeqRecord(
                    seq_16s,
                    id=assembly_acc,
                    description=f"16S rRNA from {organism} ({assembly_acc}) [BLAST fallback]"
                )

                meta = {
                    'accession': assembly_acc,
                    'organism': organism,
                    'hit_def': f"16S rRNA from {organism} [BLAST fallback]",
                    'host': assembly_metadata.get('host', ''),
                    'isolation_source': assembly_metadata.get('isolation_source', ''),
                    'strain': assembly_metadata.get('strain', ''),
                    'country': '',
                    'identity': pident,
                    'alignment_length': aln_length,
                    'evalue': evalue,
                    'bit_score': bitscore,
                    'source': source,
                }
                return (seq_record, meta)

        return None

    except Exception as e:
        print(f"    BLAST fallback failed for {assembly_acc}: {e}")
        return None

    finally:
        # Clean up temp files
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


def retrieve_gtdb_16s(accessions_path: str, fasta_output: str, metadata_output: str,
                      email: str = "user@example.com", api_key: str = None,
                      query_fasta: str = None, source: str = "gtdb") -> None:
    """Fetch 16S rRNA sequences from NCBI for genome assembly accessions.

    For each genome assembly accession:
    1. Resolve to WGS contig IDs via assembly DB + elink
    2. Fetch GenBank records in batches and parse for 16S rRNA features
    3. Extract the 16S sequences (keep longest per assembly)
    4. If no annotation found and *query_fasta* is provided, fall back to
       BLASTing the query 16S against the genome's contigs

    The *source* label (default ``'gtdb'``) is written into each metadata row.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    rate_delay = 0.1 if api_key else 0.4

    with open(accessions_path) as fh:
        accessions = [line.strip() for line in fh if line.strip()]

    # Load existing results to enable incremental retrieval
    existing_accessions = set()
    sequences = []
    metadata_rows = []

    if os.path.isfile(fasta_output) and os.path.getsize(fasta_output) > 0:
        for rec in SeqIO.parse(fasta_output, "fasta"):
            existing_accessions.add(rec.id)
            sequences.append(rec)
    if os.path.isfile(metadata_output) and os.path.getsize(metadata_output) > 0:
        existing_meta = pd.read_csv(metadata_output, sep='\t')
        metadata_rows = existing_meta.to_dict('records')

    remaining = [a for a in accessions if a not in existing_accessions]

    if not remaining:
        print(f"All {len(accessions)} {source.upper()} genomes already retrieved ({len(sequences)} sequences).")
        return

    print(f"Retrieving 16S sequences for {len(remaining)} {source.upper()} genomes "
          f"({len(existing_accessions)} already on disk)...")
    if query_fasta:
        print(f"  BLAST fallback enabled (query: {query_fasta})")

    no_16s_genomes = []

    # Transient error keywords for retry logic
    _transient = ('TXCLIENT', 'IncompleteRead', 'ConnectionReset',
                  'timeout', 'EOF', 'URLError', 'RemoteDisconnected')

    for idx, assembly_acc in enumerate(remaining):
        print(f"  [{idx+1}/{len(remaining)}] Searching {assembly_acc}...")

        max_retries = 5
        for attempt in range(max_retries):
            try:
                # Fetch assembly-level metadata (organism, host, etc.)
                assembly_metadata = _get_assembly_metadata(assembly_acc, rate_delay)

                contig_ids = _get_wgs_contig_ids(assembly_acc, rate_delay)

                if not contig_ids:
                    print(f"    No nucleotide records found for {assembly_acc}")
                    no_16s_genomes.append(assembly_acc)
                    break

                print(f"    Found {len(contig_ids)} contigs, scanning for 16S features...")

                found_16s = False
                # Process in batches
                batch_size = 20
                for batch_start in range(0, len(contig_ids), batch_size):
                    batch_ids = contig_ids[batch_start:batch_start + batch_size]
                    time.sleep(rate_delay)

                    handle = Entrez.efetch(db="nucleotide", id=",".join(batch_ids),
                                           rettype="gb", retmode="text")
                    gb_text = handle.read()
                    handle.close()

                    for gb_record in SeqIO.parse(StringIO(gb_text), "genbank"):
                        # Look for 16S rRNA features
                        for feature in gb_record.features:
                            if feature.type != "rRNA":
                                continue
                            product = feature.qualifiers.get("product", [""])[0]
                            if "16S" not in product:
                                continue

                            seq_16s = feature.extract(gb_record.seq)

                            # Some WGS/RefSeq records use CONTIG joins
                            # and return UndefinedSequence; fall back to
                            # fetching the subsequence via Entrez.
                            try:
                                str(seq_16s)
                            except Exception:
                                loc = feature.location
                                time.sleep(rate_delay)
                                sub_handle = Entrez.efetch(
                                    db="nucleotide", id=gb_record.id,
                                    rettype="fasta", retmode="text",
                                    seq_start=str(int(loc.start) + 1),
                                    seq_stop=str(int(loc.end)),
                                )
                                sub_text = sub_handle.read()
                                sub_handle.close()
                                sub_rec = SeqIO.read(StringIO(sub_text), "fasta")
                                seq_16s = sub_rec.seq
                                if loc.strand == -1:
                                    seq_16s = seq_16s.reverse_complement()

                            # Use assembly accession as the sequence ID
                            seq_id = assembly_acc
                            if found_16s:
                                # Multiple 16S copies: keep the longest
                                existing = [s for s in sequences if s.id == seq_id]
                                if existing and len(seq_16s) <= len(existing[0].seq):
                                    continue
                                # Remove shorter copy
                                sequences = [s for s in sequences if s.id != seq_id]
                                metadata_rows = [m for m in metadata_rows if m['accession'] != seq_id]

                            # Get metadata from contig, with assembly metadata as fallback
                            organism = gb_record.annotations.get("organism", "") or assembly_metadata['organism']
                            source_features = [f for f in gb_record.features if f.type == "source"]

                            host = ""
                            isolation_source = ""
                            strain = ""
                            country = ""

                            if source_features:
                                src = source_features[0]
                                host = src.qualifiers.get("host", [""])[0]
                                isolation_source = src.qualifiers.get("isolation_source", [""])[0]
                                strain = src.qualifiers.get("strain", [""])[0]
                                country = src.qualifiers.get("country", [""])[0]

                            # Use assembly metadata as fallback if contig metadata is empty
                            if not host and assembly_metadata['host']:
                                host = assembly_metadata['host']
                            if not isolation_source and assembly_metadata['isolation_source']:
                                isolation_source = assembly_metadata['isolation_source']
                            if not strain and assembly_metadata['strain']:
                                strain = assembly_metadata['strain']

                            rec = SeqRecord(
                                seq_16s,
                                id=seq_id,
                                description=f"16S rRNA from {organism} ({assembly_acc})"
                            )
                            sequences.append(rec)

                            metadata_rows.append({
                                'accession': seq_id,
                                'organism': organism,
                                'hit_def': f"16S rRNA from {organism}",
                                'host': host,
                                'isolation_source': isolation_source,
                                'strain': strain,
                                'country': country,
                                'identity': '',
                                'alignment_length': '',
                                'evalue': '',
                                'bit_score': '',
                                'source': source,
                            })

                            found_16s = True
                            print(f"    Found 16S: {len(seq_16s)} bp from {gb_record.id}")

                # Fallback: BLAST query 16S against contigs
                if not found_16s and query_fasta:
                    blast_result = _blast_16s_from_contigs(
                        contig_ids, query_fasta, assembly_acc, rate_delay, source,
                        assembly_metadata
                    )
                    if blast_result:
                        seq_rec, meta = blast_result
                        sequences.append(seq_rec)
                        metadata_rows.append(meta)
                        found_16s = True

                if not found_16s:
                    print(f"    No 16S found for {assembly_acc} (annotation + BLAST)")
                    no_16s_genomes.append(assembly_acc)

                break  # success, move to next genome

            except Exception as e:
                err_str = str(e)
                is_transient = any(kw in err_str for kw in _transient)
                if attempt < max_retries - 1 and is_transient:
                    delay = 15 * (2 ** attempt)
                    print(f"    Transient error: {e}")
                    print(f"    Retrying in {delay}s (attempt {attempt + 2}/{max_retries})...")
                    time.sleep(delay)
                else:
                    print(f"    Warning: Failed to process {assembly_acc}: {e}")
                    no_16s_genomes.append(assembly_acc)
                    break

    SeqIO.write(sequences, fasta_output, "fasta")
    print(f"\nWrote {len(sequences)} {source.upper()} 16S sequences to: {fasta_output}")

    metadata_df = pd.DataFrame(metadata_rows)
    metadata_df.to_csv(metadata_output, sep='\t', index=False)
    print(f"Wrote {source.upper()} metadata to: {metadata_output}")

    if no_16s_genomes:
        print(f"\nGenomes without 16S ({len(no_16s_genomes)}):")
        for g in no_16s_genomes:
            print(f"  {g}")


# ---------------------------------------------------------------------------
# Merge sources
# ---------------------------------------------------------------------------

def merge_sources(fasta_paths: list, metadata_paths: list, query_fasta: str,
                  fasta_output: str, metadata_output: str) -> None:
    """Merge FASTA and metadata from multiple sources, deduplicating by accession.

    When the same accession appears in multiple sources, we keep one copy of
    the sequence and combine the source labels (e.g. ``blast,silva``).
    The query sequence is always included and tagged with ``source=query``.
    """
    # Collect all sequences keyed by accession
    seq_by_acc = {}    # accession -> SeqRecord
    source_by_acc = {} # accession -> set of source labels

    # Collect all metadata rows keyed by accession
    meta_by_acc = {}   # accession -> dict (first seen wins for non-source fields)

    for fasta_path, meta_path in zip(fasta_paths, metadata_paths):
        if not Path(fasta_path).exists():
            print(f"  Skipping missing FASTA: {fasta_path}")
            continue

        # Read metadata
        meta_df = pd.DataFrame()
        if Path(meta_path).exists() and Path(meta_path).stat().st_size > 1:
            try:
                meta_df = pd.read_csv(meta_path, sep='\t')
            except pd.errors.EmptyDataError:
                pass

        meta_lookup = {}
        if not meta_df.empty:
            for _, row in meta_df.iterrows():
                acc = str(row.get('accession', ''))
                if acc:
                    meta_lookup[acc] = row.to_dict()

        # Read sequences
        for record in SeqIO.parse(fasta_path, "fasta"):
            acc = record.id
            meta_row = meta_lookup.get(acc, {})
            src = str(meta_row.get('source', 'unknown'))

            if acc not in seq_by_acc:
                seq_by_acc[acc] = record
                source_by_acc[acc] = {src}
                meta_by_acc[acc] = meta_row
            else:
                source_by_acc[acc].add(src)
                # Keep longer sequence
                if len(record.seq) > len(seq_by_acc[acc].seq):
                    seq_by_acc[acc] = record

    # Add query sequence
    if Path(query_fasta).exists():
        query_rec = SeqIO.read(query_fasta, "fasta")
        seq_by_acc[query_rec.id] = query_rec
        source_by_acc[query_rec.id] = {'query'}
        meta_by_acc[query_rec.id] = {
            'accession': query_rec.id,
            'organism': '',
            'hit_def': query_rec.description,
            'host': '',
            'isolation_source': '',
            'strain': '',
            'country': '',
            'source': 'query',
        }

    # Write merged FASTA
    all_records = list(seq_by_acc.values())
    SeqIO.write(all_records, fasta_output, "fasta")

    # Build merged metadata with combined source column
    merged_meta_rows = []
    for acc, meta_row in meta_by_acc.items():
        row = dict(meta_row)
        row['source'] = ','.join(sorted(source_by_acc.get(acc, {'unknown'})))
        merged_meta_rows.append(row)

    merged_df = pd.DataFrame(merged_meta_rows)
    merged_df.to_csv(metadata_output, sep='\t', index=False)

    # Report
    n_total_input = sum(
        sum(1 for _ in SeqIO.parse(fp, "fasta"))
        for fp in fasta_paths if Path(fp).exists()
    )
    # Count non-query merged
    n_merged = len([a for a in seq_by_acc if a not in source_by_acc or 'query' not in source_by_acc[a]])
    n_query = len(all_records) - n_merged

    print(f"Merged sources:")
    print(f"  Input sequences (excl. query): {n_total_input}")
    print(f"  After deduplication: {n_merged} reference + {n_query} query")
    print(f"  Total output: {len(all_records)}")

    # Per-source counts
    source_counts = {}
    for acc, sources in source_by_acc.items():
        if 'query' in sources:
            continue
        for s in sources:
            source_counts[s] = source_counts.get(s, 0) + 1
    for s, c in sorted(source_counts.items()):
        print(f"  {s}: {c} accessions")

    multi = sum(1 for sources in source_by_acc.values() if len(sources) > 1 and 'query' not in sources)
    if multi:
        print(f"  Shared across sources: {multi} accessions")

    print(f"Output FASTA: {fasta_output}")
    print(f"Output metadata: {metadata_output}")


# ---------------------------------------------------------------------------
# Metadata and taxonomy helpers
# ---------------------------------------------------------------------------

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


def get_microbial_lineage_entrez(organism: str, email: str, api_key: str = None,
                                 rate_delay: float = 0.4) -> tuple:
    """Look up NCBI taxonomy for a microbial organism and return formatted lineage and taxon ID.

    Returns a tuple (lineage_string, taxon_id) where lineage_string is like:
      d__Bacteria:2;p__Pseudomonadota:1224;c__Alphaproteobacteria:28211;...

    Returns ("NA", "NA") on failure.
    """
    if not organism:
        return ("NA", "NA")

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Try full organism name, then first two words, then first word
    words = organism.split()
    candidates = [organism]
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
        except Exception:
            continue

    if not tax_id:
        return ("NA", "NA")

    try:
        time.sleep(rate_delay)
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        if not records:
            return ("NA", "NA")

        rec = records[0]
        lineage_ex = rec.get("LineageEx", [])

        # Also include the organism itself if it has a species-level rank
        own_rank = rec.get("Rank", "")
        own_name = rec.get("ScientificName", "")
        own_id = rec.get("TaxId", "")

        rank_map = {}
        for taxon in lineage_ex:
            rank_map[taxon.get("Rank", "")] = (taxon.get("ScientificName", ""), taxon.get("TaxId", ""))

        if own_rank:
            rank_map[own_rank] = (own_name, own_id)

        parts = []
        for rank_name, prefix in MICROBIAL_RANKS:
            if rank_name in rank_map:
                name, tid = rank_map[rank_name]
                parts.append(f"{prefix}{name}:{tid}")
            else:
                parts.append(f"{prefix}:{rank_name}")

        lineage = ";".join(parts)
        # Return the organism's own taxon ID (not the domain/phylum/etc)
        return (lineage, own_id if own_id else tax_id)

    except Exception:
        return ("NA", "NA")


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


def parse_metadata(raw_metadata_path: str, output_path: str, email: str,
                   api_key: str = None, survivors_fasta: str = None) -> None:
    """Parse and annotate sequence metadata.

    If *survivors_fasta* is provided (path to a FASTA file of sequences that
    survived deduplication), only accessions present in that FASTA are kept
    in the output metadata.
    """
    df = pd.read_csv(raw_metadata_path, sep='\t')

    # Filter to survivors if provided
    if survivors_fasta and Path(survivors_fasta).exists():
        survivor_ids = {rec.id for rec in SeqIO.parse(survivors_fasta, "fasta")}
        before = len(df)
        df = df[df['accession'].astype(str).isin(survivor_ids)]
        print(f"Filtering to deduplication survivors: {before} -> {len(df)} sequences")

    rate_delay = 0.1 if api_key else 0.4
    annotated_rows = []

    # Cache for host taxonomy lookups to avoid redundant API calls
    host_taxonomy_cache = {}

    # Cache for microbial lineage lookups
    microbial_lineage_cache = {}

    for idx, row in df.iterrows():
        organism = _str_or_empty(row.get('organism', ''))
        host = _str_or_empty(row.get('host', ''))
        isolation_source = _str_or_empty(row.get('isolation_source', ''))
        source = _str_or_empty(row.get('source', ''))

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
                print(f"  [{idx+1}/{len(df)}] Looking up host taxonomy for: {host_name}")
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

        # Get microbial lineage with caching
        if organism and organism not in microbial_lineage_cache:
            print(f"  [{idx+1}/{len(df)}] Looking up microbial lineage for: {organism}")
            microbial_lineage_cache[organism] = get_microbial_lineage_entrez(
                organism, email, api_key, rate_delay
            )
        lineage_tuple = microbial_lineage_cache.get(organism, ("NA", "NA"))
        # Handle both old string format and new tuple format for backwards compatibility
        if isinstance(lineage_tuple, str):
            microbial_lineage = lineage_tuple
            microbial_taxon_id = "NA"
        else:
            microbial_lineage, microbial_taxon_id = lineage_tuple

        # Clean taxon name
        taxon_name = clean_taxon_name(organism)

        annotated_rows.append({
            'accession': row['accession'],
            'taxon_name': taxon_name,
            'source': source,
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
            'microbial_lineage': microbial_lineage,
            'microbial_taxon_id': microbial_taxon_id,
        })

    annotated_df = pd.DataFrame(annotated_rows)
    annotated_df.to_csv(output_path, sep='\t', index=False)

    # Print summary
    symbiont_count = (annotated_df['symbiont_status'] == 'symbiont').sum()
    host_known_count = (annotated_df['host_name'] != 'NA').sum()
    host_order_known = (annotated_df['host_order'] != 'NA').sum()
    host_phylum_known = (annotated_df['host_phylum'] != 'NA').sum()
    lineage_known = (annotated_df['microbial_lineage'] != 'NA').sum()

    print(f"Metadata annotation complete:")
    print(f"  Total sequences: {len(annotated_df)}")
    print(f"  Symbionts: {symbiont_count}")
    print(f"  Non-symbionts: {len(annotated_df) - symbiont_count}")
    print(f"  With known host: {host_known_count}")
    print(f"  With host order resolved: {host_order_known}")
    print(f"  With host phylum resolved: {host_phylum_known}")
    print(f"  With microbial lineage resolved: {lineage_known}")
    print(f"  Unique hosts looked up: {len(host_taxonomy_cache)}")
    print(f"  Unique organisms looked up: {len(microbial_lineage_cache)}")
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

    # Source counts from merged metadata
    merged_meta_path = output_dir / "all_sources_raw_metadata.tsv"
    if merged_meta_path.exists():
        merged_df = pd.read_csv(merged_meta_path, sep='\t')
        # Count sequences per source (source column may contain comma-separated values)
        source_counts = {}
        for _, row in merged_df.iterrows():
            src = _str_or_empty(row.get('source', ''))
            if not src or src == 'query':
                continue
            for s in src.split(','):
                s = s.strip()
                if s and s != 'query':
                    source_counts[s] = source_counts.get(s, 0) + 1

        total_ref = len(merged_df[merged_df['accession'].astype(str) != query_name])
        report_lines.append(f"- **Total reference sequences (merged)**: {total_ref}")
        for s, c in sorted(source_counts.items()):
            report_lines.append(f"  - from {s}: {c}")
        report_lines.append("")

    # Deduplication stats
    all_seqs_path = output_dir / "all_sequences.fasta"
    dedup_path = output_dir / "deduplicated.fasta"
    if all_seqs_path.exists() and dedup_path.exists():
        n_before = sum(1 for _ in SeqIO.parse(all_seqs_path, "fasta"))
        n_after = sum(1 for _ in SeqIO.parse(dedup_path, "fasta"))
        report_lines.append(f"- **Before cd-hit deduplication**: {n_before} sequences")
        report_lines.append(f"- **After cd-hit deduplication (99% identity)**: {n_after} sequences")
        report_lines.append(f"- **Sequences removed**: {n_before - n_after}")
        report_lines.append("")

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

        # Source breakdown
        if 'source' in df.columns:
            report_lines.extend([
                "## Source Breakdown",
                "",
                "| Source | Count |",
                "|--------|-------|",
            ])
            source_col_counts = df['source'].value_counts()
            for src, cnt in source_col_counts.items():
                report_lines.append(f"| {src} | {cnt} |")
            report_lines.append("")

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

    # Reverse complement stats
    alignment_path = output_dir / "alignment.fasta"
    alignment_raw_path = output_dir / "alignment_raw.fasta"
    if alignment_raw_path.exists():
        rc_count = 0
        for record in SeqIO.parse(alignment_raw_path, "fasta"):
            if record.id.startswith("_R_"):
                rc_count += 1
        if rc_count > 0:
            report_lines.append(f"- **Sequences reverse-complemented by MAFFT**: {rc_count}")
            report_lines.append("")

    # Check for alignment
    if alignment_path.exists():
        from Bio import AlignIO
        alignment = AlignIO.read(str(alignment_path), "fasta")
        report_lines.extend([
            "## Alignment",
            "",
            f"- **Total sequences**: {len(alignment)} (query + references)",
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
        "| `blast_hits.fasta` | Retrieved BLAST sequences |",
        "| `blast_hits_raw_metadata.tsv` | Raw BLAST metadata from NCBI |",
        "| `silva_hits.tsv` | Parsed SILVA nearest-neighbour accessions |",
        "| `silva_hits.fasta` | Retrieved SILVA sequences |",
        "| `silva_hits_raw_metadata.tsv` | Raw SILVA metadata from NCBI |",
        "| `gtdb_wrau01_accessions.txt` | WRAU01 genome accessions from GTDB |",
        "| `gtdb_hits.fasta` | Retrieved GTDB 16S sequences |",
        "| `gtdb_hits_raw_metadata.tsv` | Raw GTDB metadata from NCBI |",
        "| `castelli_accessions.txt` | Castelli et al. 2025 genome accessions |",
        "| `castelli_hits.fasta` | Retrieved Castelli et al. 2025 16S sequences |",
        "| `castelli_hits_raw_metadata.tsv` | Raw Castelli metadata from NCBI |",
        "| `all_sequences.fasta` | All sources merged + query |",
        "| `all_sources_raw_metadata.tsv` | Merged raw metadata |",
        "| `deduplicated.fasta` | After cd-hit-est 99% clustering |",
        "| `sequence_metadata.tsv` | Annotated metadata (survivors only) |",
        "| `alignment.fasta` | MAFFT alignment |",
        "| `16s_tree.treefile` | Best ML tree (Newick format) |",
        "| `16s_tree.iqtree` | IQ-TREE report for best model |",
        "| `16s_tree_rev.iqtree` | IQ-TREE report (reversible models) |",
        "| `16s_tree_nonrev.iqtree` | IQ-TREE report (non-reversible models) |",
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

    # parse-silva command
    silva_parser = subparsers.add_parser('parse-silva', help='Parse SILVA nearest-neighbour FASTA')
    silva_parser.add_argument('--silva-fasta', required=True, help='SILVA FASTA file')
    silva_parser.add_argument('--output', required=True, help='Output TSV path')

    # retrieve-silva-sequences command
    silva_retrieve = subparsers.add_parser('retrieve-silva-sequences', help='Fetch SILVA hit sequences from GenBank')
    silva_retrieve.add_argument('--silva-tsv', required=True, help='TSV from parse-silva')
    silva_retrieve.add_argument('--fasta-output', required=True, help='Output FASTA')
    silva_retrieve.add_argument('--metadata-output', required=True, help='Output metadata TSV')
    silva_retrieve.add_argument('--email', required=True, help='Email for NCBI Entrez')
    silva_retrieve.add_argument('--api-key', help='NCBI API key')

    # parse-gtdb-wrau01 command
    gtdb_parse = subparsers.add_parser('parse-gtdb-wrau01', help='Extract WRAU01 accessions from GTDB MSA')
    gtdb_parse.add_argument('--msa', required=True, help='Path to GTDB MSA FASTA (.gz supported)')
    gtdb_parse.add_argument('--output', required=True, help='Output accessions file')

    # retrieve-gtdb-16s command
    gtdb_retrieve = subparsers.add_parser('retrieve-gtdb-16s', help='Fetch 16S from GTDB WRAU01 genomes')
    gtdb_retrieve.add_argument('--accessions', required=True, help='Accessions file from parse-gtdb-wrau01')
    gtdb_retrieve.add_argument('--fasta-output', required=True, help='Output FASTA')
    gtdb_retrieve.add_argument('--metadata-output', required=True, help='Output metadata TSV')
    gtdb_retrieve.add_argument('--email', required=True, help='Email for NCBI Entrez')
    gtdb_retrieve.add_argument('--api-key', help='NCBI API key')
    gtdb_retrieve.add_argument('--query-fasta', help='Query 16S FASTA for BLAST fallback on unannotated genomes')
    gtdb_retrieve.add_argument('--source', default='gtdb', help='Source label for metadata (default: gtdb)')

    # merge-sources command
    merge_parser = subparsers.add_parser('merge-sources', help='Merge FASTA and metadata from multiple sources')
    merge_parser.add_argument('--fasta-inputs', nargs='+', required=True, help='Input FASTA files')
    merge_parser.add_argument('--metadata-inputs', nargs='+', required=True, help='Input metadata TSV files')
    merge_parser.add_argument('--query-fasta', required=True, help='Query FASTA file')
    merge_parser.add_argument('--fasta-output', required=True, help='Output merged FASTA')
    merge_parser.add_argument('--metadata-output', required=True, help='Output merged metadata TSV')

    # parse-metadata command
    metadata_parser = subparsers.add_parser('parse-metadata', help='Parse and annotate metadata')
    metadata_parser.add_argument('--raw-metadata', required=True, help='Raw metadata TSV')
    metadata_parser.add_argument('--output', required=True, help='Annotated metadata output')
    metadata_parser.add_argument('--email', required=True, help='Email for NCBI Entrez')
    metadata_parser.add_argument('--api-key', help='NCBI API key (optional, increases rate limit)')
    metadata_parser.add_argument('--survivors-fasta', help='FASTA of sequences surviving deduplication')

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

    elif args.command == 'parse-silva':
        parse_silva(args.silva_fasta, args.output)

    elif args.command == 'retrieve-silva-sequences':
        retrieve_silva_sequences(args.silva_tsv, args.fasta_output, args.metadata_output,
                                 args.email, args.api_key)

    elif args.command == 'parse-gtdb-wrau01':
        parse_gtdb_wrau01(args.msa, args.output)

    elif args.command == 'retrieve-gtdb-16s':
        retrieve_gtdb_16s(args.accessions, args.fasta_output, args.metadata_output,
                          args.email, args.api_key, args.query_fasta, args.source)

    elif args.command == 'merge-sources':
        merge_sources(args.fasta_inputs, args.metadata_inputs, args.query_fasta,
                      args.fasta_output, args.metadata_output)

    elif args.command == 'parse-metadata':
        parse_metadata(args.raw_metadata, args.output, args.email, args.api_key,
                       args.survivors_fasta)

    elif args.command == 'find-outgroup':
        outgroup = find_outgroup(args.alignment)
        # Print just the outgroup ID for shell script capture
        print(outgroup)

    elif args.command == 'summary':
        generate_summary(args.output_dir, args.query_name)


if __name__ == '__main__':
    main()
