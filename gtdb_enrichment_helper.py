#!/usr/bin/env python3
"""
Helper script for enriching the phylogeny with GTDB WRAU01 genomes.

Subcommands:
  query-gtdb         Search GTDB API for WRAU01 genomes, build manifest
  add-manual-genome  Download a genome from a WGS accession and add to manifest
  download-proteins  Download protein FASTAs from NCBI for new genomes
  run-eggnog         Run eggNOG-mapper (diamond + emapper) for new genomes
  build-og-mappings  Map proteins to Castelli OGs for new genomes
  download-genomic   Download genomic nucleotide FASTAs (needed for CheckM)
  run-checkm         Run CheckM v1 lineage_wf on genomic FASTAs

Usage:
    python gtdb_enrichment_helper.py query-gtdb --output-dir gtdb_genomes
    python gtdb_enrichment_helper.py add-manual-genome --manifest gtdb_genomes/genome_manifest.json --wgs-accession JBEXBW000000000 --tip-name Symbiont_of_S_maritima ...
    python gtdb_enrichment_helper.py download-proteins --manifest gtdb_genomes/genome_manifest.json
    python gtdb_enrichment_helper.py run-eggnog --manifest gtdb_genomes/genome_manifest.json --eggnog-db eggnog_data --threads 20
    python gtdb_enrichment_helper.py build-og-mappings --manifest gtdb_genomes/genome_manifest.json --castelli-og-dir castelli_et_al/single_ogs
    python gtdb_enrichment_helper.py download-genomic --manifest gtdb_genomes/genome_manifest.json
    python gtdb_enrichment_helper.py run-checkm --manifest gtdb_genomes/genome_manifest.json --threads 20
"""

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path


# ---------------------------------------------------------------------------
# Known existing accessions (from merge_checkm_stats.py TIP_TO_ACCESSION)
# ---------------------------------------------------------------------------

EXISTING_ACCESSIONS = {
    "GCA_002632265.1", "GCA_001510075.1", "GCA_001830425.1",
    "GCA_002327565.1", "GCA_002687515.1", "GCA_009694195.1",
    "GCA_013204045.1", "GCA_014859895.1", "GCA_018662225.1",
    "GCA_009649675.1", "GCA_009780035.1", "GCA_000688235.1",
    "GCA_000300275.1",
}

# Known existing tip names (Castelli-only genomes without GCA accessions)
EXISTING_TIP_NAMES = {
    "Hepatincola_Av", "Hepatincola_Pdp", "Hepatincola_Pp",
    "Tardigradibacter_bertolanii", "s7_ctg000008c",
    "Symbiont_of_S_maritima",
}

# Explicit duplicates: GTDB accession -> existing tip name it duplicates
KNOWN_DUPLICATES = {
    "GCA_023518375.1": "Hepatincola_Av",
}

# Host family mappings for new symbiont genomes
NEW_HOST_FAMILIES = {
    "Pipizella viduata": "Syrphidae",
    "Ecdyonurus torrentis": "Heptageniidae",
}


def _fetch_url(url: str) -> str:
    """Fetch URL content using system curl (works around conda SSL issues)."""
    result = subprocess.run(
        ["/usr/bin/curl", "-s", "-L", url],
        capture_output=True, text=True, check=True,
    )
    return result.stdout


def _fetch_url_binary(url: str, output_path: str) -> None:
    """Download binary file using system curl."""
    subprocess.run(
        ["/usr/bin/curl", "-s", "-L", "-o", output_path, url],
        check=True,
    )


# ---------------------------------------------------------------------------
# query-gtdb
# ---------------------------------------------------------------------------

def _fetch_gtdb_search(order: str) -> list[dict]:
    """Search GTDB API for genomes in an order, return list of genome dicts."""
    url = (
        f"https://gtdb-api.ecogenomic.org/search/gtdb"
        f"?search=o__{order}&page=1&itemsPerPage=100&searchField=all"
    )
    print(f"Querying GTDB API: {url}")
    text = _fetch_url(url)
    data = json.loads(text)

    rows = data.get("rows", [])
    print(f"  Found {len(rows)} genomes in GTDB for o__{order}")
    return rows


def _fetch_gtdb_card(accession: str) -> dict:
    """Fetch detailed card data for a genome from GTDB API."""
    # GTDB API uses accessions with 'RS_' or 'GB_' prefix
    # Try GB_ prefix first (GenBank), then RS_ (RefSeq)
    for prefix in ["GB_", "RS_"]:
        url = f"https://gtdb-api.ecogenomic.org/genome/{prefix}{accession}/card"
        text = _fetch_url(url)
        if text.strip() and not text.startswith("<!"):
            try:
                return json.loads(text)
            except json.JSONDecodeError:
                continue
    print(f"  WARNING: Could not fetch GTDB card for {accession}")
    return {}


def _extract_accession_from_row(row: dict) -> str | None:
    """Extract GCA accession (with version) from a GTDB search row."""
    # The 'gid' field typically has format like 'GB_GCA_031256515.1'
    gid = row.get("gid", "")
    if "GCA_" in gid:
        return gid.split("GCA_", 1)[1]
        # Return full accession
    # Try accession field
    acc = row.get("accession", "")
    if "GCA_" in acc:
        return acc.split("GCA_", 1)[1]
    return None


def _normalize_accession(raw: str) -> str:
    """Normalize to GCA_XXXXXXXXX.X format."""
    if raw.startswith("GCA_"):
        return raw
    # Remove GB_ or RS_ prefix
    for prefix in ["GB_", "RS_"]:
        if raw.startswith(prefix):
            return raw[len(prefix):]
    return raw


def _is_existing(accession: str) -> bool:
    """Check if accession matches an existing genome."""
    return accession in EXISTING_ACCESSIONS


def _is_duplicate(accession: str) -> str | None:
    """Check if genome is a known duplicate. Returns tip name if so."""
    return KNOWN_DUPLICATES.get(accession)


def _generate_tip_name(accession: str, organism: str, existing_names: set[str]) -> str:
    """Generate a clean tip name for a new genome.

    Rules:
    - 'endosymbiont of X' -> 'Symbiont_of_X' (abbreviated genus)
    - Environmental MAGs -> 'WRAU01_MAG_XXXXXXXXX' (using 9-digit accession)
    """
    organism_lower = organism.lower() if organism else ""

    # Endosymbiont pattern
    if "endosymbiont of" in organism_lower:
        # Extract host name
        host = organism.split("endosymbiont of")[-1].strip()
        if host.startswith("*"):
            host = host.strip("*")
        # Abbreviate: "Pipizella viduata" -> "P_viduata"
        parts = host.split()
        if len(parts) >= 2:
            name = f"Symbiont_of_{parts[0][0]}_{parts[1]}"
        else:
            name = f"Symbiont_of_{host.replace(' ', '_')}"
    else:
        # Environmental MAG: use 9-digit accession number
        digits = accession.replace("GCA_", "").split(".")[0]
        name = f"WRAU01_MAG_{digits}"

    # Ensure uniqueness
    if name in existing_names:
        suffix = 2
        while f"{name}_{suffix}" in existing_names:
            suffix += 1
        name = f"{name}_{suffix}"

    return name


def query_gtdb(output_dir: str) -> None:
    """Search GTDB for WRAU01 genomes and build a genome manifest."""
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    manifest_path = out_path / "genome_manifest.json"

    # Fetch GTDB search results
    rows = _fetch_gtdb_search("WRAU01")
    print()

    # Process each genome
    genomes = []
    used_tip_names = set(EXISTING_TIP_NAMES)

    for row in rows:
        gid = row.get("gid", "")
        accession = _normalize_accession(gid)

        if not accession.startswith("GCA_"):
            print(f"  Skipping non-GCA genome: {gid}")
            continue

        print(f"Processing {accession}...")

        # Fetch detailed card data
        card = _fetch_gtdb_card(accession)
        time.sleep(0.3)  # rate limit

        # Extract metadata from card
        metadata = card.get("metadata", {}) if card else {}
        ncbi_organism = (
            metadata.get("ncbi_organism_name", "")
            or row.get("ncbiOrgName", "")
            or ""
        )
        gtdb_taxonomy = row.get("gtdbTaxonomy", "")

        checkm_comp = metadata.get("checkm_completeness", row.get("checkm_completeness", None))
        checkm_cont = metadata.get("checkm_contamination", row.get("checkm_contamination", None))
        gc_raw = metadata.get("gc_percentage", row.get("gc_percentage", None))
        # Normalize GC to fraction (0-1) for consistency with rest of pipeline
        if gc_raw is not None and gc_raw > 1:
            gc_content = gc_raw / 100.0
        else:
            gc_content = gc_raw
        genome_size = metadata.get("genome_size", row.get("genome_size", None))
        n_contigs = metadata.get("contig_count", row.get("contig_count", None))
        biosample_accn = metadata.get("ncbi_biosample", "")

        # Determine status
        dup_of = _is_duplicate(accession)
        if dup_of:
            status = "duplicate"
            tip_name = None
            print(f"  -> DUPLICATE of {dup_of}")
        elif _is_existing(accession):
            status = "existing"
            tip_name = None
            print(f"  -> EXISTING in current tree")
        else:
            status = "new"
            tip_name = _generate_tip_name(accession, ncbi_organism, used_tip_names)
            used_tip_names.add(tip_name)
            print(f"  -> NEW: tip_name = {tip_name}")

        genome_entry = {
            "accession": accession,
            "gid": gid,
            "ncbi_organism": ncbi_organism,
            "gtdb_taxonomy": gtdb_taxonomy,
            "status": status,
            "tip_name": tip_name,
            "duplicate_of": dup_of,
            "checkm_completeness": checkm_comp,
            "checkm_contamination": checkm_cont,
            "gc_content": gc_content,
            "genome_size": genome_size,
            "n_contigs": n_contigs,
            "biosample_accn": biosample_accn,
            "card_metadata": metadata,
        }
        genomes.append(genome_entry)
        print()

    # Summary
    n_existing = sum(1 for g in genomes if g["status"] == "existing")
    n_duplicate = sum(1 for g in genomes if g["status"] == "duplicate")
    n_new = sum(1 for g in genomes if g["status"] == "new")

    print("=" * 60)
    print(f"GTDB WRAU01 genome inventory:")
    print(f"  Total: {len(genomes)}")
    print(f"  Existing: {n_existing}")
    print(f"  Duplicate (skipped): {n_duplicate}")
    print(f"  New: {n_new}")
    print("=" * 60)

    # Write manifest
    manifest = {
        "order": "WRAU01",
        "query_date": time.strftime("%Y-%m-%d"),
        "genomes": genomes,
    }
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"\nManifest written to: {manifest_path}")


# ---------------------------------------------------------------------------
# add-manual-genome
# ---------------------------------------------------------------------------

def _download_wgs_proteins(accession: str, output_dir: Path) -> str | None:
    """Download proteins for a WGS accession (no GCA assembly).

    Three-tier fallback:
      (a) Download pre-computed proteins from NCBI WGS trace archive
      (b) efetch fasta_cds_aa on the WGS master record
      (c) Download nucleotide FASTA + run Prodigal
    """
    faa_path = output_dir / "proteins.faa"
    if faa_path.exists() and faa_path.stat().st_size > 0:
        n_seqs = sum(1 for line in open(faa_path) if line.startswith(">"))
        print(f"  Protein FASTA already exists: {faa_path} ({n_seqs} proteins)")
        return str(faa_path)

    # Derive WGS prefix from accession (e.g. JBEXBW000000000 -> JBEXBW01)
    # WGS master accessions: 6-letter prefix + 9 zeros; version is 01
    wgs_prefix = accession[:6] + "01"
    gz_path = output_dir / "proteins.faa.gz"

    # (a) Try NCBI WGS trace archive (pre-computed protein FASTA)
    # URL pattern: https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/JB/EX/BW/JBEXBW01/JBEXBW01.1.fsa_aa.gz
    p = wgs_prefix  # e.g. JBEXBW01
    wgs_url = (
        f"https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/"
        f"{p[0:2]}/{p[2:4]}/{p[4:6]}/{p}/{p}.1.fsa_aa.gz"
    )
    print(f"  Trying WGS trace archive: {wgs_url}")
    _fetch_url_binary(wgs_url, str(gz_path))

    if gz_path.exists() and gz_path.stat().st_size > 100:
        try:
            with gzip.open(str(gz_path), 'rb') as f_in:
                with open(str(faa_path), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            gz_path.unlink()
            with open(faa_path) as fh:
                first_line = fh.readline()
            if first_line.startswith(">"):
                n_seqs = sum(1 for line in open(faa_path) if line.startswith(">"))
                print(f"  Downloaded {n_seqs} proteins from WGS trace archive")
                return str(faa_path)
            else:
                faa_path.unlink(missing_ok=True)
        except gzip.BadGzipFile:
            gz_path.unlink(missing_ok=True)
            faa_path.unlink(missing_ok=True)
    else:
        gz_path.unlink(missing_ok=True)

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # (b) Try efetch fasta_cds_aa directly on WGS master
    print(f"  Trying efetch fasta_cds_aa on WGS master {accession}...")
    url = f"{base_url}/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta_cds_aa&retmode=text"
    text = _fetch_url(url)
    time.sleep(0.5)
    if text.strip() and text.strip().startswith(">"):
        with open(faa_path, "w") as fh:
            fh.write(text)
        n_seqs = sum(1 for line in open(faa_path) if line.startswith(">"))
        if n_seqs > 10:
            print(f"  Downloaded {n_seqs} proteins via efetch (WGS master)")
            return str(faa_path)
        else:
            faa_path.unlink(missing_ok=True)

    # (c) Download nucleotide FASTA + Prodigal
    print(f"  Falling back to nucleotide download + Prodigal...")
    fna_path = output_dir / "genomic.fna"

    if not (fna_path.exists() and fna_path.stat().st_size > 0):
        # Try WGS trace archive for nucleotide
        nuc_gz_path = output_dir / "genomic.fna.gz"
        nuc_url = (
            f"https://sra-download.ncbi.nlm.nih.gov/traces/wgs04/wgs_aux/"
            f"{p[0:2]}/{p[2:4]}/{p[4:6]}/{p}/{p}.1.fsa_nt.gz"
        )
        print(f"  Trying WGS trace nucleotide: {nuc_url}")
        _fetch_url_binary(nuc_url, str(nuc_gz_path))

        if nuc_gz_path.exists() and nuc_gz_path.stat().st_size > 100:
            try:
                with gzip.open(str(nuc_gz_path), 'rb') as f_in:
                    with open(str(fna_path), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                nuc_gz_path.unlink()
                print(f"  Downloaded nucleotide FASTA ({fna_path.stat().st_size:,} bytes)")
            except gzip.BadGzipFile:
                nuc_gz_path.unlink(missing_ok=True)
        else:
            nuc_gz_path.unlink(missing_ok=True)

    if not (fna_path.exists() and fna_path.stat().st_size > 0):
        print(f"  WARNING: Could not download nucleotide FASTA")
        return None

    # Run Prodigal
    prodigal_bin = shutil.which("prodigal")
    if not prodigal_bin:
        prodigal_bin = os.path.expanduser("~/miniconda3/envs/checkm2/bin/prodigal")
    if not os.path.exists(prodigal_bin):
        print(f"  WARNING: prodigal not found")
        return None

    print(f"  Running prodigal -p meta...")
    try:
        subprocess.run(
            [prodigal_bin, "-i", str(fna_path), "-a", str(faa_path),
             "-p", "meta", "-q"],
            capture_output=True, text=True, check=True,
        )
        n_seqs = sum(1 for line in open(faa_path) if line.startswith(">"))
        print(f"  Prodigal predicted {n_seqs} proteins")
        return str(faa_path)
    except subprocess.CalledProcessError as e:
        print(f"  WARNING: Prodigal failed: {e.stderr[:200]}")
        return None


def add_manual_genome(manifest_path: str, wgs_accession: str, tip_name: str,
                      organism_name: str, eggnog_db: str, castelli_og_dir: str,
                      threads: int, checkm_completeness: float | None = None,
                      checkm_contamination: float | None = None,
                      gc_content: float | None = None,
                      genome_size: int | None = None,
                      n_contigs: int | None = None) -> None:
    """Download a genome from a WGS accession and integrate it into the manifest.

    This handles genomes that don't have a GCA assembly accession (e.g. the
    Strigamia maritima endosymbiont which only has WGS accession JBEXBW000000000).
    """
    with open(manifest_path) as fh:
        manifest = json.load(fh)

    manifest_dir = Path(manifest_path).parent
    genome_dir = manifest_dir / wgs_accession
    genome_dir.mkdir(parents=True, exist_ok=True)

    # Check if already fully processed
    faa_path = genome_dir / "proteins.faa"
    eggnog_annotations = str(genome_dir / "eggnog" / "eggnog.emapper.annotations")
    og_mapping_file = str(genome_dir / "og_mapping.tsv")

    for genome in manifest["genomes"]:
        if genome.get("tip_name") == tip_name:
            if (os.path.exists(og_mapping_file) and os.path.getsize(og_mapping_file) > 0
                    and faa_path.exists() and faa_path.stat().st_size > 0):
                print(f"Genome {tip_name} already fully processed, ensuring manifest paths are correct.")
                # Repair manifest entry paths in case a previous step clobbered them
                genome["protein_fasta"] = str(faa_path)
                genome["eggnog_annotations"] = eggnog_annotations
                genome["og_mapping"] = og_mapping_file
                with open(manifest_path, "w") as fh:
                    json.dump(manifest, fh, indent=2)
                return
            break
    else:
        genome = None

    # Step 1: Download proteins
    print(f"=== Adding manual genome: {tip_name} ({wgs_accession}) ===")
    protein_path = _download_wgs_proteins(wgs_accession, genome_dir)
    if not protein_path:
        print(f"ERROR: Could not obtain proteins for {wgs_accession}")
        return

    # Step 2: Run eggNOG-mapper
    eggnog_dir = genome_dir / "eggnog"
    eggnog_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = str(eggnog_dir / "eggnog")
    annotations_file = f"{output_prefix}.emapper.annotations"

    if os.path.exists(annotations_file) and os.path.getsize(annotations_file) > 0:
        print(f"  eggNOG annotations already exist: {annotations_file}")
    else:
        diamond_bin = _conda_bin("diamond")
        emapper_bin = _conda_bin("emapper.py")
        print(f"  Running diamond blastp...")
        diamond_tmp = f"{output_prefix}.diamond.tmp"
        seed_orthologs = f"{output_prefix}.seed_orthologs"
        try:
            subprocess.run(
                [diamond_bin, "blastp",
                 "-d", f"{eggnog_db}/eggnog_proteins.dmnd",
                 "-q", protein_path,
                 "--threads", str(threads),
                 "-o", diamond_tmp,
                 "-e", "0.001",
                 "--max-target-seqs", "10",
                 "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                 "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                 "evalue", "bitscore", "qcovhsp", "scovhsp"],
                capture_output=True, text=True, check=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"  diamond failed: {e.stderr[:300]}")
            return

        with open(diamond_tmp) as fin, open(seed_orthologs, "w") as fout:
            for line in fin:
                fields = line.strip().split("\t")
                if len(fields) >= 12:
                    fout.write(f"{fields[0]}\t{fields[1]}\t{fields[10]}\t{fields[11]}\n")
        os.unlink(diamond_tmp)

        print(f"  Running emapper.py --no_search...")
        try:
            subprocess.run(
                [emapper_bin,
                 "-m", "no_search",
                 "--annotate_hits_table", seed_orthologs,
                 "--no_file_comments",
                 "--output", output_prefix,
                 "--data_dir", eggnog_db,
                 "--cpu", str(threads)],
                capture_output=True, text=True, check=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"  emapper failed: {e.stderr[:300]}")
            return

    if not os.path.exists(annotations_file):
        print(f"  WARNING: eggNOG annotations not produced")
        return

    n_annotated = sum(1 for line in open(annotations_file)
                      if not line.startswith("#") and line.strip())
    print(f"  Annotated {n_annotated} proteins")

    # Step 3: Build OG mapping
    sys.path.insert(0, str(Path(manifest_path).resolve().parent.parent))
    from phylo_genome_helper import map_ogs

    og_mapping_path = str(genome_dir / "og_mapping.tsv")
    if os.path.exists(og_mapping_path) and os.path.getsize(og_mapping_path) > 0:
        print(f"  OG mapping already exists: {og_mapping_path}")
    else:
        map_ogs(annotations_file, castelli_og_dir, og_mapping_path)

    n_mapped = sum(1 for i, line in enumerate(open(og_mapping_path))
                   if i > 0 and line.strip())
    print(f"  Mapped to {n_mapped} Castelli OGs")

    # Step 4: Add/update manifest entry
    entry = {
        "accession": wgs_accession,
        "gid": wgs_accession,
        "ncbi_organism": organism_name,
        "gtdb_taxonomy": "",
        "status": "new",
        "tip_name": tip_name,
        "duplicate_of": None,
        "checkm_completeness": checkm_completeness,
        "checkm_contamination": checkm_contamination,
        "gc_content": gc_content,
        "genome_size": genome_size,
        "n_contigs": n_contigs,
        "biosample_accn": "",
        "protein_fasta": protein_path,
        "eggnog_annotations": annotations_file,
        "og_mapping": og_mapping_path,
    }

    # Replace existing entry or append
    replaced = False
    for i, g in enumerate(manifest["genomes"]):
        if g.get("tip_name") == tip_name:
            manifest["genomes"][i] = entry
            replaced = True
            break
    if not replaced:
        manifest["genomes"].append(entry)

    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\n  Manifest updated with {tip_name}")
    print(f"  Proteins: {protein_path}")
    print(f"  OG mapping: {og_mapping_path} ({n_mapped} OGs)")


# ---------------------------------------------------------------------------
# download-proteins
# ---------------------------------------------------------------------------

def _get_ncbi_ftp_path(accession: str) -> str | None:
    """Query NCBI Entrez for the GenBank FTP path of an assembly."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    term = urllib.parse.quote(f"{accession}[Assembly Accession]")
    search_url = f"{base_url}/esearch.fcgi?db=assembly&term={term}&retmode=xml"

    xml_text = _fetch_url(search_url)
    tree = ET.fromstring(xml_text)
    id_list = tree.findall(".//Id")
    if not id_list:
        return None

    assembly_id = id_list[0].text
    time.sleep(0.35)

    summary_url = f"{base_url}/esummary.fcgi?db=assembly&id={assembly_id}&retmode=xml"
    xml_text = _fetch_url(summary_url)
    tree = ET.fromstring(xml_text)

    # Look for FtpPath_GenBank
    ftp_elem = tree.find(".//FtpPath_GenBank")
    if ftp_elem is not None and ftp_elem.text:
        return ftp_elem.text

    # Fallback: FtpPath_RefSeq
    ftp_elem = tree.find(".//FtpPath_RefSeq")
    if ftp_elem is not None and ftp_elem.text:
        return ftp_elem.text

    return None


def _download_protein_fasta(accession: str, ftp_path: str, output_dir: Path) -> str | None:
    """Download protein FASTA from NCBI FTP. Returns path or None."""
    assembly_name = ftp_path.rstrip("/").split("/")[-1]
    protein_url = f"{ftp_path}/{assembly_name}_protein.faa.gz"

    # Convert ftp:// to https://
    protein_url = protein_url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")

    gz_path = output_dir / "proteins.faa.gz"
    faa_path = output_dir / "proteins.faa"

    if faa_path.exists() and faa_path.stat().st_size > 0:
        print(f"    Protein FASTA already exists: {faa_path}")
        return str(faa_path)

    print(f"    Downloading: {protein_url}")
    _fetch_url_binary(protein_url, str(gz_path))

    if gz_path.exists() and gz_path.stat().st_size > 100:
        # Decompress
        try:
            with gzip.open(str(gz_path), 'rb') as f_in:
                with open(str(faa_path), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            gz_path.unlink()
        except gzip.BadGzipFile:
            print(f"    WARNING: Downloaded file is not a valid gzip")
            gz_path.unlink(missing_ok=True)
            faa_path.unlink(missing_ok=True)
            return None

        # Verify it's a valid FASTA
        with open(faa_path) as fh:
            first_line = fh.readline()
        if first_line.startswith(">"):
            n_seqs = sum(1 for line in open(faa_path) if line.startswith(">"))
            print(f"    Downloaded {n_seqs} proteins")
            return str(faa_path)
        else:
            print(f"    WARNING: Downloaded file is not a valid FASTA")
            faa_path.unlink(missing_ok=True)
    else:
        gz_path.unlink(missing_ok=True)

    return None


def _run_prodigal(accession: str, ftp_path: str, output_dir: Path) -> str | None:
    """Download genomic FASTA and run prodigal to predict proteins."""
    assembly_name = ftp_path.rstrip("/").split("/")[-1]
    genomic_url = f"{ftp_path}/{assembly_name}_genomic.fna.gz"
    genomic_url = genomic_url.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")

    gz_path = output_dir / "genomic.fna.gz"
    fna_path = output_dir / "genomic.fna"
    faa_path = output_dir / "proteins.faa"

    # Skip download if genomic FASTA already exists
    if not (fna_path.exists() and fna_path.stat().st_size > 0):
        print(f"    Downloading genomic FASTA: {genomic_url}")
        _fetch_url_binary(genomic_url, str(gz_path))

        if not gz_path.exists() or gz_path.stat().st_size < 100:
            print(f"    WARNING: Failed to download genomic FASTA")
            gz_path.unlink(missing_ok=True)
            return None

        try:
            with gzip.open(str(gz_path), 'rb') as f_in:
                with open(str(fna_path), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            gz_path.unlink()
        except gzip.BadGzipFile:
            print(f"    WARNING: Genomic FASTA is not a valid gzip")
            gz_path.unlink(missing_ok=True)
            return None

    # Find prodigal: try PATH first, then known conda env location
    prodigal_bin = shutil.which("prodigal")
    if not prodigal_bin:
        prodigal_bin = os.path.expanduser("~/miniconda3/envs/checkm2/bin/prodigal")
    if not os.path.exists(prodigal_bin):
        print(f"    WARNING: prodigal not found")
        fna_path.unlink(missing_ok=True)
        return None

    print(f"    Running prodigal -p meta...")
    try:
        subprocess.run(
            [prodigal_bin, "-i", str(fna_path), "-a", str(faa_path),
             "-p", "meta", "-q"],
            capture_output=True, text=True, check=True,
        )
        # Clean up genomic file
        fna_path.unlink(missing_ok=True)

        n_seqs = sum(1 for line in open(faa_path) if line.startswith(">"))
        print(f"    Prodigal predicted {n_seqs} proteins")
        return str(faa_path)
    except subprocess.CalledProcessError as e:
        print(f"    WARNING: Prodigal failed: {e.stderr[:200]}")
        fna_path.unlink(missing_ok=True)
        return None


def download_proteins(manifest_path: str) -> None:
    """Download protein FASTAs for new genomes in the manifest."""
    with open(manifest_path) as fh:
        manifest = json.load(fh)

    new_genomes = [g for g in manifest["genomes"] if g["status"] == "new"]
    print(f"Downloading proteins for {len(new_genomes)} new genomes...\n")

    for genome in new_genomes:
        accession = genome["accession"]
        tip_name = genome["tip_name"]
        output_dir = Path(manifest_path).parent / accession
        output_dir.mkdir(parents=True, exist_ok=True)

        # Skip if protein FASTA already exists (e.g. from add-manual-genome)
        existing_faa = genome.get("protein_fasta")
        faa_on_disk = output_dir / "proteins.faa"
        if existing_faa and os.path.exists(existing_faa) and os.path.getsize(existing_faa) > 0:
            print(f"[{accession}] {tip_name}: proteins already exist, skipping.")
            continue
        if faa_on_disk.exists() and faa_on_disk.stat().st_size > 0:
            print(f"[{accession}] {tip_name}: proteins found on disk, updating manifest.")
            genome["protein_fasta"] = str(faa_on_disk)
            continue

        print(f"[{accession}] {tip_name}")

        # Get FTP path from NCBI
        ftp_path = _get_ncbi_ftp_path(accession)
        time.sleep(0.35)

        if not ftp_path:
            print(f"    WARNING: No FTP path found for {accession}")
            genome["protein_fasta"] = None
            continue

        genome["ftp_path"] = ftp_path

        # Try downloading protein FASTA
        protein_path = _download_protein_fasta(accession, ftp_path, output_dir)

        if not protein_path:
            # Fallback: prodigal
            print(f"    Protein FASTA not available, trying prodigal...")
            protein_path = _run_prodigal(accession, ftp_path, output_dir)

        genome["protein_fasta"] = protein_path
        print()

    # Update manifest
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    n_success = sum(1 for g in new_genomes if g.get("protein_fasta"))
    n_fail = sum(1 for g in new_genomes if not g.get("protein_fasta"))
    print(f"\nDownload summary: {n_success} succeeded, {n_fail} failed")


# ---------------------------------------------------------------------------
# run-eggnog
# ---------------------------------------------------------------------------

def _conda_bin(tool: str) -> str:
    """Resolve a tool in the current Python interpreter's conda env bin dir."""
    env_bin = str(Path(sys.executable).parent)
    candidate = os.path.join(env_bin, tool)
    if os.path.exists(candidate):
        return candidate
    # Fallback: try PATH
    found = shutil.which(tool)
    if found:
        return found
    return tool  # hope for the best


def run_eggnog(manifest_path: str, eggnog_db: str, threads: int) -> None:
    """Run eggNOG-mapper for new genomes."""
    with open(manifest_path) as fh:
        manifest = json.load(fh)

    diamond_bin = _conda_bin("diamond")
    emapper_bin = _conda_bin("emapper.py")
    print(f"Using diamond: {diamond_bin}")
    print(f"Using emapper: {emapper_bin}")

    new_genomes = [g for g in manifest["genomes"]
                   if g["status"] == "new" and g.get("protein_fasta")]
    print(f"Running eggNOG-mapper for {len(new_genomes)} new genomes...\n")

    for genome in new_genomes:
        accession = genome["accession"]
        tip_name = genome["tip_name"]
        protein_fasta = genome["protein_fasta"]
        eggnog_dir = Path(manifest_path).parent / accession / "eggnog"
        eggnog_dir.mkdir(parents=True, exist_ok=True)
        output_prefix = str(eggnog_dir / "eggnog")

        annotations_file = f"{output_prefix}.emapper.annotations"
        if os.path.exists(annotations_file) and os.path.getsize(annotations_file) > 0:
            print(f"[{accession}] {tip_name}: eggNOG output already exists, skipping.")
            genome["eggnog_annotations"] = annotations_file
            continue

        print(f"[{accession}] {tip_name}")

        # Step 1: diamond blastp
        diamond_tmp = f"{output_prefix}.diamond.tmp"
        seed_orthologs = f"{output_prefix}.seed_orthologs"
        print(f"  Running diamond blastp...")
        try:
            subprocess.run(
                [diamond_bin, "blastp",
                 "-d", f"{eggnog_db}/eggnog_proteins.dmnd",
                 "-q", protein_fasta,
                 "--threads", str(threads),
                 "-o", diamond_tmp,
                 "-e", "0.001",
                 "--max-target-seqs", "10",
                 "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                 "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                 "evalue", "bitscore", "qcovhsp", "scovhsp"],
                capture_output=True, text=True, check=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"  diamond failed: {e.stderr[:300]}")
            continue

        # Reformat to seed_orthologs: query, hit, evalue, score
        with open(diamond_tmp) as fin, open(seed_orthologs, "w") as fout:
            for line in fin:
                fields = line.strip().split("\t")
                if len(fields) >= 12:
                    fout.write(f"{fields[0]}\t{fields[1]}\t{fields[10]}\t{fields[11]}\n")
        os.unlink(diamond_tmp)

        # Step 2: emapper annotation
        print(f"  Running emapper.py --no_search...")
        try:
            subprocess.run(
                [emapper_bin,
                 "-m", "no_search",
                 "--annotate_hits_table", seed_orthologs,
                 "--no_file_comments",
                 "--output", output_prefix,
                 "--data_dir", eggnog_db,
                 "--cpu", str(threads)],
                capture_output=True, text=True, check=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"  emapper failed: {e.stderr[:300]}")
            continue

        if os.path.exists(annotations_file):
            n_annotated = sum(1 for line in open(annotations_file)
                              if not line.startswith("#") and line.strip())
            print(f"  Annotated {n_annotated} proteins")
            genome["eggnog_annotations"] = annotations_file
        else:
            print(f"  WARNING: Annotations file not produced")

        print()

    # Update manifest
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    n_annotated = sum(1 for g in new_genomes if g.get("eggnog_annotations"))
    print(f"\neggNOG summary: {n_annotated}/{len(new_genomes)} genomes annotated")


# ---------------------------------------------------------------------------
# build-og-mappings
# ---------------------------------------------------------------------------

def build_og_mappings(manifest_path: str, castelli_og_dir: str) -> None:
    """Map proteins to Castelli OGs for new genomes, reusing phylo_genome_helper.map_ogs."""
    # Import map_ogs from phylo_genome_helper
    sys.path.insert(0, str(Path(manifest_path).resolve().parent.parent))
    from phylo_genome_helper import map_ogs

    with open(manifest_path) as fh:
        manifest = json.load(fh)

    new_genomes = [g for g in manifest["genomes"]
                   if g["status"] == "new" and g.get("eggnog_annotations")]
    print(f"Building OG mappings for {len(new_genomes)} new genomes...\n")

    for genome in new_genomes:
        accession = genome["accession"]
        tip_name = genome["tip_name"]
        annotations = genome["eggnog_annotations"]
        genome_dir = Path(manifest_path).parent / accession
        output_file = str(genome_dir / "og_mapping.tsv")

        print(f"[{accession}] {tip_name}")
        map_ogs(annotations, castelli_og_dir, output_file)
        genome["og_mapping"] = output_file

        # Count mapped OGs
        n_mapped = sum(1 for i, line in enumerate(open(output_file))
                       if i > 0 and line.strip())
        print(f"  Mapped to {n_mapped} Castelli OGs")
        print()

    # Update manifest
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    n_mapped = sum(1 for g in new_genomes if g.get("og_mapping"))
    print(f"\nOG mapping summary: {n_mapped}/{len(new_genomes)} genomes mapped")


# ---------------------------------------------------------------------------
# download-genomic
# ---------------------------------------------------------------------------

def download_genomic(manifest_path: str) -> None:
    """Download genomic nucleotide FASTAs for new genomes (needed for CheckM2)."""
    with open(manifest_path) as fh:
        manifest = json.load(fh)

    new_genomes = [g for g in manifest["genomes"]
                   if g["status"] == "new" and g.get("ftp_path")]
    print(f"Downloading genomic FASTAs for {len(new_genomes)} new genomes...\n")

    for genome in new_genomes:
        accession = genome["accession"]
        tip_name = genome["tip_name"]
        ftp_path = genome["ftp_path"]
        output_dir = Path(manifest_path).parent / accession
        output_dir.mkdir(parents=True, exist_ok=True)

        fna_path = output_dir / "genomic.fna"

        if fna_path.exists() and fna_path.stat().st_size > 0:
            print(f"[{accession}] {tip_name}: genomic.fna already exists")
            genome["genomic_fasta"] = str(fna_path)
            continue

        assembly_name = ftp_path.rstrip("/").split("/")[-1]
        genomic_url = f"{ftp_path}/{assembly_name}_genomic.fna.gz"
        genomic_url = genomic_url.replace(
            "ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov"
        )
        gz_path = output_dir / "genomic.fna.gz"

        print(f"[{accession}] {tip_name}")
        print(f"  Downloading: {genomic_url}")
        _fetch_url_binary(genomic_url, str(gz_path))

        if not gz_path.exists() or gz_path.stat().st_size < 100:
            print(f"  WARNING: Download failed")
            gz_path.unlink(missing_ok=True)
            continue

        try:
            with gzip.open(str(gz_path), 'rb') as f_in:
                with open(str(fna_path), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            gz_path.unlink()
            print(f"  Downloaded: {fna_path.stat().st_size:,} bytes")
            genome["genomic_fasta"] = str(fna_path)
        except gzip.BadGzipFile:
            print(f"  WARNING: Not a valid gzip file")
            gz_path.unlink(missing_ok=True)

        print()

    # Update manifest
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    n_success = sum(1 for g in new_genomes if g.get("genomic_fasta"))
    print(f"\nDownload summary: {n_success}/{len(new_genomes)} genomic FASTAs")


# ---------------------------------------------------------------------------
# run-checkm2
# ---------------------------------------------------------------------------

def run_checkm(manifest_path: str, threads: int) -> None:
    """Run CheckM v1 (lineage_wf) on genomic FASTAs for new genomes.

    Uses the same workflow as 06_genome_phylogenetic_analysis.sh to ensure
    comparable completeness/contamination estimates with Castelli et al.
    """
    with open(manifest_path) as fh:
        manifest = json.load(fh)

    new_genomes = [g for g in manifest["genomes"]
                   if g["status"] == "new" and g.get("genomic_fasta")]
    if not new_genomes:
        print("No genomes with genomic FASTAs to run CheckM on.")
        return

    # Prepare input directory with symlinks
    checkm_input = Path(manifest_path).parent / "checkm_input"
    checkm_output = Path(manifest_path).parent / "checkm_output"
    checkm_input.mkdir(parents=True, exist_ok=True)

    for genome in new_genomes:
        accession = genome["accession"]
        fna_path = Path(genome["genomic_fasta"])
        link_path = checkm_input / f"{accession}.fna"
        if link_path.exists():
            link_path.unlink()
        link_path.symlink_to(fna_path.resolve())

    # Find CheckM v1 binary
    checkm_bin = os.path.expanduser("~/miniconda3/envs/checkm/bin/checkm")
    if not os.path.exists(checkm_bin):
        checkm_bin = shutil.which("checkm")
    if not checkm_bin or not os.path.exists(checkm_bin):
        print("ERROR: CheckM v1 not found")
        return

    quality_report = checkm_output / "quality_report.tsv"
    lineage_ms = checkm_output / "lineage.ms"

    # Build env with CheckM's conda env bin dir on PATH so it finds prodigal, hmmer, etc.
    checkm_env_bin = str(Path(checkm_bin).parent)
    env = dict(os.environ)
    env["PATH"] = checkm_env_bin + ":" + env.get("PATH", "")

    # Step 1: Run lineage_wf if not already done
    if lineage_ms.exists():
        print(f"CheckM lineage_wf already completed (found {lineage_ms})")
    else:
        print(f"Running CheckM v1 lineage_wf on {len(new_genomes)} genomes...")
        print(f"  Using: {checkm_bin}")
        print(f"  PATH prepend: {checkm_env_bin}")

        try:
            result = subprocess.run(
                [checkm_bin, "lineage_wf",
                 "-x", "fna",
                 "-t", str(threads),
                 "--reduced_tree",
                 str(checkm_input),
                 str(checkm_output)],
                capture_output=True, text=True, check=True,
                env=env,
            )
            print(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"CheckM lineage_wf failed: {e.stderr[-500:]}")
            return

    # Step 2: Run checkm qa with extended output (-o 2) to get GC, genome size, etc.
    if quality_report.exists() and quality_report.stat().st_size > 0:
        print(f"CheckM quality report already exists: {quality_report}")
    else:
        print(f"Running CheckM qa with extended output...")
        try:
            result = subprocess.run(
                [checkm_bin, "qa",
                 "-o", "2",
                 "--tab_table",
                 "-f", str(quality_report),
                 "-t", str(threads),
                 str(lineage_ms),
                 str(checkm_output)],
                capture_output=True, text=True, check=True,
                env=env,
            )
            print(result.stdout[-500:] if len(result.stdout) > 500 else result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"CheckM qa failed: {e.stderr[-500:]}")
            return

    # Parse CheckM v1 output
    if not quality_report.exists():
        print("WARNING: CheckM quality_report.tsv not found")
        return

    import csv
    with open(quality_report) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        checkm_results = {row['Bin Id']: row for row in reader}

    print(f"\nCheckM v1 results:")
    for genome in new_genomes:
        accession = genome["accession"]
        if accession in checkm_results:
            row = checkm_results[accession]
            comp = float(row['Completeness'])
            cont = float(row['Contamination'])
            gc = float(row['GC']) / 100.0  # Convert % to fraction
            genome_size = int(row['Genome size (bp)'])
            n_contigs = int(row['# contigs'])
            genome["checkm_completeness"] = comp
            genome["checkm_contamination"] = cont
            genome["gc_content"] = gc
            genome["genome_size"] = genome_size
            genome["n_contigs"] = n_contigs
            print(f"  {accession} ({genome['tip_name']}): "
                  f"comp={comp:.1f}%, cont={cont:.1f}%, GC={gc:.3f}, "
                  f"size={genome_size:,}bp, contigs={n_contigs}")
        else:
            print(f"  {accession}: not found in CheckM output")

    # Update manifest
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\nCheckM v1 complete. Results saved to manifest.")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Enrich phylogeny with GTDB WRAU01 genomes"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # query-gtdb
    p = subparsers.add_parser("query-gtdb", help="Search GTDB for WRAU01 genomes")
    p.add_argument("--output-dir", required=True,
                   help="Output directory for manifest and genome data")

    # add-manual-genome
    p = subparsers.add_parser("add-manual-genome",
                              help="Download a WGS genome and add to manifest")
    p.add_argument("--manifest", required=True, help="Path to genome_manifest.json")
    p.add_argument("--wgs-accession", required=True, help="WGS master accession (e.g. JBEXBW000000000)")
    p.add_argument("--tip-name", required=True, help="Tip name for phylogeny")
    p.add_argument("--organism-name", default="", help="Organism name")
    p.add_argument("--eggnog-db", required=True, help="Path to eggNOG database directory")
    p.add_argument("--castelli-og-dir", required=True, help="Directory with Castelli single-copy OG files")
    p.add_argument("--threads", type=int, default=4, help="Number of threads")
    p.add_argument("--checkm-completeness", type=float, default=None, help="Pre-computed CheckM completeness")
    p.add_argument("--checkm-contamination", type=float, default=None, help="Pre-computed CheckM contamination")
    p.add_argument("--gc-content", type=float, default=None, help="Pre-computed GC content (fraction 0-1)")
    p.add_argument("--genome-size", type=int, default=None, help="Pre-computed genome size (bp)")
    p.add_argument("--n-contigs", type=int, default=None, help="Pre-computed number of contigs")

    # download-proteins
    p = subparsers.add_parser("download-proteins",
                              help="Download protein FASTAs for new genomes")
    p.add_argument("--manifest", required=True, help="Path to genome_manifest.json")

    # run-eggnog
    p = subparsers.add_parser("run-eggnog",
                              help="Run eggNOG-mapper for new genomes")
    p.add_argument("--manifest", required=True, help="Path to genome_manifest.json")
    p.add_argument("--eggnog-db", required=True, help="Path to eggNOG database directory")
    p.add_argument("--threads", type=int, default=4, help="Number of threads")

    # build-og-mappings
    p = subparsers.add_parser("build-og-mappings",
                              help="Map proteins to Castelli OGs")
    p.add_argument("--manifest", required=True, help="Path to genome_manifest.json")
    p.add_argument("--castelli-og-dir", required=True,
                   help="Directory with Castelli single-copy OG files")

    # download-genomic
    p = subparsers.add_parser("download-genomic",
                              help="Download genomic nucleotide FASTAs for new genomes")
    p.add_argument("--manifest", required=True, help="Path to genome_manifest.json")

    # run-checkm
    p = subparsers.add_parser("run-checkm",
                              help="Run CheckM v1 on genomic FASTAs for new genomes")
    p.add_argument("--manifest", required=True, help="Path to genome_manifest.json")
    p.add_argument("--threads", type=int, default=4, help="Number of threads")

    args = parser.parse_args()

    if args.command == "query-gtdb":
        query_gtdb(args.output_dir)
    elif args.command == "add-manual-genome":
        add_manual_genome(
            args.manifest, args.wgs_accession, args.tip_name,
            args.organism_name, args.eggnog_db, args.castelli_og_dir,
            args.threads, args.checkm_completeness, args.checkm_contamination,
            args.gc_content, args.genome_size, args.n_contigs,
        )
    elif args.command == "download-proteins":
        download_proteins(args.manifest)
    elif args.command == "run-eggnog":
        run_eggnog(args.manifest, args.eggnog_db, args.threads)
    elif args.command == "build-og-mappings":
        build_og_mappings(args.manifest, args.castelli_og_dir)
    elif args.command == "download-genomic":
        download_genomic(args.manifest)
    elif args.command == "run-checkm":
        run_checkm(args.manifest, args.threads)


if __name__ == "__main__":
    main()
