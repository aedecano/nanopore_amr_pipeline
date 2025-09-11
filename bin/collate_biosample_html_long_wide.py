#!/usr/bin/env python3
"""
Collate NCBI BioSample esummary XML files (saved as .html) into wide and long tables,
and extract 'sample_source' and 'collection_year'.

Outputs:
- Wide TSV: one row per BioSample (core fields + attributes + sample_source + collection_year)
- Long TSV: one row per (BioSample, attribute_name, value)

Usage:
  python collate_biosample_html_long_wide.py <input_dir> \
      --wide-out biosample_wide.tsv \
      --long-out biosample_long.tsv
"""

import argparse
import csv
import re
from pathlib import Path
from html import unescape
import xml.etree.ElementTree as ET

# Core fields commonly present in BioSample esummary DocumentSummary
CORE_FIELDS = [
    "Id",             # NCBI integer ID
    "Accession",      # e.g., SAMN12345678
    "Title",
    "Organism",
    "TaxID",
    "BioProject",
    "SubmissionDate",
    "LastUpdate",
    "Owner",
    "Model",
]

# Attribute names that likely describe "sample source"
SOURCE_KEYS = [
    "isolation_source",
    "sample source",
    "sample_source",
    "sample_type",
    "isolation source",
    "source",
    "specimen_source",
    "host_body_site",
    "body_site",
    "body site",
]

YEAR_REGEXES = [
    re.compile(r"\b(19\d{2}|20\d{2}|210\d)\b"),  # four-digit year 1900–2109
    re.compile(r"(\d{4})-\d{2}-\d{2}"),          # ISO date
]

def findall_docsummaries(root):
    docs = root.findall(".//DocumentSummary")
    if docs:
        return docs
    if root.tag == "DocumentSummary":
        return [root]
    return []

def text_or_none(elem, tag):
    node = elem.find(tag)
    return node.text.strip() if (node is not None and node.text) else None

def attr_name(attr_elem):
    """Prefer harmonized_name, then attribute_name, otherwise 'Attribute'."""
    return (
        (attr_elem.get("harmonized_name") or "").strip()
        or (attr_elem.get("attribute_name") or "").strip()
        or "Attribute"
    )

def normalize_key(s):
    """Lowercase, collapse spaces/underscores to single underscore for matching."""
    return re.sub(r"[\s]+", "_", s.strip().lower())

def extract_collection_year(values):
    """
    Given a list of strings from collection_date-like fields, return the first year found.
    Tries several regexes. Returns '' if none.
    """
    for v in values:
        if not v:
            continue
        for rx in YEAR_REGEXES:
            m = rx.search(v)
            if m:
                # The first capturing group should be the year
                return m.group(1)
    return ""

def parse_biosample_file(path):
    """
    Parse one esummary XML file (even if saved with .html extension).
    Returns:
      - rows_wide: list[dict]  (usually length 1) → core fields + attributes
      - rows_long: list[dict]  (sample, attribute, value triplets)
    """
    raw = path.read_text(encoding="utf-8", errors="ignore")
    raw = unescape(raw)

    try:
        root = ET.fromstring(raw)
    except ET.ParseError as e:
        raise RuntimeError(f"XML parse error in {path.name}: {e}") from e

    rows_wide = []
    rows_long = []

    for ds in findall_docsummaries(root):
        row = {}

        # Core fields
        for f in CORE_FIELDS:
            row[f] = text_or_none(ds, f)

        # Fall back accession from filename if missing
        if not row.get("Accession"):
            row["Accession"] = path.stem

        # Collect attributes
        # Typical structure: <Attributes><Attribute attribute_name="host">Homo sapiens</Attribute>...</Attributes>
        attr_values_by_name = {}  # normalized name → list of values
        for attr in ds.findall(".//Attributes/Attribute"):
            key_raw = attr_name(attr)
            key_norm = normalize_key(key_raw)
            value = (attr.text or "").strip()
            if not value:
                continue
            attr_values_by_name.setdefault(key_norm, []).append(value)

            # Add to long format (use the raw name for readability)
            rows_long.append({
                "Accession": row["Accession"],
                "Attribute": key_raw,
                "Value": value
            })

        # Store attributes into the wide row (join multiple values with " | ")
        for key_norm, values in attr_values_by_name.items():
            # Keep a readable column name: normalized but not screaming caps
            col = key_norm
            row[col] = " | ".join(values)

        # Extract sample_source (best-effort from likely fields)
        sample_source = ""
        for k in SOURCE_KEYS:
            k_norm = normalize_key(k)
            if k_norm in attr_values_by_name and attr_values_by_name[k_norm]:
                sample_source = " | ".join(attr_values_by_name[k_norm])
                break
        row["sample_source"] = sample_source

        # Extract collection_year from any collection_date-like fields
        collection_candidates = []
        for k in list(attr_values_by_name.keys()):
            if "collection_date" in k or k == "collection_date":
                collection_candidates.extend(attr_values_by_name[k])
        row["collection_year"] = extract_collection_year(collection_candidates)

        # Also capture multiple BioProject accessions if present as separate elements
        projects = [el.text.strip() for el in ds.findall(".//BioProjectAccn") if el.text]
        if projects:
            row.setdefault("bioprojectaccn_list", " | ".join(projects))
            if not row.get("BioProject"):
                row["BioProject"] = projects[0]

        rows_wide.append(row)

    return rows_wide, rows_long

def write_tsv(path, rows, columns=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        # Write header if columns provided, else do nothing
        if columns:
            with path.open("w", newline="", encoding="utf-8") as out:
                w = csv.DictWriter(out, fieldnames=columns, delimiter="\t", extrasaction="ignore")
                w.writeheader()
        return

    if columns is None:
        # infer all keys observed
        keys = set()
        for r in rows:
            keys.update(r.keys())
        columns = sorted(keys)

    with path.open("w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

def main():
    p = argparse.ArgumentParser(description="Collate BioSample esummary XML/HTML files into wide and long TSVs, extracting sample_source and collection_year.")
    p.add_argument("input_dir", type=Path, help="Directory containing *.html (esummary) files")
    p.add_argument("--wide-out", type=Path, default=Path("biosample_wide.tsv"), help="Output TSV (wide format)")
    p.add_argument("--long-out", type=Path, default=Path("biosample_long.tsv"), help="Output TSV (long/tidy format)")
    args = p.parse_args()

    if not args.input_dir.is_dir():
        raise SystemExit(f"Input directory not found: {args.input_dir}")

    files = sorted(args.input_dir.glob("*.html"))
    if not files:
        raise SystemExit(f"No .html files found in {args.input_dir}")

    all_wide = []
    all_long = []

    for f in files:
        try:
            w, l = parse_biosample_file(f)
            all_wide.extend(w)
            all_long.extend(l)
        except Exception as e:
            print(f"[WARN] Skipping {f.name}: {e}")

    if not all_wide:
        raise SystemExit("No records parsed successfully.")

    # Prepare columns: core fields first, then standard extras, then the rest alphabetically
    core = list(CORE_FIELDS)
    # Ensure Accession is present even if not in XML core
    if "Accession" not in core:
        core.insert(0, "Accession")

    standard_extras = ["sample_source", "collection_year", "bioprojectaccn_list"]

    # Gather dynamic attribute columns from observed rows
    dynamic_cols = sorted({k for r in all_wide for k in r.keys()} - set(core) - set(standard_extras))

    wide_cols = core + standard_extras + dynamic_cols

    write_tsv(args.wide_out, all_wide, wide_cols)
    write_tsv(args.long_out, all_long, columns=["Accession", "Attribute", "Value"])

    print(f"Wrote wide table: {args.wide_out}")
    print(f"Wrote long table: {args.long_out}")

if __name__ == "__main__":
    main()