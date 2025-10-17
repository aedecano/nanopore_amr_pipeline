#!/usr/bin/env python3
"""
Fetch BioSample ESummary XML (saved with .html extension) for accessions listed in a TSV.

Usage:
  python3 fetch_biosample_htmls_from_tsv.py <input.tsv> [output_dir]

- Auto-detects a BioSample accession column (common names tried).
- Accepts accessions like SAMN*, SAMD*, SAME* or numeric UIDs.
- For accessions, resolves via: esearch -db biosample -query <acc> | esummary
- Skips files that already exist.
- Requires Entrez Direct in PATH: esearch, esummary.
- Set env NCBI_API_KEY to improve rate limits (optional).

Optional (commented-out) section provided to fetch assemblies from ENA.
"""

import csv
import os
import sys
import time
import shutil
import subprocess
from pathlib import Path
from typing import Optional, List

CANDIDATE_COLS = [
    "biosample", "BioSample", "biosample_accession", "BioSampleAccession",
    "sample_accession", "sample", "SAMPLE", "biosample_id", "BioSample ID"
]

def die(msg: str, code: int = 1) -> None:
    print(msg, file=sys.stderr); sys.exit(code)

def have_tool(name: str) -> bool:
    return shutil.which(name) is not None

def is_numeric_uid(s: str) -> bool:
    return s.isdigit()

def run_cmd(args, *, stdin_data: Optional[str] = None, timeout: int = 120) -> subprocess.CompletedProcess:
    return subprocess.run(
        args, input=stdin_data, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, check=True, timeout=timeout
    )

def fetch_biosample_xml(token: str, retries: int = 3, backoff: float = 0.6) -> Optional[str]:
    """Return ESummary XML text for a BioSample token (UID or accession)."""
    for attempt in range(1, retries + 1):
        try:
            if is_numeric_uid(token):
                p = run_cmd(["esummary", "-db", "biosample", "-id", token])
                return p.stdout
            p1 = run_cmd(["esearch", "-db", "biosample", "-query", token])
            p2 = run_cmd(["esummary"], stdin_data=p1.stdout)
            return p2.stdout
        except subprocess.CalledProcessError as e:
            if attempt < retries:
                time.sleep(backoff * attempt)
                continue
            sys.stderr.write(f"[WARN] fetch failed for {token}: {e.stderr or e}\n")
            return None
        except Exception as e:
            if attempt < retries:
                time.sleep(backoff * attempt)
                continue
            sys.stderr.write(f"[WARN] fetch error for {token}: {e}\n")
            return None
    return None

def detect_accession_column(header: List[str]) -> Optional[str]:
    lower = {h.lower(): h for h in header}
    for cand in CANDIDATE_COLS:
        if cand.lower() in lower:
            return lower[cand.lower()]
    for h in header:
        if h.lower().startswith("bio") or h.lower().endswith("accession"):
            return h
    return None

def split_multi(val: str) -> List[str]:
    parts = []
    for sep in (";", ","):
        if sep in val:
            parts = [p.strip() for p in val.split(sep)]
            break
    return parts or [val.strip()]

def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        die(f"Usage: {Path(sys.argv[0]).name} <input.tsv> [output_dir]", 2)

    tsv_path = Path(sys.argv[1])
    out_dir = Path(sys.argv[2]) if len(sys.argv) == 3 else Path("BioSample_HTMLs")

    if not tsv_path.is_file():
        die(f"Error: TSV not found: {tsv_path}")

    if not have_tool("esearch") or not have_tool("esummary"):
        die("Error: Requires Entrez Direct (esearch, esummary) in PATH.")

    out_dir.mkdir(parents=True, exist_ok=True)

    api_key = os.environ.get("NCBI_API_KEY")
    if not api_key:
        print("[INFO] NCBI_API_KEY not set; consider setting it for better throughput.", file=sys.stderr)

    total, written, skipped, failed = 0, 0, 0, 0

    with tsv_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            die("Error: TSV has no header row.")
        acc_col = detect_accession_column(reader.fieldnames)
        if not acc_col:
            die(f"Error: Could not find a BioSample accession column. "
                f"Tried: {', '.join(CANDIDATE_COLS)}. Columns present: {', '.join(reader.fieldnames)}")

        print(f"[INFO] Using accession column: '{acc_col}'", file=sys.stderr)

        for row in reader:
            raw = (row.get(acc_col) or "").strip()
            if not raw:
                continue
            for token in split_multi(raw):
                if not token:
                    continue
                total += 1
                out_path = out_dir / f"{token}.html"
                if out_path.exists():
                    skipped += 1
                    print(f"[SKIP] exists: {out_path}", file=sys.stderr)
                    continue

                xml = fetch_biosample_xml(token)
                if not xml:
                    failed += 1
                    try:
                        out_path.write_text(f"<!-- fetch failed for {token} -->\n", encoding="utf-8")
                    except Exception as e:
                        print(f"[WARN] Stub write failed for {token}: {e}", file=sys.stderr)
                    continue

                # Save even if there's an <ERROR> in payload (for inspection)
                try:
                    out_path.write_text(xml, encoding="utf-8")
                    written += 1
                    print(f"[OK] {token} -> {out_path}", file=sys.stderr)
                except Exception as e:
                    failed += 1
                    print(f"[WARN] Write failed for {token}: {e}", file=sys.stderr)

                # --------------------------------------------------------
                # ⬇️ OPTIONAL: download ENA assembly (commented out)
                # --------------------------------------------------------
                # If you later want to download assemblies directly from ENA,
                # uncomment and adjust the following block.
                #
                # Example assumes you have a column named 'run_accession' or 'assembly_accession'
                # in your TSV. You can modify as needed (e.g., using FTP links).
                #
                # ena_acc = (row.get("run_accession") or row.get("assembly_accession") or "").strip()
                # if ena_acc:
                #     ena_outdir = Path("ENA_assemblies")
                #     ena_outdir.mkdir(exist_ok=True)
                #     try:
                #         ena_url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{ena_acc}?download=true"
                #         out_fasta = ena_outdir / f"{ena_acc}.fasta.gz"
                #         if not out_fasta.exists():
                #             print(f"[ENA] Downloading {ena_acc} -> {out_fasta}", file=sys.stderr)
                #             subprocess.run(
                #                 ["curl", "-L", "-o", str(out_fasta), ena_url],
                #                 check=True, text=True
                #             )
                #         else:
                #             print(f"[ENA] Skipped existing {out_fasta}", file=sys.stderr)
                #     except subprocess.CalledProcessError as e:
                #         print(f"[WARN] ENA download failed for {ena_acc}: {e}", file=sys.stderr)
                # --------------------------------------------------------

    print(f"\nDone. total={total} written={written} skipped={skipped} failed={failed}", file=sys.stderr)

if __name__ == "__main__":
    main()