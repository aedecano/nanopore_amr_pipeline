#!/usr/bin/env python3

#Usage: python3 collate_parse_biosampleIDs_metadata.py Blackwell2021_Biosample_HTMLs/   -o biosample.tsv --out-csv biosample.csv --verbose

import argparse, csv, re, time, sys
from pathlib import Path
from html import unescape
import xml.etree.ElementTree as ET
from urllib.parse import urlencode
from urllib.request import urlopen, Request

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
USER_AGENT = "biosample_collate/1.0 (+https://example.org)"
DEFAULT_EMAIL = None  # optionally fill or pass --email

CORE_FIELDS = [
    "Id","Accession","Title","Organism","TaxID","BioProject",
    "SubmissionDate","LastUpdate","Owner","Model"
]

SOURCE_KEYS = {"isolation_source","sample_source","source","sample_type","isolation"}
DATE_KEYS = {"collection_date","collection year","date","year"}
YEAR_RE = re.compile(r"\b(19|20)\d{2}\b")

def vprint(enabled, *a, **k):
    if enabled: print(*a, **k, file=sys.stderr, flush=True)

def fetch_url(path, params, email=None, tool="biosample_collate", verbose=False):
    q = dict(params)
    if email: q["email"] = email
    if tool: q["tool"] = tool
    url = EUTILS_BASE + path + "?" + urlencode(q)
    req = Request(url, headers={"User-Agent": USER_AGENT})
    vprint(verbose, f"[HTTP] GET {url}")
    with urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "ignore")

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

def get_attribute_name(attr_elem):
    return (attr_elem.get("harmonized_name")
            or attr_elem.get("attribute_name")
            or "Attribute")

def parse_biosample_xml(raw_text, filename_hint=None):
    raw = unescape(raw_text)
    try:
        root = ET.fromstring(raw)
    except ET.ParseError as e:
        raise RuntimeError(f"XML parse error in {filename_hint or '<memory>'}: {e}")
    rows = []
    for ds in findall_docsummaries(root):
        row = {}
        for f in CORE_FIELDS:
            row[f] = text_or_none(ds, f)
        for attr in ds.findall(".//Attributes/Attribute"):
            name = get_attribute_name(attr)
            value = (attr.text or "").strip()
            if name in row and row[name]:
                if value: row[name] = f"{row[name]} | {value}"
            else:
                row[name] = value
        projects = [el.text.strip() for el in ds.findall(".//BioProjectAccn") if el.text]
        if projects:
            row["BioProjectAccn_list"] = " | ".join(projects)
            if not row.get("BioProject"):
                row["BioProject"] = projects[0]
        # derived
        sample_source, collection_year = None, None
        for k, v in list(row.items()):
            if not v: continue
            kl = k.lower()
            if sample_source is None and any(key in kl for key in SOURCE_KEYS):
                sample_source = v
            if collection_year is None and any(key in kl for key in DATE_KEYS):
                m = YEAR_RE.search(v)
                if m: collection_year = m.group(0)
        row["SampleSource"] = sample_source
        row["CollectionYear"] = collection_year
        rows.append(row)
    return rows

def parse_biosample_file(path: Path):
    return parse_biosample_xml(path.read_text(encoding="utf-8", errors="ignore"), filename_hint=path.name)

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

def fetch_biosample_summaries(ids, email, batch, delay, verbose):
    all_rows = []
    for i, chunk_ids in enumerate(chunks(ids, batch), start=1):
        xml = fetch_url("esummary.fcgi", {
            "db":"biosample",
            "id":",".join(chunk_ids),
            "retmode":"xml"
        }, email=email, verbose=verbose)
        rows = parse_biosample_xml(xml, filename_hint=f"batch_{i}")
        # Try to fill missing accession from Title
        for r in rows:
            if not r.get("Accession"):
                title = (r.get("Title") or "")
                m = re.search(r"(SAM[NES]\w+\d+)", title)
                if m: r["Accession"] = m.group(1)
        all_rows.extend(rows)
        vprint(verbose, f"[INFO] fetched batch {i} ({len(chunk_ids)} IDs); rows so far: {len(all_rows)}")
        time.sleep(delay)
    return all_rows

def elink_biosample_to_assembly(ids, email, batch, delay, verbose):
    mapping = {acc: [] for acc in ids}
    for i, chunk_ids in enumerate(chunks(ids, batch), start=1):
        xml = fetch_url("elink.fcgi", {
            "dbfrom":"biosample",
            "db":"assembly",
            "id":",".join(chunk_ids),
            "retmode":"xml"
        }, email=email, verbose=verbose)
        try:
            root = ET.fromstring(xml)
        except ET.ParseError:
            time.sleep(delay); continue
        # fallback: query one-by-one to ensure mapping
        for acc in chunk_ids:
            one = fetch_url("elink.fcgi", {
                "dbfrom":"biosample","db":"assembly","id":acc,"retmode":"xml"
            }, email=email, verbose=verbose)
            try:
                r = ET.fromstring(one)
                uids = [n.text.strip() for n in r.findall(".//LinkSetDb[DbTo='assembly']/Link/Id") if n.text]
                mapping[acc] = uids
            except ET.ParseError:
                mapping[acc] = []
            time.sleep(delay)
        vprint(verbose, f"[INFO] elink batch {i}: mapped {len(chunk_ids)} BioSamples")
        time.sleep(delay)
    return mapping

def esummary_assembly_ftps(uids, email, batch, delay, verbose):
    out = {}
    for i, chunk_uids in enumerate(chunks(uids, batch), start=1):
        xml = fetch_url("esummary.fcgi", {
            "db":"assembly","id":",".join(chunk_uids),"retmode":"xml"
        }, email=email, verbose=verbose)
        try:
            root = ET.fromstring(xml)
        except ET.ParseError:
            time.sleep(delay); continue
        for ds in root.findall(".//DocumentSummary"):
            uid = ds.get("uid") or text_or_none(ds, "Uid")
            if not uid: continue
            refseq = text_or_none(ds, "FtpPath_RefSeq")
            genbank = text_or_none(ds, "FtpPath_GenBank")
            out[uid] = {"FtpPath_RefSeq": refseq, "FtpPath_GenBank": genbank}
        vprint(verbose, f"[INFO] assembly esummary batch {i}: +{len(chunk_uids)}")
        time.sleep(delay)
    return out

def fasta_from_ftp(ftp_dir):
    if not ftp_dir: return None
    tail = ftp_dir.rstrip("/").split("/")[-1]
    return f"{ftp_dir.rstrip('/')}/{tail}_genomic.fna.gz"

def write_tables(rows, out_tsv: Path, out_csv: Path|None, verbose=False):
    if not rows:
        raise SystemExit("No records parsed successfully.")
    cols = list(CORE_FIELDS)
    preferred = [
        "BioProjectAccn_list","SampleSource","CollectionYear",
        "AssemblyFTP_RefSeq","AssemblyFTP_GenBank",
        "AssemblyFTP_FASTA_RefSeq","AssemblyFTP_FASTA_GenBank"
    ]
    extra = sorted({k for r in rows for k in r.keys()} - set(cols))
    ordered_pref = [k for k in preferred if k in extra]
    rest = [k for k in extra if k not in set(preferred)]
    cols.extend(ordered_pref + rest)

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows: w.writerow(r)
    vprint(verbose, f"[OK] Wrote TSV: {out_tsv}")

    if out_csv:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        with out_csv.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for r in rows: w.writerow(r)
        vprint(verbose, f"[OK] Wrote CSV: {out_csv}")

def main():
    p = argparse.ArgumentParser(description="Collate NCBI BioSample summaries from a folder or IDs file (pure Python).")
    p.add_argument("input_path", type=Path, help="Directory with .html/.xml OR text file of BioSample accessions.")
    p.add_argument("-o","--out-tsv", type=Path, default=Path("biosample_collated.tsv"))
    p.add_argument("--out-csv", type=Path, default=None)
    p.add_argument("--fetch-ftp", action="store_true", help="Also fetch Assembly FTP paths (slower).")
    p.add_argument("--batch-size", type=int, default=100, help="Batch size for network calls (default 100).")
    p.add_argument("--delay", type=float, default=0.35, help="Delay between calls (seconds).")
    p.add_argument("--email", type=str, default=DEFAULT_EMAIL, help="Contact email (E-utilities etiquette).")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    verbose = args.verbose
    vprint(verbose, f"[INFO] Mode check for {args.input_path}")

    all_rows = []

    if args.input_path.is_dir():
        files = sorted(list(args.input_path.glob("*.html")) + list(args.input_path.glob("*.xml")))
        if not files:
            raise SystemExit(f"No .html/.xml files found in {args.input_path}")
        vprint(verbose, f"[INFO] Found {len(files)} files")
        for f in files:
            try:
                rows = parse_biosample_file(f)
                for r in rows:
                    if not r.get("Accession"): r["Accession"] = f.stem
                all_rows.extend(rows)
            except Exception as e:
                print(f"[WARN] Skipping {f.name}: {e}", file=sys.stderr)
    elif args.input_path.is_file():
        ids = [ln.strip() for ln in args.input_path.read_text().splitlines() if ln.strip()]
        if not ids:
            raise SystemExit(f"No accessions found in {args.input_path}")
        # de-dup preserve order
        seen=set(); ids=[x for x in ids if not (x in seen or seen.add(x))]
        vprint(verbose, f"[INFO] Will fetch {len(ids)} BioSamples via E-utilities")
        all_rows = fetch_biosample_summaries(ids, email=args.email, batch=args.batch_size,
                                             delay=args.delay, verbose=verbose)
        if not all_rows:
            raise SystemExit("No records parsed successfully from IDs.")
    else:
        raise SystemExit(f"Input path not found: {args.input_path}")

    if args.fetch_ftp and all_rows:
        accessions = sorted({r.get("Accession") for r in all_rows if r.get("Accession", "").startswith("SAM")})
        if accessions:
            vprint(verbose, f"[INFO] Linking {len(accessions)} BioSamples → Assembly")
            mapping = elink_biosample_to_assembly(accessions, email=args.email,
                                                  batch=min(args.batch_size, 50),
                                                  delay=args.delay, verbose=verbose)
            all_uids = sorted({u for lst in mapping.values() for u in lst})
            vprint(verbose, f"[INFO] Fetching FTP for {len(all_uids)} assemblies")
            uid2ftp = esummary_assembly_ftps(all_uids, email=args.email,
                                             batch=min(args.batch_size, 200),
                                             delay=args.delay, verbose=verbose)
            for r in all_rows:
                acc = r.get("Accession")
                if not acc: continue
                uids = mapping.get(acc, [])
                refseq = genbank = None
                for u in uids:
                    info = uid2ftp.get(u, {})
                    refseq = refseq or info.get("FtpPath_RefSeq")
                    genbank = genbank or info.get("FtpPath_GenBank")
                    if refseq and genbank: break
                r["AssemblyFTP_RefSeq"] = refseq
                r["AssemblyFTP_GenBank"] = genbank
                r["AssemblyFTP_FASTA_RefSeq"] = fasta_from_ftp(refseq) if refseq else None
                r["AssemblyFTP_FASTA_GenBank"] = fasta_from_ftp(genbank) if genbank else None

    write_tables(all_rows, args.out_tsv, args.out_csv, verbose=verbose)
    # Also echo to STDOUT so you always see *something*
    print(f"Wrote TSV: {args.out_tsv}")
    if args.out_csv: print(f"Wrote CSV: {args.out_csv}")

if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        # surface the reason prominently
        print(str(e), file=sys.stderr)
        raise