#!/usr/bin/env bash
set -euo pipefail

input_list="bakta_list.txt"
outdir="merged_gff"
mkdir -p "$outdir"

echo "[i] Merging GFF + FASTA pairs listed in: $input_list"
echo "[i] Output directory: $outdir"
count=0

# Read tab-separated lines: <gff> <fasta>
while IFS=$'\t' read -r gff fa; do
  # skip empty or commented lines
  [[ -z "$gff" || "$gff" =~ ^# ]] && continue

  if [[ ! -f "$gff" ]]; then
    echo "[WARN] Missing GFF: $gff" >&2
    continue
  fi
  if [[ ! -f "$fa" ]]; then
    echo "[WARN] Missing FASTA: $fa" >&2
    continue
  fi

  base="$(basename "$gff" .gff3)"
  out="${outdir}/${base}.merged.gff"

  echo "[+] Processing: $base"

  # Write GFF section only (strip any previous FASTA)
  awk 'BEGIN{inF=0} /^##FASTA/{inF=1; next} !inF{print}' "$gff" > "$out"

  # Append new FASTA
  echo '##FASTA' >> "$out"
  case "$fa" in
    *.gz) gunzip -c "$fa" >> "$out" ;;
    *)    cat "$fa" >> "$out" ;;
  esac

  ((count++))
done < "$input_list"

echo "[✓] Completed merging $count files into: $outdir"
ls -lh "$outdir" | head

# make a new list for Panaroo
ls merged_gff/*.gff > merged_list.txt