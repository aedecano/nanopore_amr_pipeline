#!/usr/bin/env bash
set -euo pipefail

# Usage: ./bam2fastq_batch.sh /path/to/bams /path/to/out
IN_DIR="${1:-.}"
OUT_DIR="${2:-fastq_out}"

mkdir -p "$OUT_DIR"

# Choose compressor: pigz (faster, multi-core) if available, else gzip
if command -v pigz >/dev/null 2>&1; then
  COMP="pigz -c"
  EXT="fastq.gz"
else
  COMP="gzip -c"
  EXT="fastq.gz"
fi

shopt -s nullglob
found=0
for bam in "$IN_DIR"/*.bam; do
  found=1
  base=$(basename "$bam" .bam)
  out="$OUT_DIR/${base}.${EXT}"
  log="$OUT_DIR/${base}.log"

  echo "[`date '+%F %T'`] Converting $bam -> $out"
  # For ONT data (single-end), one FASTQ per BAM:
  # -n keeps original read names; adjust/remove if you prefer renumbered names.
  { samtools fastq -n "$bam" 2> "$log" | $COMP > "$out"; } || {
    echo "ERROR converting $bam — see $log" >&2
    exit 1
  }
done

if [[ $found -eq 0 ]]; then
  echo "No BAM files found in: $IN_DIR" >&2
  exit 2
fi

echo "Done. FASTQs saved to: $OUT_DIR"