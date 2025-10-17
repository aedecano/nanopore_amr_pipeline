#!/bin/bash
# Usage: ./download_from_ftp.sh input_file.txt [BASE_URL]
# input file: <sample_id><whitespace><ftp_or_path>
# If links start with /ebi/ftp/... we'll map them to https://ftp.ebi.ac.uk/pub/...

set -euo pipefail

input_file="${1:-}"
BASE_URL_DEFAULT="https://ftp.ebi.ac.uk"
BASE_URL="${2:-$BASE_URL_DEFAULT}"

if [[ -z "$input_file" || ! -f "$input_file" ]]; then
  echo "Usage: $0 input_file [BASE_URL]"
  echo "Example: $0 test_ftpLinks.txt https://ftp.ebi.ac.uk"
  exit 1
fi

# Trim CRLF helper (in case the file came from Windows)
trim_cr() {
  awk '{ sub(/\r$/,""); print }'
}

while IFS=$' \t' read -r sample raw_link || [[ -n "${sample:-}" ]]; do
  # Skip blank/comment lines
  [[ -z "${sample:-}" ]] && continue
  [[ "${sample:0:1}" == "#" ]] && continue

  # Clean CR if present
  raw_link="$(printf '%s' "${raw_link:-}" | tr -d '\r')"

  # Normalize the link
  link="$raw_link"
  case "$link" in
    http://*|https://*|ftp://*)
      # already a full URL
      ;;
    /ebi/ftp/*)
      # Map /ebi/ftp/... to https://ftp.ebi.ac.uk/...
      # On the host, the path is typically /pub/...
      # Strip the /ebi/ftp prefix so /ebi/ftp/pub/... -> /pub/...
      path="${link#/ebi/ftp}"
      link="${BASE_URL%/}${path}"
      ;;
    /*)
      # Generic absolute path -> append to BASE_URL
      link="${BASE_URL%/}${link}"
      ;;
    *)
      # Relative path -> BASE_URL + / + path
      link="${BASE_URL%/}/${link}"
      ;;
  esac

  echo "Downloading for sample ${sample}: ${link}"
  # Resume enabled; save with original filename (or customize -O if you want sample-prefixed files)
  wget -c "${link}"
done < <(trim_cr < "$input_file")