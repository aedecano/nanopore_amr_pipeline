while read -r gff fa; do
  sid=$(basename "$gff")
  g_ids=$(grep -v '^#' "$gff" | awk -F'\t' 'NF>1{print $1}' | sort -u)
  f_ids=$(grep '^>' "$fa" | sed 's/^>//; s/ .*//' | sort -u)
  missing=$(comm -23 <(printf "%s\n" $g_ids | sort -u) <(printf "%s\n" $f_ids | sort -u) | head -n1)
  if [ -n "$missing" ]; then
    echo "MISMATCH: $sid -> seqid not in FASTA: $missing"
  else
    echo "OK: $sid"
  fi
done < bakta_list.txt