#!/bin/bash

DIR="${1:-.}"  # Default to current directory if not specified

echo -e "File\tNum_Seqs\tLength(s)\tUniform"

for fasta in "$DIR"/*.fasta; do
  if [[ ! -f "$fasta" ]]; then
    continue
  fi

  lengths=$(awk '
    /^>/ { if (len > 0) print len; len = 0; next }
    { len += length($0) }
    END { if (len > 0) print len }
  ' "$fasta")

  count=$(echo "$lengths" | wc -l)
  uniq_lengths=$(echo "$lengths" | sort -n | uniq | paste -sd "," -)
  n_uniq=$(echo "$uniq_lengths" | awk -F"," '{print NF}')
  uniform="YES"
  [[ "$n_uniq" -gt 1 ]] && uniform="NO"

  fname=$(basename "$fasta")
  echo -e "${fname}\t${count}\t${uniq_lengths}\t${uniform}"
done

