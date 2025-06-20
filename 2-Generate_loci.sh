#!/bin/bash

# Usage check
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <mappable_regions.bed> <window_size> <num_loci>"
  exit 1
fi

MAPPABLE="$1"
WINDOW_SIZE="$2"
NUM_LOCI="$3"

SORTED_MAPPABLE="sorted_mappable.bed"
OUTPUT_TEMP="temp_selected_loci.bed"

# 1. Sort mappable BED
sort -k1,1 -k2,2n "$MAPPABLE" > "$SORTED_MAPPABLE"

# 2. Get total mappable length
TOTAL_MAPPABLE=$(awk '{sum += $3 - $2} END {print sum}' "$SORTED_MAPPABLE")

# 3. Calculate step size (bp between loci)
STEP=$(awk -v total="$TOTAL_MAPPABLE" -v n="$NUM_LOCI" 'BEGIN { printf "%.0f", total / n }')

echo "Total mappable length: $TOTAL_MAPPABLE bp"
echo "Target loci: $NUM_LOCI"
echo "Step size: $STEP bp"

# 4. Select windows at step intervals through mappable space
awk -v step="$STEP" -v win="$WINDOW_SIZE" -v max="$NUM_LOCI" '
BEGIN { OFS="\t"; next_target = 0; count = 0; genome_pos = 0 }
{
  chr = $1; start = $2; end = $3
  region_len = end - start

  while (next_target >= genome_pos && next_target < genome_pos + region_len) {
    offset = next_target - genome_pos
    win_start = start + offset
    win_end = win_start + win
    if (win_end <= end) {
      print chr, win_start, win_end
      count++
    }
    next_target += step
    if (count >= max) exit
  }
  genome_pos += region_len
}' "$SORTED_MAPPABLE" > "$OUTPUT_TEMP"

# 5. Count actual loci and rename output
ACTUAL=$(wc -l < "$OUTPUT_TEMP")
STEP_KB=$((STEP / 1000))
WIN_KB=$((WINDOW_SIZE / 1000))
OUTPUT_NAME="LOCI_N${ACTUAL}_len${WIN_KB}kb_mindist${STEP_KB}kb_mappable.bed"

mv "$OUTPUT_TEMP" "$OUTPUT_NAME"
echo "Output written to $OUTPUT_NAME"

