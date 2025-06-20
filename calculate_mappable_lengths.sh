#!/bin/bash

BEDFILE="$1"

if [[ ! -f "$BEDFILE" ]]; then
  echo "❌ Error: File '$BEDFILE' not found."
  echo "Usage: $0 path/to/mappable_regions.bed"
  exit 1
fi

echo "📄 Calculating lengths from: $BEDFILE"

# Corrected awk line
awk '{ print $0 "\t" ($3 - $2) }' "$BEDFILE" > lengths.bed
# Summary stats
echo "📊 Summary:"
awk '{ sum += ($3 - $2); count++ } END { print "  Regions: " count "\n  Total length: " sum "\n  Average length: " (sum/count) }' "$BEDFILE"

