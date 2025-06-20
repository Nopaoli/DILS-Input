#!/bin/bash

# ===============================
# Mappability Mask Generator
# ===============================
# Arguments:
#   -genome    <path to genome FASTA>
#   -dir       <output directory>
#   -readlen   <read length (e.g. 150)>
#   -seqbility <path to seqbility tools directory>
# ===============================

# ========== Parse Arguments ==========
while [[ "$#" -gt 0 ]]; do
  case $1 in
    -genome) GENOME="$(readlink -f "$2")"; shift ;;
    -dir) WORKDIR="$2"; shift ;;
    -readlen) READLEN="$2"; shift ;;
    -seqbility) SEQBILITY_DIR="$2"; shift ;;
    -h|--help)
      echo "Usage: $0 -genome genome.fa -dir output_dir -readlen 150 -seqbility /path/to/seqbility"
      exit 0 ;;
    *) echo "❌ Unknown parameter: $1"; exit 1 ;;
  esac
  shift
done

# ========== Validate Inputs ==========
if [[ -z "$GENOME" || -z "$WORKDIR" || -z "$READLEN" || -z "$SEQBILITY_DIR" ]]; then
  echo "❌ Missing required arguments."
  echo "Usage: $0 -genome genome.fa -dir output_dir -readlen 150 -seqbility /path/to/seqbility"
  exit 1
fi

# ========== Set Tool Paths ==========
SPLITFA="$(readlink -f "$SEQBILITY_DIR/splitfa")"
GEN_RAW_MASK="$(readlink -f "$SEQBILITY_DIR/gen_raw_mask.pl")"
GEN_MASK="$(readlink -f "$SEQBILITY_DIR/gen_mask")"
# Resolve absolute path to this script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# ========== Prepare Output Directory ==========
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

echo "📂 Working in: $WORKDIR"
echo "🧬 Using genome: $GENOME"
echo "📏 Read length: $READLEN"
echo "🔧 Seqbility tools: $SEQBILITY_DIR"

export PATH="$SEQBILITY_DIR:$PATH"

# ========== Step 1: Split genome ==========
echo "🔧 Splitting genome into $READLEN-bp fragments..."
$SPLITFA "$GENOME" "$READLEN" | split -l 20000000

# ========== Step 2: Align reads with BWA ==========
echo "🚀 Aligning fragments with BWA..."
for i in x??; do
  echo "  Aligning $i"
  bwa aln -R 1000000 -t 4 -O 3 -E 3 "$GENOME" "$i" > "${i}.sai"
done

# ========== Step 3: Convert to SAM ==========
echo "📄 Generating SAM files..."
for i in x??; do
  bwa samse "$GENOME" "${i}.sai" "$i" > "${i}.sam"
done

# ========== Step 4: Generate Mappability Mask ==========
echo "🧬 Generating mappability mask..."
cat x??.sam | "$GEN_RAW_MASK" > "rawMask_${READLEN}.fa"
"$GEN_MASK" -l "$READLEN" -r 0.5 "rawMask_${READLEN}.fa" > "mask_${READLEN}_50.fa"

# ========== Step 5: Convert Mask to BED ==========
echo "📁 Converting mask to BED..."
python "$SCRIPT_DIR/Generate_mapmask.py" "mask_${READLEN}_50.fa" mappable_regions.bed
# ========== Step 6: Cleanup ==========
echo "🧹 Cleaning up intermediate files..."
rm -f x?? x??.sai x??.sam rawMask_${READLEN}.fa

echo "✅ Done: mappable_regions.bed created in $WORKDIR"

