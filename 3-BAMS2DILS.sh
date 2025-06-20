#!/bin/bash

# --- Parse named arguments ---
while [[ "$#" -gt 0 ]]; do
  case $1 in
    -genome) GENOME="$2"; shift ;;
    -loci) BED="$2"; shift ;;
    -popmap) POP_MAP="$2"; shift ;;
    -out) OUTDIR="$2"; shift ;;
    -minDP) MIN_DP="$2"; shift ;;
    -mapQ) MAPQ="$2"; shift ;;
    -baseQ) BASEQ="$2"; shift ;;
    -h|--help)
      echo "Usage: $0 -genome <genome.fa> -loci <loci.bed> -popmap <pop_map.txt> -out <output_dir> [-minDP 5] [-mapQ 20] [-baseQ 20]"
      exit 0 ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
  shift
done

# Set defaults if not provided
MIN_DP="${MIN_DP:-5}"
MAPQ="${MAPQ:-20}"
BASEQ="${BASEQ:-20}"

# Basic check
if [[ -z "$GENOME" || -z "$BED" || -z "$POP_MAP" || -z "$OUTDIR" ]]; then
  echo "❌ Missing required argument. Use -h for help."
  exit 1
fi

# Create output directory
mkdir -p "$OUTDIR"

# Print config
echo "Running with:"
echo "  GENOME:  $GENOME"
echo "  BED:     $BED"
echo "  POP_MAP: $POP_MAP"
echo "  OUTDIR:  $OUTDIR"
echo "  MIN_DP:  $MIN_DP"
echo "  MAPQ:    $MAPQ"
echo "  BASEQ:   $BASEQ"

mkdir -p "$OUTDIR"


# ========== LOOP OVER LOCI ==========

i=0
while IFS=$'\t' read -r CHROM START END; do
    ((i++))
    LOCUS=$(printf "Locus%04d" "$i")
    REGION="${CHROM}:$((START + 1))-${END}"
    FASTA_OUT="$OUTDIR/${LOCUS}.fasta"

    echo "?? Processing $REGION as $LOCUS"

    # Extract reference sequence for the locus (shared template)
    REF_TMP="ref_${LOCUS}.tmp"
    samtools faidx "$GENOME" "$REGION" | grep -v '^>' | tr -d '\n' > "$REF_TMP"

    # LOOP OVER INDIVIDUALS
    while IFS=$'\t' read -r POP SAMPLE BAM; do
        echo "  ?? $SAMPLE ($POP)"

        # Define unique temporary files per sample-locus
        DEPTH_TMP="depth_${SAMPLE}_${LOCUS}.tmp"
        VARIANT_TSV="variants_${SAMPLE}_${LOCUS}.tsv"
        VCF_GZ="vcf_${SAMPLE}_${LOCUS}.vcf.gz"

        # Generate depth per position
        samtools mpileup -aa -f "$GENOME" -r "$REGION" -q "$MAPQ" -Q "$BASEQ" "$BAM" | \
            awk '!/^#/ { print $4 }' > "$DEPTH_TMP"

        # Call variants from BAM
        bcftools mpileup -f "$GENOME" -r "$REGION" -a FORMAT/DP  -Q "$BASEQ" -q "$MAPQ" "$BAM" -Ou | \
        bcftools call -mv -Ou -f GQ -A | \
        bcftools view -v snps -m2 -M2 -Ou | \
        bcftools filter -e 'strlen(REF)!=1 || strlen(ALT)!=1' -Ou | \
        bcftools norm -f "$GENOME" -Oz -o "$VCF_GZ"
        bcftools index "$VCF_GZ"

        # Extract position-wise genotype and DP
        bcftools query -f '%POS\t%REF\t%ALT\t[%GT]\t[%DP]\n' -r "$REGION" "$VCF_GZ" > "$VARIANT_TSV"

        # === BUILD HAPLOTYPES ===
        awk -v start="$START" -v end="$END" -v min_dp="$MIN_DP" -v POP="$POP" -v SAMPLE="$SAMPLE" -v LOCUS="$LOCUS" \
            -v ref_file="$REF_TMP" -v depth_file="$DEPTH_TMP" '
        BEGIN {
            # Load reference
            getline ref_seq < ref_file
            ref_length = length(ref_seq)
            for (i = 1; i <= ref_length; i++) {
                hap1[i] = substr(ref_seq, i, 1)
                hap2[i] = substr(ref_seq, i, 1)
            }

            # Load per-base depth
            i = 1
            while ((getline d < depth_file) > 0) {
                base_dp[i] = d
                i++
            }
        }
        {
            pos = $1 - (start + 1) + 1
            ref = $2
            alt = $3
            gt  = $4
            dp  = $5

            if (pos < 1 || pos > ref_length) next

            if (dp == "." || dp < min_dp || gt == "./.") {
                hap1[pos] = "N"
                hap2[pos] = "N"
            } else if (gt == "0/0") {
                hap1[pos] = ref
                hap2[pos] = ref
            } else if (gt == "1/1") {
                hap1[pos] = alt
                hap2[pos] = alt
            } else if (gt ~ /0\/1|1\/0/) {
                if (int(rand()*2) == 0) {
                    hap1[pos] = ref
                    hap2[pos] = alt
                } else {
                    hap1[pos] = alt
                    hap2[pos] = ref
                }
            }
        }
        END {
            for (i = 1; i <= ref_length; i++) {
                if (base_dp[i] < min_dp) {
                    hap1[i] = "N"
                    hap2[i] = "N"
                }
            }

            printf ">%s|%s|%s|allele_1\n", LOCUS, POP, SAMPLE
            for (i = 1; i <= ref_length; i++) printf "%s", hap1[i]
            print ""

            printf ">%s|%s|%s|allele_2\n", LOCUS, POP, SAMPLE
            for (i = 1; i <= ref_length; i++) printf "%s", hap2[i]
            print ""
        }
        ' "$VARIANT_TSV" >> "$FASTA_OUT"

        # Clean up
        rm -f "$VCF_GZ" "$VCF_GZ.csi" "$VARIANT_TSV" "$DEPTH_TMP"
    done < "$POP_MAP"

    rm -f "$REF_TMP"
done < "$BED"

