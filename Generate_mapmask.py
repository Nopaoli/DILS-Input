import argparse
from Bio import SeqIO

def fasta_to_bed_intervals(fasta_file, output_bed):
    """
    Extract contiguous regions with mappability value `3` from a multi-chromosome FASTA file and write to a BED file.

    Parameters:
    - fasta_file: Path to the input FASTA file.
    - output_bed: Path to the output BED file.
    """
    with open(output_bed, "w") as bed:
        for record in SeqIO.parse(fasta_file, "fasta"):
            chrom_name = record.id  # Use the FASTA record ID as the chromosome name
            sequence = record.seq
            start = None
            for idx, value in enumerate(sequence):
                if value == "3":  # Check for `c3` values
                    if start is None:
                        start = idx  # Mark the start of a `c3` region
                elif start is not None:
                    # Write the interval when a `c3` region ends
                    bed.write(f"{chrom_name}\t{start}\t{idx}\n")
                    start = None
            # Handle the case where a `c3` region extends to the end of the sequence
            if start is not None:
                bed.write(f"{chrom_name}\t{start}\t{len(sequence)}\n")

def main():
    parser = argparse.ArgumentParser(description="Convert a mappability FASTA file to a BED file with uniquely mappable (c3) regions.")
    parser.add_argument("fasta_file", help="Input mappability FASTA file")
    parser.add_argument("output_bed", help="Output BED file with uniquely mappable regions")
    args = parser.parse_args()

    fasta_to_bed_intervals(args.fasta_file, args.output_bed)
    print(f"BED file with uniquely mappable regions written to: {args.output_bed}")

if __name__ == "__main__":
    main()
