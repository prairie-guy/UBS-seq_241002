#!/usr/bin/env python3
import pysam, argparse
from xopen import xopen
from Bio.Seq import Seq
from seqpy import revcomp

def find_motif(fasta_file, chrom, pos, strand, motif_length):
    try:
        if strand == "+":
            seq = fasta_file.fetch(chrom, pos-1, pos+motif_length-1)
        else:
            seq = revcomp(fasta_file.fetch(chrom, pos-motif_length, pos))
        return(seq)
    except (ValueError, KeyError):
        return "N" * motif_length

def extract_motifs(input_file, output_file, motif_length, reference_fasta):
    with xopen(input_file) as input_file, xopen(output_file, "w") as output_file, pysam.FastaFile(reference_fasta) as fasta_file:
        header = next(input_file).strip()
        print(f"{header}\tMotif", file=output_file)
        for line in input_file:
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue
            sample, chrom, pos, strand = fields[0], fields[1], int(fields[2]), fields[3]
            seq = find_motif(fasta_file, chrom, pos, strand, motif_length)
            print(f"{line.strip()}\t{seq}", file=output_file)

def main():
    parser = argparse.ArgumentParser(description="Extract motifs from a file using a reference genome.")
    parser.add_argument("input_file", help="Input file (can be gzipped)")
    parser.add_argument("output_file", help="Output file (will be gzipped if ends with .gz)")
    parser.add_argument("motif_length", type=int, help="Length of the motif to extract")
    parser.add_argument("reference_fasta", help="Reference genome FASTA file")

    args = parser.parse_args()
    extract_motifs(args.input_file, args.output_file, args.motif_length, args.reference_fasta)

if __name__ == "__main__":
    main()
