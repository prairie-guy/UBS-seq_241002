#!/usr/bin/env python3
import argparse
from pysam import FastaFile
from xopen import xopen
from seqpy import revcomp

def find_motif(fasta_file, chrom, pos, strand, motif_length):
    try:
        if strand == "+":
            seq = fasta_file.fetch(chrom, pos-1, pos+motif_length-1)
        else:
            seq = revcomp(fasta_file.fetch(chrom, pos-motif_length, pos))
        return seq
    except (ValueError, KeyError):
        return "N" * motif_length

def extract_motif(input_file, output_file, reference_fasta, motif_len, field_idxs):
    chrom_idx, pos_idx, strand_idx = [index - 1 for index in field_idxs] # Convert 1-based indices to 0-based indices
    with xopen(input_file) as input_file, xopen(output_file, "w") as output_file, FastaFile(reference_fasta) as fasta_file:
        header = next(input_file).strip()
        print(f"{header}\tMotif", file=output_file)
        for line in input_file:
            line = line.strip()
            fields = line.split('\t')
            chrom, pos, strand = fields[chrom_idx], int(fields[pos_idx]), fields[strand_idx]
            motif = find_motif(fasta_file, chrom, pos, strand, motif_len)
            print(f"{line}\t{motif}", file=output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract motifs from a file using a reference genome.")
    parser.add_argument("input_file", help="Input file (can be gzipped)")
    parser.add_argument("output_file", help="Output file (will be gzipped if ends with .gz)")
    parser.add_argument("reference_fasta", help="Reference genome FASTA file")
    parser.add_argument("--motif_length", type=int, default=3, help="Length of the motif to extract (default: 3)")
    parser.add_argument("--field-indices", nargs=3, type=int, default=[2, 3, 4],
                        help="Indices for chromosome, position, and strand fields (1-based, default: 2, 3, 4)")

    args = parser.parse_args()
    extract_motif(args.input_file, args.output_file, args.reference_fasta, args.motif_length, args.field_indices)
