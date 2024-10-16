#!/usr/bin/env python3

# C. Bryan Daniels
# 8/28/24

# Based on code by Chang Ye: https://github.com/y9c/variant

import argparse
from pysam import FastaFile
from xopen import xopen
from seqpy import revcomp

def find_motif(fasta_idx, chrom, pos, strand, motif_len):
    try:
        if strand == "+":
            seq = fasta_idx.fetch(chrom, pos-1, pos+motif_len-1)
        else:
            seq = revcomp(fasta_idx.fetch(chrom, pos-motif_len, pos))
        return seq
    except (ValueError, KeyError):
        return "N" * motif_len

def append_motif(in_tsv, out_tsv, fasta_reference, motif_len, field_idxs):
    # Uncomment to return file.tsv.gz vs file.tsv
    # out_tsv = check_ext(in_tsv, out_tsv)
    chrom_idx, pos_idx, strand_idx = [index - 1 for index in field_idxs] # Convert 1-based indices to 0-based indices
    with xopen(in_tsv) as in_tsv, xopen(out_tsv, "w") as out_tsv, FastaFile(fasta_reference) as fasta_idx:
        header = next(in_tsv).strip()
        print(f"{header}\tMotif", file=out_tsv)
        for line in in_tsv:
            line = line.strip()
            fields = line.split('\t')
            chrom, pos, strand = fields[chrom_idx], int(fields[pos_idx]), fields[strand_idx]
            motif = find_motif(fasta_idx, chrom, pos, strand, motif_len)
            print(f"{line}\t{motif}", file=out_tsv)

def is_gzipped(file):
    with open(file, 'rb') as f: return f.read(2) == b'\x1f\x8b'

def check_ext(in_tsv,out_tsv):
    "Use to return file.tsv.gz vs file.tsv"
    if is_gzipped(in_tsv):
        return out_tsv if out_tsv.endswith('.gz') else f'{out_tsv}.gz'
    return out_tsv.removesuffix('.gz') if out_tsv.endswith('.gz') else out_tsv

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the motif in each line of in_tsv and append it to the line based upon a fasta_reference file")
    parser.add_argument("in_tsv", help="Input file (can be gzipped)")
    parser.add_argument("out_tsv", help="Output file (will be gzipped if ends with .gz)")
    parser.add_argument("fasta_reference", help="Reference genome FASTA file")
    parser.add_argument("--motif_len", type=int, default=3, help="Length of the motif to extract (default: 3)")
    parser.add_argument("--fields", type=lambda x: [int(i) for i in x.split(',')], default=[2,3,4],
                        help="Indices for chromosome, position, and strand fields (1-based, comma-separated, default: 2,3,4)")
    args = parser.parse_args()
    append_motif(args.in_tsv, args.out_tsv, args.fasta_reference, args.motif_len, args.fields)
