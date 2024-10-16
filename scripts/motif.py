#!/usr/bin/env python3

# C. Bryan Daniels
# 8/28/24

# Based on code by Chang Ye: https://github.com/y9c/variant

import argparse, os
from pysam import FastaFile
from xopen import xopen
from seqpy import revcomp

def find_motif(fasta_idx: FastaFile, chrom: str, pos: int, strand: str, motif_len: int) -> str:
    """
    Takes a 'fasta_idx: FastaFile' and returns a 'motif', keyed by ['chrom', 'pos', 'strand']
    Depends upon seqpy.revcomp compiled from seqpy.c
    """
    try:
        if strand == "+":
            seq = fasta_idx.fetch(chrom, pos-1, pos+motif_len-1)
        else:
            seq = revcomp(fasta_idx.fetch(chrom, pos-motif_len, pos))
        return seq
    except (ValueError, KeyError):
        return "N" * motif_len

def append_motif(in_tsv: str, out_tsv: str, fasta_reference: str, motif_len=3, field_idxs=[2,3,4])-> bool:
    """
    Append a column field 'Motif' to 'in_tsv' file and saves as 'out_tsv'
        - 'in_tsv'  can be tsv or tsv.gz, where the ext  does not matter as function determines 'in_tsv' format
        - 'out_tsv' can be tsv or tsv.gz, where gz ext determines gzipped format of 'out_tsv'
        - 'fasta_reference' is a well formatted fasta file with single or multiple chromosomes
        - `motif_len` is motif length, starting with each 'C' on (+) or (-) strand, default=3
        - 'field_indx' is a 1-based position list of fields for 'chrom','pos','strand' in 'in_tsv', default=[2,3,4]
          Note: Actual name of fields do not matter, only field positions
        - Returns True is sucessful and False otherwise
    """
    # Uncomment to return file.tsv.gz vs file.tsv
    # out_tsv = check_ext(in_tsv, out_tsv)
    if in_tsv == out_tsv:
        print(f'Warning: {in_tsv} and {out_tsv} are the same file')
        return(False)
    chrom_idx, pos_idx, strand_idx = [index - 1 for index in field_idxs] # Convert 1-based indices to 0-based indices
    with xopen(in_tsv) as in_stream, xopen(out_tsv, "w") as out_stream, FastaFile(fasta_reference) as fasta_idx:
        header = next(in_stream).strip()
        if 'Motif' in header or 'motif' in header:
            print(f'Warning: Motif previously appended to {in_tsv}' )
            os.remove(out_tsv)
            return(None)
        print(f"{header}\tMotif", file=out_stream)
        for line in in_stream:
            line = line.strip()
            fields = line.split('\t')
            chrom, pos, strand = fields[chrom_idx], int(fields[pos_idx]), fields[strand_idx]
            motif = find_motif(fasta_idx, chrom, pos, strand, motif_len)
            print(f"{line}\t{motif}", file=out_stream)
    return(True)

def is_gzipped(file):
    with open(file, 'rb') as f: return f.read(2) == b'\x1f\x8b'

def check_ext(in_tsv,out_tsv):
    "Use to return file.tsv.gz vs file.tsv"
    if is_gzipped(in_tsv):
        return out_tsv if out_tsv.endswith('.gz') else f'{out_tsv}.gz'
    return out_tsv.removesuffix('.gz') if out_tsv.endswith('.gz') else out_tsv

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Append a column field 'Motif' to 'in_tsv' file and save as 'out_tsv'")
    parser.add_argument("in_tsv", help="'in_tsv'  can be tsv or tsv.gz, where the ext  does not matter as function determines 'in_tsv' format")
    parser.add_argument("out_tsv", help="'out_tsv' can be tsv or tsv.gz, where gz ext determines gzipped format of 'out_tsv'")
    parser.add_argument("fasta_reference", help="'fasta_reference' is a well formatted fasta file with single or multiple chromosomes")
    parser.add_argument("--motif_len", type=int, default=3, help="`motif_len` is motif length, starting with each 'C' on (+) or (-) strand, default=3")
    parser.add_argument("--fields", type=lambda x: [int(i) for i in x.split(',')], default=[2,3,4],
                        help="'field_indx' is a comma delimited position list of fields for 'chrom','pos','strand' in 'in_tsv', default=2,3,4")
    args = parser.parse_args()
    append_motif(args.in_tsv, args.out_tsv, args.fasta_reference, args.motif_len, args.fields)
