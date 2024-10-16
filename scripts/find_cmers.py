#!/usr/bin/env python3

# Example Usage: Find all CG positions on pUC19.fa for both positive and negative strands.

import argparse, re

def rev_comp(dna):
    "Find reverse complement of DNA with only: 'A','T','C','G'"
    #rc ={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    rc = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "[^G]",
    "[^G]": "H",
    "B": "V",
    "N": "N",
    ".": ".",
    "-": "-",
}
    return ''.join([rc[nt.capitalize()] for nt in dna])


def fasta2seq(fasta_file, chr=None):
    sections = open(fasta_file, 'r').read().split('>')
    if chr is not None:
        matched_section = next((s for s in sections if re.match(rf"^{chr}\b", s)), None)
        if matched_section:
            return '\n'.join(matched_section.splitlines()[1:]).replace('\n', '')
        print(f"Error: No sequence found matching chr: {chr}")
        return None
    else:
        if len(sections) > 2:  # Check if there are multiple sections
            print("Error: Multiple sections found in the FASTA file; please specify a chromosome header with --chr.")
            return None
        return '\n'.join(sections[1].splitlines()[1:]).replace('\n', '')

def find_cmers(seq, cmer, c_offset = 0, std='both'):
    """
    Returns the positions (1-based) of the 'C' in a cmer of either 'pos', 'neg' or both (default) strands.
    If the target 'C' is not at start of cmer, then set `c_offset` (default=0)
    """
    if not all([n in 'ATCG' for n in cmer]) or 'C' not in cmer:
        return('usage: cmer must contain only "ATCG" and at least on "C" ')
    if cmer[c_offset] != 'C':
        return(f'usage: "No C" located in cmer={cmer} at "c_offset={c_offset}"')
    #print(seq)
    #print(rev_comp(seq))

    # Note: (?=...) is a lookahead assertion matching characters w/o consuming, allowing for overlapping sequences
    # cmer[::-1] -> rev(cmer)
    pos = set([1 + c_offset + match.start() for match in re.finditer(f'(?={cmer})',seq)])
    neg = set([len(cmer) - c_offset + match.start() for match in re.finditer(f'(?={cmer[::-1]})',rev_comp(seq))])
    if std == 'pos': return pos
    if std == 'neg': return neg
    return pos | neg


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Returns positions of cmers positive, negative or both strands')
    parser.add_argument('fasta_file', type=str, help='Path to the FASTA file')
    parser.add_argument('cmer', type=str, help='Cmer sequence to search for')
    parser.add_argument('--c_offset', type=int, default=0, help='0-base Offset value for C position (default: 0)')
    parser.add_argument('--std', type=str, default='Both', help='Find cmers on pos, neg or both strands. (default: both)')
    parser.add_argument('--chr', type=str, default=None, help='Specify chromosome header to match (default: None)')
    args = parser.parse_args()

    fasta_file, cmer, c_offset, std, chr = args.fasta_file, args.cmer, args.c_offset, args.std, args.chr
    sequence = fasta2seq(fasta_file, chr)
    if sequence is not None:
        for position in sorted(find_cmers(sequence, cmer, c_offset, std)):
            print(position)
