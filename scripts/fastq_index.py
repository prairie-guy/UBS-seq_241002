#!/usr/bin/env python3

import argparse, gzip, re, sys
import pandas as pd
from pathlib import Path
from typing import Union, List


# Dual Indexes
# https://raw.githubusercontent.com/chaochungkuo/GPM/436ba85d2b94f1b79b15beb5a58b9d5dd6c65654/demultiplex/index_sequences/NEB_Dual_Index_E6440S.txt

# Listing of more indexes
# https://github.com/chaochungkuo/GPM/tree/436ba85d2b94f1b79b15beb5a58b9d5dd6c65654/demultiplex/index_sequences

adapter_index = {
    # NEB Single Index
    # "ATCACG": 1, "CGATGT": 2, "TTAGGC": 3, "TGACCA": 4, "ACAGTG": 5,
    # "GCCAAT": 6, "CAGATC": 7, "ACTTGA": 8, "GATCAG": 9, "TAGCTT": 10,
    # "GGCTAC": 11, "CTTGTA": 12, "AGTCAA": 13, "AGTTCC": 14, "ATGTCA": 15,
    # "CCGTCC": 16, "GTAGAG": 17, "GTCCGC": 18, "GTGAAA": 19, "GTGGCC": 20,
    # "GTTTCG": 21, "CGTACG": 22, "GAGTGG": 23, "GGTAGC": 24, "ACTGAT": 25,
    # "ATGAGC": 26, "ATTCCT": 27, "CAAAAG": 28, "CAACTA": 29, "CACCGG": 30,
    # "CACGAT": 31, "CACTCA": 32, "CAGGCG": 33, "CATGGC": 34, "CATTTT": 35,
    # "CCAACA": 36, "CGGAAT": 37, "CTAGCT": 38, "CTATAC": 39, "GTGATC": 40,
    # "GACGAC": 41, "TAATCG": 42, "TACAGC": 43, "TATAAT": 44, "TCATTC": 45,
    # "TCCCGA": 46, "TCGAAG": 47, "TCGGCA": 48,
    # NEB Dual Index: Reverse Complement Workflow
    # "CGTATTCG": "A1", "TCAAGGAC": "B1", "AAGCACTG": "C1", "GCAATGGA": "D1",
    # "CAATCGAC": "E1", "GGCGTTAT": "F1", "GTTAAGGC": "G1", "CCTATACC": "H1",
    # "CTCCTAGA": "A2", "GTTACGCA": "B2", "CTAGCAAG": "C2", "ATCTCGCT": "D2",
    # "GTGCCATA": "E2", "GGTGATTC": "F2", "CACCTTAC": "G2", "TTCTCTCG": "H2",
    # "TAGTTGCG": "A3", "AGTCTGTG": "B3", "TGCTTCCA": "C3", "GGCTATTG": "D3",
    # "TGTTCGAG": "E3", "AACTTGCC": "F3", "TGGTAGCT": "G3", "GTATGCTG": "H3",
    # "GAGATACG": "A4", "GCACGTAA": "B4", "GCTTAGCT": "C4", "GGTGTCTT": "D4",
    # "TGGAGTTG": "E4", "GCAAGATC": "F4", "CAGTGAAG": "G4", "AAGTCGAG": "H4",
    # "AGGTGTAC": "A5", "AACCTTGG": "B5", "AACCGTTC": "C5", "TCAACTGG": "D5",
    # "ACGATGAC": "E5", "TCGCATTG": "F5", "GTTCAACC": "G5", "AACCGAAG": "H5",
    # "TAATGCCG": "A6", "ATTGCGTG": "B6", "GACATTCC": "C6", "CTTCACCA": "D6",
    # "TGATGTCC": "E6", "TGTACACC": "F6", "TGGCTATC": "G6", "TGTTGTGG": "H6",
    # "GTCGGTAA": "A7", "TCAGACGA": "B7", "ACCTGGAA": "C7", "AGACCGTA": "D7",
    # "ACGGTCTT": "E7", "TGAACCTG": "F7", "AGCTCCTA": "G7", "CTGGAGTA": "H7",
    # "AGGTCACT": "A8", "GATAGGCT": "B8", "GGAGATGA": "C8", "GATACTGG": "D8",
    # "TCTCGCAA": "E8", "CTTCGTTC": "F8", "GCAATTCG": "G8", "TCTCTTCC": "H8",
    # "GAATCCGA": "A9", "TGGTACAG": "B9", "GTACTCTC": "C9", "TGCGTAGA": "D9",
    # "GGAATTGC": "E9", "CTTCTGAG": "F9", "CTTAGGAC": "G9", "TCTAACGC": "H9",
    # "GTACCTTG": "A10", "CAAGGTCT": "B10", "GTAACGAC": "C10", "TCGGTTAC": "D10",
    # "ACGGATTC": "E10", "TGCTCATG": "F10", "GTCCTAAG": "G10", "GGTCAGAT": "H10",
    # "CATGAGGA": "A11", "GCTATCCT": "B11", "ATTCCTCC": "C11", "ATGACGTC": "D11",
    # "TTAAGCGG": "E11", "AGTTCGTC": "F11", "AACGTGGA": "G11", "CTCTGGTT": "H11",
    # "TGACTGAC": "A12", "ATGGAAGG": "B12", "GTGTTCCT": "C12", "GCTGTAAG": "D12",
    # "TGCAGGTA": "E12", "TAGCGTCT": "F12", "CTGTGTTG": "G12", "TGTGGTAC": "H12",
    # NEB Single Index
    "ATCACG": '1', "CGATGT": '2', "TTAGGC": '3', "TGACCA": '4', "ACAGTG": '5',
    "GCCAAT": '6', "CAGATC": '7', "ACTTGA": '8', "GATCAG": '9', "TAGCTT": '10',
    "GGCTAC": '11', "CTTGTA": '12', "AGTCAA": '13', "AGTTCC": '14', "ATGTCA": '15',
    "CCGTCC": '16', "GTAGAG": '17', "GTCCGC": '18', "GTGAAA": '19', "GTGGCC": '20',
    "GTTTCG": '21', "CGTACG": '22', "GAGTGG": '23', "GGTAGC": '24', "ACTGAT": '25',
    "ATGAGC": '26', "ATTCCT": '27', "CAAAAG": '28', "CAACTA": '29', "CACCGG": '30',
    "CACGAT": '31', "CACTCA": '32', "CAGGCG": '33', "CATGGC": '34', "CATTTT": '35',
    "CCAACA": '36', "CGGAAT": '37', "CTAGCT": '38', "CTATAC": '39', "GTGATC": '40',
    "GACGAC": '41', "TAATCG": '42', "TACAGC": '43', "TATAAT": '44', "TCATTC": '45',
    "TCCCGA": '46', "TCGAAG": '47', "TCGGCA": '48',
    # NEB Dual Index: i7 Index
    "TTACCGAC": "A1", "TCGTCTGA": "B1", "TTCCAGGT": "C1", "TACGGTCT": "D1",
    "AAGACCGT": "E1", "CAGGTTCA": "F1", "TAGGAGCT": "G1", "TACTCCAG": "H1",
    "AGTGACCT": "A2", "AGCCTATC": "B2", "TCATCTCC": "C2", "CCAGTATC": "D2",
    "TTGCGAGA": "E2", "GAACGAAG": "F2", "CGAATTGC": "G2", "GGAAGAGA": "H2",
    "TCGGATTC": "A3", "CTGTACCA": "B3", "GAGAGTAC": "C3", "TCTACGCA": "D3",
    "GCAATTCC": "E3", "CTCAGAAG": "F3", "GTCCTAAG": "G3", "GCGTTAGA": "H3",
    "CAAGGTAC": "A4", "AGACCTTG": "B4", "GTCGTTAC": "C4", "GTAACCGA": "D4",
    "GAATCCGT": "E4", "CATGAGCA": "F4", "CTTAGGAC": "G4", "ATCTGACC": "H4",
    "TCCTCATG": "A5", "AGGATAGC": "B5", "GGAGGAAT": "C5", "GACGTCAT": "D5",
    "CCGCTTAA": "E5", "GACGAACT": "F5", "TCCACGTT": "G5", "AACCAGAG": "H5",
    "GTCAGTCA": "A6", "CCTTCCAT": "B6", "AGGAACAC": "C6", "CTTACAGC": "D6",
    "TACCTGCA": "E6", "AGACGCTA": "F6", "CAACACAG": "G6", "GTACCACA": "H6",
    "CGAATACG": "A7", "GTCCTTGA": "B7", "CAGTGCTT": "C7", "TCCATTGC": "D7",
    "GTCGATTG": "E7", "ATAACGCC": "F7", "GCCTTAAC": "G7", "GGTATAGG": "H7",
    "TCTAGGAG": "A8", "TGCGTAAC": "B8", "CTTGCTAG": "C8", "AGCGAGAT": "D8",
    "TATGGCAC": "E8", "GAATCACC": "F8", "GTAAGGTG": "G8", "CGAGAGAA": "H8",
    "CGCAACTA": "A9", "CACAGACT": "B9", "TGGAAGCA": "C9", "CAATAGCC": "D9",
    "CTCGAACA": "E9", "GGCAAGTT": "F9", "AGCTACCA": "G9", "CAGCATAC": "H9",
    "CGTATCTC": "A10", "TTACGTGC": "B10", "AGCTAAGC": "C10", "AAGACACC": "D10",
    "CAACTCCA": "E10", "GATCTTGC": "F10", "CTTCACTG": "G10", "CTCGACTT": "H10",
    "GTACACCT": "A11", "CCAAGGTT": "B11", "GAACGGTT": "C11", "CCAGTTGA": "D11",
    "GTCATCGT": "E11", "CAATGCGA": "F11", "GGTTGAAC": "G11", "CTTCGGTT": "H11",
    "CGGCATTA": "A12", "CACGCAAT": "B12", "GGAATGTC": "C12", "TGGTGAAG": "D12",
    "GGACATCA": "E12", "GGTGTACA": "F12", "GATAGCCA": "G12", "CCACAACA": "H12"
}


def fastq_index(path: Union[str, Path], suffix: str = 'gz') -> Union[List[str], pd.DataFrame]:
    """
    Examine fastq.gz formatted files and determine the Index used during library preparation.
    Handles both single and dual indexes from NEB.
    Returns a list for single files or a pandas DataFrame for multiple files.
    """
    def process_file(file_path):
        try:
            match = re.search(r'-(.+?)_.*_(R[12])_', file_path.name)
            sample_id, std = match.groups() if match else ('0', '0')
            std = std[-1]  # Ensure std is only '1' or '2'
            full_id = f'{sample_id}_R{std}' if sample_id != '0' else '0'

            with gzip.open(file_path, 'rt') as f:
                # LH and RH keys are both 8 ch long
                # NEB Single Index is first 6 ch of LH key
                # NEB Dual Indexes are either full LH or RH
                # There should only be a single idx returned
                key1, key2 = next(f).strip().split(':')[-1].split('+')
                key0 = key1[0:6]
                idx0, idx1, idx2 = adapter_index.get(key0, '0') ,adapter_index.get(key1, '0'),adapter_index.get(key2, '0')
                #print([idx0, idx1, idx2])
                indx = idx1 if idx1 != '0' else idx2
                if indx == '0': indx = idx0
            return [sample_id, indx, full_id, str(file_path)]
        except Exception:
            return ['0', '0', '0', str(file_path)]

    path = Path(path)
    if path.is_file():
        return process_file(path)
    files = list(path.glob(f'*.{suffix}'))
    if not files:
        return pd.DataFrame(columns=['SampleID', 'Index', 'FullID', 'Path'])
    results = [process_file(f) for f in files]
    return pd.DataFrame(results, columns=['SampleID', 'Index', 'FullID', 'Path']).sort_values('FullID')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FASTQ files and determine the Index used during library preparation.")
    parser.add_argument("path", help="Path to a single FASTQ file or a directory containing FASTQ files")
    parser.add_argument("--suffix", default="gz", help="Suffix for FASTQ files (default: gz)")
    args = parser.parse_args()

    result = fastq_index(args.path, args.suffix)

    if isinstance(result, list):
        print("\t".join(result))
    elif isinstance(result, pd.DataFrame):
        if not result.empty:
            result.to_csv(sys.stdout, sep='\t', index=False)
        else:
            print("No matching files found.")
