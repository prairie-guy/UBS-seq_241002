#!/usr/bin/env python3

from Bio import SeqIO

def abi2peaks(file_abi, file_out=None):
    ab = SeqIO.read(file_abi,'abi')
    ab.annotations['abif_raw']['P1AM1']
    return list(enumerate( zip(map(chr, ab.annotations['abif_raw']['PBAS1']),
                               map(chr, ab.annotations['abif_raw']['P2BA1']),
                               ab.annotations['abif_raw']['P1AM1'],
                               ab.annotations['abif_raw']['P2AM1'])))
