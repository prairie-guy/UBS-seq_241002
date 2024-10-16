#!/usr/bin/env python3

import os
import re
from collections import defaultdict
from pathlib import Path
from itertools import dropwhile
from subprocess import run

# Search up path for the first 'reference' location that contains 'fasta'
# Expected structure: reference/fasta/ reference/hisat3n reference/meth
ref_paths = ['reference/fasta', '../reference/fasta', '../../reference/fasta', '../../../reference/fasta', '../../../../reference/fasta']
ref_path = next(dropwhile(lambda p: not Path(p).exists(), ref_paths), None)
if not ref_path: raise RuntimeError('error: Could not locate a valid reference dir')
ref_path = Path(ref_path).parent

def get_ref(ref, type):
    """
    Takes ref or [refs] and a type ={'fa', 'hisat3n', 'meth'} and if ref exists returns ref_path
    Example: get_ref('lambda','fa') -> PosixPath('../../reference/fasta/lambda.fa')
    """
    if isinstance(ref,list): ref = mkrefs(ref)
    if type == 'fa':
        ref = refs2fasta.get(ref) if refs2fasta.get(ref) else ref
        path = ref_path/f'fasta/{ref}.fa'
        if path.exists(): return path
    if type == 'meth':
        path = ref_path/f'meth/{ref}.meth'
        if path.exists(): return path
        with open(path,'w') as fh: print(0,file=fh)
        return path
    if type == 'hisat3n':
        path = ref_path/f'hisat3n/{ref}'
        if path.exists(): return path/ref
    if type == 'gtf':
        path = ref_path/f'gtf/{ref}.gtf'
        if path.exists(): return path
    if type == 'bed':
        path = ref_path/f'bed/{ref}.bed'
        if path.exists(): return path

    return None

def get_chr(fasta_file, as_str=False):
    'Returns a list of chr in fastafile'
    chrs = [line[1:].split()[0] for line in open(fasta_file) if line.startswith('>')]
    return ' '.join(chrs) if as_str else chrs

def mkrefs(refs):
    """
    Takes a list of refs and returns an alpha sorted string conjoining the refs with '_'
    Example: mkrefs(['lambda','pUCI9','5mC164']) -> '5mC164_lambda_pUCI9'
    """
    return '_'.join(sorted(refs))

def split_refs(refs):
    """
    Takes a string of refs conjoined with '_' and returns a list of refs
    Example: split_refs('5mC164_lambda_pUCI9') -> ['5mC164', 'lambda', 'pUCI9']
    """
    return refs.split('_')

def index_refs(refs):
    """
    Takes a list of refs and creates new hisat3n index and saves it to ref_path/hisat-3n/
    The respective fasta files must exist.
    Example: index_refs(['lambda','pUC19','5mC164']) ->
        Saved to ../../reference/fasta/5mC164_lambda_pUC19.fa
    """
    ref_fns = [Path(f'{ref_path}/fasta/{ref}.fa') for ref in sorted(refs)]
    if not all([fn.exists() for fn in ref_fns]):
        raise RuntimeError(f'{refs} fasta files not all found in {ref_path}/fasta')
    idx_ref = mkrefs(refs)
    idx_path =  ref_path/f'hisat3n/{idx_ref}'
    idx_fn = idx_path/idx_ref
    fasta_fns = ','.join(list(map(str,ref_fns)))
    if idx_path.exists():
        return(f'{idx_path} hisat-3n index exists')
    idx_path.mkdir(parents=True, exist_ok=False)
    command = ["hisat-3n-build", "--base-change", "C,T", fasta_fns, idx_fn]
    run(command)
    return(f'{idx_ref} hisat-3n index created: {idx_path}')

def fasta_map(reference):
    """
    Create a global mapping of fasta header names to fasta fiilename
    """
    fasta_path = Path(reference)
    filters = {}
    pattern = re.compile(r'^>(\S*)')
    for filepath in fasta_path.glob('*.fa'):
        with filepath.open('r') as file:
            for line in file:
                match = pattern.match(line)
                if match:
                    key = match.group(1)
                    filters[key] = filepath.stem
    return filters


# Create a global mapping of fasta header names to fasta fiilename
refs2fasta = fasta_map(ref_path/'fasta')

# Creates a global mapping of fasta filenames to a list of headers
fasta2refs = {}
for k, v in refs2fasta.items():
    fasta2refs.setdefault(v, []).append(k)
