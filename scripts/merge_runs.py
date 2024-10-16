#!/usr/bin/env python3

# Use merge_runs to conbine muutiple runs of sequencing data from same experiment
# Matches sample and strand across multiple runs (comprised of dirs) to create a
# single file for downstream processing

from pathlib import Path
import re, os, subprocess
from typing import Dict, List, Union
import argparse

def map_run(run_dir: Union[str, Path]) -> Dict[str, str]:
    """
    Extract sample and strand from filename to create sample_key and returns a mapping of sample_strand to filename.
    
    Args:
    run_dir (Union[str, Path]): The run directory to parse.
    
    Returns:
    Dict[str, str]: A dictionary where keys are 'sample_strand' and values are filenames.
    """
    pattern = r'[^-]+-([A-Za-z0-9]+)_[^_]+_(R[12])_[^.]+\.(?:fastq|fq)\.gz$'
    run_dir = Path(run_dir) if isinstance(run_dir, str) else run_dir
    result = {}
    for file in run_dir.iterdir():
        if file.is_file() and file.suffix == '.gz':
            filename = file.name
            match = re.match(pattern, filename)
            if match:
                sample, strand = match.groups()
                key = f"{sample}_{strand}"
                result[key] = filename
    return result

def map_all_runs(*run_dirs: Union[str, Path]) -> Dict[str, List[str]]:
    """
    Match files across multiple run directories based on their sample_strand key.
    
    Args:
    *run_dirs: Variable number of run directory paths.
    
    Returns:
    Dict[str, List[str]]: A dictionary where keys are 'sample_strand' and values are lists of full file paths.
    """
    all_keys = set()
    dir_results = []
    for run_dir in run_dirs:
        result = map_run(run_dir)
        all_keys.update(result.keys())
        dir_results.append((run_dir, result))
    matched_files = {}
    for key in all_keys:
        matched_files[key] = []
        for run_dir, result in dir_results:
            if key in result:
                matched_files[key].append(str(Path(run_dir) / result[key]))
    return matched_files

def merge_runs(out_dir: Union[str, Path], *run_dirs: Union[str, Path], compress: bool = True) -> None:
    """
    Merge files from multiple run directories based on their sample_strand key.
    
    Args:
    out_dir (Union[str, Path]): The output directory for merged files.
    *run_dirs: Variable number of input run directory paths.
    compress (bool): Whether to compress the output files (default: True).
    
    Raises:
    FileExistsError: If the output directory already exists.
    
    Returns:
    None
    """
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise FileExistsError(f"The '{out_dir}' directory already exists.")
    out_dir.mkdir()
    matched_files = map_all_runs(*run_dirs)
    for key, file_list in matched_files.items():
        if file_list:
            ext = '.'.join(Path(file_list[0]).suffixes[:-1])
            output_file = out_dir / f"{key}{ext}"
            with output_file.open('wb') as outfile:
                for file in file_list:
                    subprocess.run(['zcat', file], stdout=outfile, check=True)
            if compress:
                compressed_file = output_file.with_suffix(output_file.suffix + '.gz')
                subprocess.run(['pigz', '-c', output_file], stdout=compressed_file.open('wb'), check=True)
                output_file.unlink()
                print(f"Merged and compressed: {compressed_file}")
            else:
                print(f"Merged: {output_file}")
    print(f"Files merged successfully in the '{out_dir}' directory.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge sequencing run files across multiple directories.")
    parser.add_argument("out_dir", type=str, help="Output directory for merged files")
    parser.add_argument("run_dirs", nargs="+", type=str, help="Input run directories")
    parser.add_argument("--no-compress", action="store_true", help="Do not compress output files")
    args = parser.parse_args()

    merge_runs(args.out_dir, *args.run_dirs, compress=not args.no_compress)
