#!/usr/bin/env python3

import sys,re,argparse
import pandas as pd
from pathlib import Path

def merge_sample_fname(in_csv: str, dir_path: str, out_csv: str = 'sample_fname.csv') -> None:
    """
    Merge sample information from a CSV file with filename data from a specified directory.

    Args:
    in_csv (str): Path to the input CSV file containing sample information.
                  This file must include 'SampleID' and 'Content' columns.
    dir_path (str): Path to the directory containing the .fastq.gz or .fq.gz files.
                    Sequences in this directory are expected to contain the SampleID.
    out_csv (str, optional): Path where the output CSV file will be written.
                             Defaults to 'sample_fname.csv' in the current directory.

    Returns:
    None: The function writes the results to a CSV file and prints a confirmation message.
    """

    dir_path_obj = Path(dir_path)
    files = list(dir_path_obj.glob('*.[ff][aq]*.gz'))
    pattern = r'-(.+?)_.*_(R[12])_'   # Sequencer fname style
    #pattern = r'([A-Za-z0-9]+)_(R[12])\.fq\.gz$' # Standard style, exmple A1_R2.fq.gz
    files_data = [{'SampleID': match.group(1), 'Std': match.group(2), 'Sequence': f"{dir_path_obj.name}/{file.name}"}
                  for file in files if (match := re.search(pattern, file.name))]


    in_path = Path(in_csv).parent
    (pd.read_csv(in_csv)
        .merge(pd.DataFrame(files_data), on='SampleID', how='inner')
        .assign(FullID=lambda df: df['SampleID'] + '_' + df['Std'])
     [[ 'SampleID', 'Std', 'FullID', 'Content', 'Sequence']]
        .to_csv(f'{in_path}/{out_csv}', index=False))

    print(f"Results written to {in_path}/{out_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge sample information with filename data.')
    parser.add_argument('in_csv', type=str, help='Path to the input CSV file (must contain SampleID and Content columns)')
    parser.add_argument('dir_path', type=str, help='Path to the directory containing the files (filenames should include SampleID)')
    parser.add_argument('--out_csv', type=str, default='sample_fname.csv', help='Path to the output CSV file (default: sample_fname.csv)')
    args = parser.parse_args()
    if not Path(args.in_csv).is_file(): parser.error(f"Error: {args.in_csv} is not a valid file")
    if not Path(args.dir_path).is_dir(): parser.error(f"Error: {args.dir_path} is not a valid directory")

    merge_sample_fname(args.in_csv, args.dir_path, args.out_csv)
