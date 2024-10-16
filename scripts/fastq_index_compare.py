#!/usr/bin/env python3

from fastq_index import *
import warnings

def fastq_index_compare(path_dir: Union[str, Path], path_csv: Union[str, Path]) -> tuple[bool, pd.DataFrame]:
    """
    Compare actual indexes from FASTQ files to expected indexes in a CSV file.
    Returns a tuple (success, df), where success is a bool: True if all indexes match, False otherwise.
    If success is True, returns the full merged DataFrame. Otherwise, returns a DataFrame of unmatched rows.
    """
    path_dir = Path(path_dir)
    path_csv = Path(path_csv)
    if not path_dir.is_dir():
        raise ValueError(f"Directory not found: {path_dir}")
    if not path_csv.is_file():
        raise ValueError(f"CSV file not found: {path_csv}")

    df = fastq_index(path_dir)
    dfi = pd.read_csv(path_csv, dtype=str)
    dfm = df.merge(dfi, on='SampleID')

    if 'Index_y' not in dfm.columns:
        warnings.warn("'Index_y' column not found in the CSV file!! Returning the unsuccessfully merged df...")
        return False, dfm

    unmatched = dfm[dfm['Index_x'] != dfm['Index_y']]
    if unmatched.empty:
        return True, dfm
    else:
        warnings.warn("Unmatched Indexes!! Returning a df of unmatched rows...")
        return False, unmatched

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare actual indexes from FASTQ files in a directory to those in an experiment CSV file.")
    parser.add_argument("path_dir", help="Path to the directory containing FASTQ files")
    parser.add_argument("path_csv", help="Path to the experiment CSV file containing expected indexes")
    parser.add_argument("--merge", action="store_true", help="Print the full merged dataframe if successful")
    args = parser.parse_args()

    try:
        success, result = fastq_index_compare(args.path_dir, args.path_csv)

        if success:
            print("Success: All indexes match.", file=sys.stderr)
            if args.merge:
                result.to_csv(sys.stdout, sep='\t', index=False)
        else:
            if 'Index_y' not in result.columns:
                print("Warning: 'Index_y' column not found in the CSV file. Cannot perform comparison.", file=sys.stderr)
            else:
                print("Warning: Some rows have mismatched indexes.", file=sys.stderr)
                mismatched = result[result['Index_x'] != result['Index_y']]
                mismatched.to_csv(sys.stdout, sep='\t', index=False)

    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)
