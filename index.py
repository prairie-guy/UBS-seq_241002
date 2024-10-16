#!/usr/bin/env python
import os
import gzip
import sys

def process_gz_files(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".gz"):
            print(filename)
            try:
                with gzip.open(os.path.join(directory, filename), 'rt') as f:
                    print(f.readline().strip())
            except Exception as e:
                print("skipping")
            print("=")

if __name__ == "__main__":
    if len(sys.argv)!= 2:
        print("Usage: python script_name.py <directory>")
        sys.exit(1)
    process_gz_files(sys.argv[1])
