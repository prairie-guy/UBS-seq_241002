import os
import sys
import shutil
import re

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <source_directory>")
    sys.exit(1)

source_dir = os.path.abspath(sys.argv[1])
parent_dir = os.path.dirname(source_dir)
dest_dir = os.path.join(parent_dir, os.path.basename(source_dir) + "_renamed")
os.makedirs(dest_dir, exist_ok=True)

def process_identifier(identifier):
    return re.sub(r'([A-Za-z])-(\d+)', r'\1\2', identifier)

for filename in os.listdir(source_dir):
    if filename.endswith('.fastq.gz'):
        match = re.match(r'^(\d+)HS-(.+)(_R\d+)(\.fastq\.gz)$', filename)
        if match:
            number, identifier, read, extension = match.groups()
            new_identifier = process_identifier(identifier)
            new_filename = f"{number}DZ-{new_identifier}_S00{read}_001{extension}"
            shutil.copy2(os.path.join(source_dir, filename), os.path.join(dest_dir, new_filename))
            print(f"Copied and renamed: {filename} -> {new_filename}")
        else:
            print(f"Skipping file (doesn't match pattern): {filename}")

print(f"Renaming complete. Renamed files are in: {dest_dir}")