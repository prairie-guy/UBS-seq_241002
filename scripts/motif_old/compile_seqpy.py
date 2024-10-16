#!/usr/bin/env python3

import os
import sys
import shutil
import glob
from setuptools import setup, Extension

def build_extension():
    # Define the extension
    module = Extension('seqpy', sources=['seqpy.c'])

    # Temporarily redirect stdout and stderr
    old_stdout, old_stderr = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = open('build.log', 'w'), open('build.log', 'w')

    try:
        # Run the setup
        setup(
            name='seqpy',
            version='1.0',
            description='This is a demo package',
            ext_modules=[module],
            script_args=['build_ext', '--inplace']
        )
    finally:
        # Restore stdout and stderr
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout, sys.stderr = old_stdout, old_stderr

def cleanup():
    # Rename the output file
    built_files = glob.glob('seqpy*.so')
    if built_files:
        os.rename(built_files[0], 'seqpy.so')
        print(f"Renamed {built_files[0]} to seqpy.so")

    # Remove build folder
    build_folder = 'build'
    if os.path.exists(build_folder) and os.path.isdir(build_folder):
        shutil.rmtree(build_folder)
        print(f"Removed {build_folder} directory")

    # Remove build log
    if os.path.exists('build.log'):
        os.remove('build.log')
        print("Removed build.log")

if __name__ == "__main__":
    print("Building seqpy extension...")
    build_extension()
    print("Build process completed.")

    print("Performing cleanup...")
    cleanup()

    print("\nBuild process completed successfully.")
    print("You can now import seqpy in Python with: import seqpy")
