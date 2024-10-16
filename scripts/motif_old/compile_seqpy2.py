#!/usr/bin/env python3
import os
import sys
import shutil
import glob
from setuptools import setup, Extension

def write_c_file():
    c_code = '''
#define PY_SSIZE_T_CLEAN 1
#include <Python.h>
#include <stdint.h>
#include <stdlib.h>

static unsigned char seq_comp_table[256] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,
    15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
    30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,
    45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,
    60,  61,  62,  63,  64,  'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J',
    'M', 'L', 'K', 'N', 'O', 'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R',
    'Z', 91,  92,  93,  94,  95,  64,  't', 'v', 'g', 'h', 'e', 'f', 'c', 'd',
    'i', 'j', 'm', 'l', 'k', 'n', 'o', 'p', 'q', 'y', 's', 'a', 'a', 'b', 'w',
    'x', 'r', 'z', 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134,
    135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
    150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
    165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
    180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
    195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
    210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224,
    225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254,
    255};

static PyObject* seqpy_revcomp(PyObject* self, PyObject* args) {
  PyObject* r;
  char *seq, *rev;
  Py_ssize_t i, len;
  PyArg_ParseTuple(args, "s#", &seq, &len);
  rev = (char*)malloc(len);
  for (i = 0; i < len; ++i) rev[len - i - 1] = seq_comp_table[(uint8_t)seq[i]];
  r = Py_BuildValue("s#", rev, len);
  free(rev);
  return r;
}

static PyMethodDef seqpy_methods[] = {
    {"revcomp",
     seqpy_revcomp,
     METH_VARARGS,
     "Reverse complement a DNA sequence"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef seqpy_module =
    {PyModuleDef_HEAD_INIT, "seqpy", NULL, -1, seqpy_methods};

PyMODINIT_FUNC PyInit_seqpy(void) { return PyModule_Create(&seqpy_module); }
'''
    with open('seqpy.c', 'w') as f:
        f.write(c_code)

def build_extension():
    # Define the extension
    module = Extension('seqpy', 
                       sources=['seqpy.c'],
                       extra_compile_args=['-O3'])  # Keep the optimization flag
    
    # Temporarily redirect stdout and stderr
    old_stdout, old_stderr = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = open('build.log', 'w'), open('build.log', 'w')
    
    try:
        # Run the setup
        setup(
            name='seqpy',
            version='1.0',
            description='Original seqpy package',
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
    
    # Remove C file
    if os.path.exists('seqpy.c'):
        os.remove('seqpy.c')
        print("Removed seqpy.c")

if __name__ == "__main__":
    print("Writing original seqpy.c file...")
    write_c_file()
    print("Building seqpy extension...")
    build_extension()
    print("Build process completed.")
    print("Performing cleanup...")
    cleanup()
    print("\nBuild process completed successfully.")
    print("You can now import seqpy in Python with: import seqpy")