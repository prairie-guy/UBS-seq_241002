#define PY_SSIZE_T_CLEAN 1
#include <Python.h>
#include <stdint.h>
#include <string.h>

static const unsigned char seq_comp_table[256] = {
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
    255
};

static PyObject* seqpy_revcomp(PyObject* self, PyObject* args) {
    const char *seq;
    Py_ssize_t len;

    if (!PyArg_ParseTuple(args, "s#", &seq, &len)) {
        return NULL;
    }

    char *rev = (char*)PyMem_Malloc(len);
    if (rev == NULL) {
        return PyErr_NoMemory();
    }

    const char *src = seq + len - 1;
    char *dst = rev;

    while (len >= 4) {
        dst[0] = seq_comp_table[(unsigned char)src[0]];
        dst[1] = seq_comp_table[(unsigned char)src[-1]];
        dst[2] = seq_comp_table[(unsigned char)src[-2]];
        dst[3] = seq_comp_table[(unsigned char)src[-3]];
        src -= 4;
        dst += 4;
        len -= 4;
    }

    while (len-- > 0) {
        *dst++ = seq_comp_table[(unsigned char)*src--];
    }

    PyObject *result = Py_BuildValue("s#", rev, dst - rev);
    PyMem_Free(rev);
    return result;
}

static PyMethodDef seqpy_methods[] = {
    {"revcomp", seqpy_revcomp, METH_VARARGS, "Reverse complement a DNA sequence"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef seqpy_module = {
    PyModuleDef_HEAD_INIT,
    "seqpy",
    NULL,
    -1,
    seqpy_methods
};

PyMODINIT_FUNC PyInit_seqpy(void) {
    return PyModule_Create(&seqpy_module);
}
