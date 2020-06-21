#include "Python.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "common/reader.c"
#include "common/math.c"

#include "python/load.c"
#include "python/run.c"
#include "python/reset.c"

static PyMethodDef methods[] = {
    {"load", (PyCFunction) load, METH_VARARGS | METH_KEYWORDS, "Read matrix to memory"},
    {"run", (PyCFunction) run, METH_VARARGS | METH_KEYWORDS, "Run GeneQuery"},
    {"reset", (PyCFunction) reset, METH_VARARGS | METH_KEYWORDS, "Clean up memory"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef genequery = {
    PyModuleDef_HEAD_INIT, "genequery", NULL, -1, methods
};

PyMODINIT_FUNC PyInit_genequery(void)
{
    char * fname = "./data/log_factorial.dat";
    factorials(fname, 100001);
    return PyModule_Create(&genequery);
}
