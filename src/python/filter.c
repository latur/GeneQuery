static PyObject * filter(PyObject *self, PyObject *args, PyObject *keywds)
{
    PyObject * offsets;
    int index = -1;
    unsigned int gse;
    static char * kwlist[] = {"offsets", "index", "gse", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "Oii|", kwlist,
                                     &offsets, &index, &gse)) return NULL;
    /* --------------------------------------------------------------------- */
    if (index < 0 || index > dbi) return NULL;

    PyObject * result = Py_BuildValue("[]");
    gmt * db = dbset[index];

    unsigned * relation = NULL;
    relation = (unsigned *) malloc(sizeof(unsigned) * db->gse_count * 1);
    memset(relation, 0, sizeof(unsigned) * db->gse_count * 1);

    PyObject * num;
    unsigned size = PyList_Size(offsets);
    for (unsigned gene = 0; gene < size; gene++) {
        num = PyList_GetItem(offsets, gene);
        long offset = PyLong_AsLong(num);
        addRelation(relation, db, offset, 0);
        if (relation[gse] < 253) PyList_Append(result, Py_BuildValue("ii", relation[gse], offset));
    }

    free(relation);
    return result;
}
