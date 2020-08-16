static PyObject * filter(PyObject *self, PyObject *args, PyObject *keywds)
{
    PyObject * offsets;
    int index = -1;
    unsigned int gse;
    unsigned int module;
    static char * kwlist[] = {"offsets", "index", "gse", "module", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "Oiii|", kwlist,
                                     &offsets, &index, &gse, &module)) return NULL;
    /* --------------------------------------------------------------------- */
    if (index < 0 || index > dbi) return NULL;

    PyObject * result = Py_BuildValue("[]");
    gmt * db = dbset[index];

    unsigned * relation = NULL;
    relation = (unsigned *) malloc(sizeof(unsigned) * db->gse_count * 1);
    memset(relation, 0, sizeof(unsigned) * db->gse_count * 1);

    // Genes -> [modules block]:
    PyObject * num;
    unsigned size = PyList_Size(offsets);
    for (unsigned gene = 0; gene < size; gene++) {
        num = PyList_GetItem(offsets, gene);
        long offset = PyLong_AsLong(num);
        addRelation(relation, db, offset, 0);
        if (relation[gse] == module) PyList_Append(result, num);
    }

    free(relation);
    return result;
}
