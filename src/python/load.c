static PyObject * load(PyObject *self, PyObject *args, PyObject *keywds)
{
    char * dbname;
    if (!PyArg_ParseTuple(args, "s", &dbname)) return NULL;

    gmt * db = reader(dbname);
    if (db == NULL) Py_RETURN_NONE;

    dbset[dbi] = db;
    PyObject *index = Py_BuildValue("i", dbi);
    dbi += 1;

    return index;
}
