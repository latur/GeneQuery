
static PyObject * reset(PyObject *self, PyObject *args, PyObject *keywds)
{
    while (dbi > 0) {
        dbi--;
        free(dbset[dbi]->genes);
        free(dbset[dbi]->counts);
        free(dbset[dbi]);
    }
    Py_RETURN_NONE;
}
