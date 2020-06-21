static PyObject * run(PyObject *self, PyObject *args, PyObject *keywds)
{
    unsigned index = 0;
    PyObject * offsets;

    if (!PyArg_ParseTuple(args, "iO", &index, &offsets)) return NULL;
    if (index >= dbi) {
        PyErr_SetString(PyExc_RuntimeError, "Target DB.bin not loaded");
        return NULL;
    }

    gmt * db = dbset[index];
    unsigned size = PyList_Size(offsets);

    unsigned * relation = NULL;
    relation = (unsigned *) malloc(sizeof(unsigned) * db->gse_count * size);
    memset(relation, 0, sizeof(unsigned) * db->gse_count * size);

    // Genes -> [modules block]:
    PyObject * num;
    for (unsigned gene = 0; gene < size; gene++) {
        num = PyList_GetItem(offsets, gene);
        long offset = PyLong_AsLong(num);

        long cursor = 0;
        unsigned space = 0;
        unsigned value;

        while (cursor < db->gse_count) {
            value = db->genes[offset];
            offset += 1;

            if (value == 254) {
                space = db->genes[offset] + 3;
                offset += 1;
            }
            if (value == 253) {
                space = db->genes[offset] * 255;
                offset += 1;
                space += db->genes[offset] + 3;
                offset += 1;
            }
            while (space > 0) {
                relation[gene * db->gse_count + cursor] = 255;
                cursor += 1;
                space -= 1;
            }
            if (value <= 252 || value == 255) {
                relation[gene * db->gse_count + cursor] = value;
                cursor += 1;
            }
        }
    }

    // Counters matrix: [GSE x modules] -> genes found
    unsigned * counters = NULL;
    counters = (unsigned *) malloc(sizeof(unsigned) * db->gse_count * db->module_count);
    memset(counters, 0, sizeof(unsigned) * db->gse_count * db->module_count);

    // Couneters: total genes from query for each GSE
    unsigned * request = NULL;
    request = (unsigned *) malloc(sizeof(unsigned) * db->gse_count);
    memset(request, 0, sizeof(unsigned) * db->gse_count);

    for (unsigned gse = 0; gse < db->gse_count; gse++) {
        for (unsigned gene = 0; gene < size; gene++) {
            unsigned m = relation[gene * db->gse_count + gse];
            if (m >= db->module_count) continue;
            counters[gse * db->module_count + m] += 1;
            request[gse] += 1;
        }
    }

    // P-value + total
    PyObject * result = Py_BuildValue("[]");
    for (unsigned gse = 0; gse < db->gse_count; gse++) {
        for (unsigned m = 0; m < db->module_count; m++) {
            // Genes in module
            unsigned a = counters[gse * db->module_count + m];
            if (a == 0) continue;

            // Total genes in all modules for this GSE
            unsigned total = decode(&db->counts[gse * (db->module_count * 2 + 1)], 3);

            // Total genes in this module
            unsigned msize = decode(&db->counts[gse * (db->module_count * 2 + 1) + 3 + m * 2], 2);

            int b = msize - a;
            int c = request[gse] - a;
            int d = total - a - b - c;

            double pval = right(a, b, c, d, lf);
            double adj_pval = pval * db->modules_total;

            if (adj_pval > 0.01) continue;

            // printf("%i\t%i\t%i\t%i\t%f\t%f\n", gse, m, a, msize, log10(pval), log10(adj_pval));
            PyObject * line = Py_BuildValue("[i,i,i,i,d,d]", gse, m, a, msize, log10(pval), log10(adj_pval));
            PyList_Append(result, line);
        }
    }

    free(relation);
    free(counters);
    free(request);

    return result;
}
