typedef struct {
    unsigned module_count;
    unsigned gse_count;
    unsigned modules_total;
    unsigned char * genes;
    unsigned char * counts;
} gmt;

gmt * dbset[99];
unsigned dbi = 0;

unsigned block(FILE * db, unsigned length)
{
    unsigned char byte;
    unsigned result = 0;
    unsigned powers[] = { 256*256*256, 256*256, 256, 1 };
    for (unsigned n = 0; n < length; n++) {
        fread(&byte, 1, 1, db);
        result += byte * powers[n + 4 - length];
    }
    return result;
}

unsigned decode(unsigned char * start, unsigned length)
{
    unsigned result = 0;
    unsigned powers[] = { 256*256*256, 256*256, 256, 1 };
    for (unsigned n = 0; n < length; n++) {
        result += start[n] * powers[n + 4 - length];
    }
    return result;
}

gmt * reader(char * dbname)
{
    FILE *f;
    if ((f = fopen(dbname, "rb")) == NULL) return NULL;

    gmt * db = (gmt *) malloc(sizeof(gmt));

    // Header [4] [1] [4] [4]
    unsigned genes = sizeof(unsigned char) * block(f, 4);
    db->module_count = block(f, 1);
    db->gse_count = block(f, 4);
    db->modules_total = block(f, 4);

    // Genes to modules matrix
    db->genes = (unsigned char *) malloc(genes);
    fread(db->genes, 1, genes, f);

    // GSE modules counts
    unsigned counts = sizeof(unsigned char) * (db->module_count * 2 + 1) * db->gse_count;
    db->counts = (unsigned char *) malloc(counts);
    fread(db->counts, 1, counts, f);

    fclose(f);
    return db;
}

double * lf;

void loadLogFactorials(char * filename, unsigned universe)
{
    FILE * fp;
    if ((fp = fopen(filename, "r")) == NULL) {
        printf("Unable to read Log[Factorial]'s file!\n\n");
        exit(1);
    }

    lf = (double *) malloc(sizeof(double) * universe);
    char chunk[64];
    char * end;

    unsigned i = 0;
    while (fgets(chunk, sizeof(chunk), fp) != NULL && i != universe) {
        lf[i] = strtod (chunk, &end);
        i += 1;
    }

    fclose(fp);
}
