# Usage:
# python3 pack.xref.py [files dir] [species names] > [output file]
# python3 pack.xref.py ../supplement hs mm rn > default.gqdb/xref.json

import sys, json, os

supplement = sys.argv[1]
samples = sys.argv[2:]

def content(file):
    with open("%s/%s" % (supplement, file), 'r') as handle:
        for line in handle:
            if line == "": continue
            yield line.replace('\n', '').split('\t')

def reader(sample, name):
    return content("%s/%s.%s.tsv" % (sample, sample, name))

xrefs = {}

for sample in next(os.walk(supplement))[1]:

    xref = {}

    for ensg, name, entrez in reader(sample, 'fix_gene3'):
        if entrez == 'NONE': continue
        ensg = ensg.split('.')[0]
        if ensg not in xref: xref[ensg] = {}
        xref[ensg][entrez] = 1
        if name not in xref: xref[name] = {}
        xref[name][entrez] = 1

    for enst, ensg in reader(sample, 'enst2ensg'):
        if ensg in xref: xref[enst] = xref[ensg]

    for refseq, name, entrez in reader(sample, 'refseq2gene'):
        if entrez == 'NONE': continue
        if refseq not in xref: xref[refseq] = {}
        xref[refseq][entrez] = 1
        if name not in xref: xref[name] = {}
        xref[name][entrez] = 1

    for name, entrez in reader(sample, 'synonyms'):
        if name not in xref: xref[name] = {}
        xref[name][entrez] = 1

    for name in xref:
        xref[name] = [k for k in xref[name]]

    xrefs[sample] = xref

links = {}
groups = []
for values in content("orthologs.tsv"):
    val = dict(zip(samples, [e.split(',') for e in values[0:len(samples)]]))
    group = len(groups)
    groups.append(val)
    for sample in val:
        if sample not in links: links[sample] = {}
        for ID in val[sample]:
            if ID != 'NONE':
                links[sample][ID] = group

print(json.dumps({'xrefs': xrefs, 'links': links, 'groups': groups}))

# Call:
# xref['groups'][ xref['links'][from_sample][entrz_ID] ] -> {'hh': ..., 'mm': ...}
