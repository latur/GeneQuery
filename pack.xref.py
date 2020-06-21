import sys, json

base = sys.argv[1]
files = ['ensembl-to-entrez.tsv', 'refseq-to-entrez.tsv', 'symbol-to-entrez.tsv']

names = {}
for fn in files:
    handle = open(base + fn, 'r')
    for row in handle.read().split('\n'):
        if row == "": continue
        species, entrez, ref = row.split('\t')
        if species not in names: names[species] = {}
        if ref not in names[species]: names[species][ref] = []
        names[species][ref].append(entrez)

groups = {}
handle = open(base + 'orthology.tsv', 'r')
for row in handle.read().split('\n'):
    if row == "": continue
    group, species, entrez = row.split('\t')[0:3]
    if group not in groups: groups[group] = {}
    if species not in groups[group]: groups[group][species] = []
    groups[group][species].append(entrez)

refs = {}
for g in groups:
    if len(groups[g]) == 1: continue
    for sp in groups[g]:
        if sp not in refs: refs[sp] = {}
        for gene in groups[g][sp]:
            refs[sp][gene] = g

handle = open('./data/xref.json', 'w')
handle.write(json.dumps({'names': names, 'groups': groups, 'refs': refs}))
handle.close()

# python3 pack.xref.py ../gqdb/

# Usage:
# xref['groups'][xref['refs'][from_sp][entrz_id]] -> {'hh': ..., 'mm': ...}
