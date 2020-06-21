import os, sys, json, math, glob
import genequery

import tracemalloc
tracemalloc.start()

if len(sys.argv) < 4:
    print("Usage:\n  python3 "+sys.argv[0]+" [species from] [species to] [genes list,..]")
    print("Example:\n  python3 "+sys.argv[0]+" mm mm 15368 93692 104263 99929 15417 16770 15275 12576 66681\n")
    sys.exit(0)

with open('./data/xref.json', 'r') as e:
    xref = json.load(e)

species = {}
for bin in glob.glob('./data/*.bin'):
    code = os.path.basename(bin).replace('.bin', '')
    species[code] = {'id': genequery.load(bin), 'info': {}}
    with open(bin.replace('.bin', '.m.json'), 'r') as e:
        species[code]['info'] = json.load(e)

# --------------------------------------------------------------------------- #

s_from, s_to = sys.argv[1:3]
genes = sys.argv[3:]

# All names -> entrez
request = []
for name in genes:
    if name in xref['names'][s_from]:
        request.extend(xref['names'][s_from][name])
    else:
        request.append(name)

# Use orthology
if s_from != s_to:
    request_target = []
    for id in request:
        if id not in xref['refs'][s_from]: continue
        group = xref['groups'][xref['refs'][s_from][id]]
        if s_to not in group: continue
        request_target.extend(group[s_to])
    request = request_target

# All entrz -> offsets
offsets = []
for id in request:
    if id not in species[s_to]['info']['genes']: continue
    offsets.append(species[s_to]['info']['genes'][id])

gse = genequery.run(species[s_to]['id'], offsets)
gse.sort(key = lambda x: x[4])
for item in gse:
    item[0] = species[s_to]['info']['gse'][item[0]]
    print(item)

current, peak = tracemalloc.get_traced_memory()
print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
