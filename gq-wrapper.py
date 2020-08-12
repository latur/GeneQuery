import os, sys, json, math, glob
import genequery

if len(sys.argv) < 4:
    print("Usage:\n  python3 "+sys.argv[0]+" [species code from] [species code to] [genes list,..]")
    print("Example:\n  python3 "+sys.argv[0]+" mm mm 15368 93692 104263 99929 15417 16770 15275 12576 66681\n")
    sys.exit(0)

# --------------------------------------------------------------------------- #
xref = {}
with open('./data/xref.json', 'r') as e:
    xref = json.load(e)

species = {}
for bin in glob.glob('./data/*.gmt.bin'):
    code = os.path.basename(bin).replace('.gmt.bin', '')
    with open(bin.replace('.gmt.bin', '.offsets.json'), 'r') as e:
        species[code] = {'gmt': bin, 'meta': json.load(e) }

# All names -> entrez
def entrez(names, code):
    request = []
    for name in names:
        if name in xref['xrefs'][code]:
            request.extend(xref['xrefs'][code][name])
        else:
            request.append(name)
    return request

def use_orthology(request, code_from, code_to):
    request_target = []
    for id in request:
        if id not in xref['links'][code_from]: continue
        group = xref['groups'][xref['links'][code_from]][id]
        if code_to not in group: continue
        request_target.extend(group[code_to])
    return request_target

# --------------------------------------------------------------------------- #
code_from, code_to = sys.argv[1:3]
genes = sys.argv[3:]
genes = entrez(genes, code_from)
if code_from != code_to:
    genes = use_orthology(genes, code_from, code_to)

# All entrz -> offsets
offsets = []
for id in genes:
    if id not in species[code_to]['meta']['genes']: continue
    offsets.append(species[code_to]['meta']['genes'][id])

print(genes)
print(offsets)

gse = genequery.run(offsets, dbname=species[code_to]['gmt'])
gse.sort(key = lambda x: x[4])
for item in gse:
    item[0] = species[code_to]['meta']['gse'][item[0]]
    print(item)
