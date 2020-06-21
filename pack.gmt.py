import os, sys, json

genes = {}
gseModules = {}
maxModules = 1

gmtH = open(sys.argv[1], 'r')
for row in gmtH.read().split('\n'):
    if row == '': continue
    obj = row.replace(',', '\t').replace('#', '\t').split('\t')
    gse, module = obj[0:2] # [GSE_GPL, 1]
    size = len(obj[2:])    # [Genes list]

    if gse not in gseModules:
        if module != '0': continue
        gseModules[gse] = [0]

    gseModules[gse][0] += size
    gseModules[gse].append(size)

    if maxModules < len(gseModules[gse]):
        maxModules = len(gseModules[gse])

    for gene in obj[2:]:
        if gene not in genes: genes[gene] = {}
        if gse in genes[gene]:
            raise Exception('[!] One gene in multilpe modules for GSE set')
        genes[gene][gse] = module
gmtH.close()

# --------------------------------------------------------------------------- #
def reset(cache):
    if cache <= 2: return [255] * cache
    val = cache - 3
    if val <= 255: return [254, val]
    div = (cache - 3) % 255
    return [253, int((val - div)/255), div]

def compressed(vec):
    i = 0; cache = 0
    while i < len(vec):
        if vec[i] == 255: i += 1; cache += 1; continue
        if cache > 0: yield reset(cache); cache = 0
        yield vec[i]
        i += 1
    if cache > 0:
        yield reset(cache)

def flatten(tree):
    for e in tree:
        if isinstance(e, list):
            for i in flatten(e):
                yield i
        else:
            yield e

def chunks(items, n):
    for i in range(0, len(items), n):
        yield items[i:i + n]

def bytes(num, size = 2):
    if num > (256**size - 1):
        raise Exception('[!] Base convert error! Out of memory')
    hexv = "{0:x}".format(num)
    hexv = (size * 2 - len(hexv)) * '0' + hexv
    return [int(e, base=16) for e in chunks(hexv, 2)]

index = {}
dbset = []
for gene in genes:
    index[gene] = len(dbset)
    vec = [(255 if gse not in genes[gene] else int(genes[gene][gse])) for gse in gseModules]
    dbset.extend([i for i in flatten(compressed(vec))])

# Header:  [modules offset]x4  [max modules]  [gse count]*4  [modules total]x4
header = []
header.extend(bytes(len(dbset), 4))
header.append(maxModules)
header.extend(bytes(len(gseModules), 4))
header.extend(bytes(sum([(len(gseModules[gse]) - 1) for gse in gseModules]), 4))

# --------------------------------------------------------------------------- #
# GSE modules: Countrs matrix

counters = []
for gse in gseModules:
    line = bytes(gseModules[gse][0], 3)
    for counter in gseModules[gse][1:]:
        line.extend(bytes(counter))
    line.extend((maxModules * 2 + 1 - len(line)) * [0])
    dbset.extend(line)

header.extend(dbset)

# ---------------------------------------------------------------------------- #
prefix = sys.argv[2]

tcgmtH = open(prefix + '.bin', 'wb')
tcgmtH.write(bytearray(header))
tcgmtH.close()

metaH = open(prefix + '.m.json', 'w')
metaH.write(json.dumps({
    'genes': index,
    'gse': [gse for gse in gseModules]
}))
metaH.close()

# python3 pack.gmt.py ../gqdb/hs.modules.gmt data/hs
# python3 pack.gmt.py ../gqdb/mm.modules.gmt data/mm
# python3 pack.gmt.py ../gqdb/rt.modules.gmt data/rt
