import os, sys, json

genes = {}
gseModules = {}
maxModules = 1

gmt_file = sys.argv[1]
export_name = sys.argv[2]

gmtH = open(gmt_file, 'r')
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

# Sort GSE modules for better compression
gseStrings = []
for gse in gseModules:
    chrs = [chr(255 if gse not in genes[gene] else int(genes[gene][gse])) for gene in genes]
    gseStrings.append([''.join(chrs), gse])
gseModulesList = [gse for s, gse in sorted(gseStrings, key=lambda x:x[0])]

index = {}
dbset = []
for gene in genes:
    index[gene] = len(dbset)
    vec = [(255 if gse not in genes[gene] else int(genes[gene][gse])) for gse in gseModulesList]
    dbset.extend([i for i in flatten(compressed(vec))])

# Header:  [modules offset]x4  [max modules]  [gse count]*4  [modules total]x4
header = []
header.extend(bytes(len(dbset), 4))
header.append(maxModules)
header.extend(bytes(len(gseModules), 4))
header.extend(bytes(sum([(len(gseModules[gse]) - 1) for gse in gseModulesList]), 4))

# --------------------------------------------------------------------------- #
# GSE modules: Countrs matrix

counters = []
for gse in gseModulesList:
    line = bytes(gseModules[gse][0], 3)
    for counter in gseModules[gse][1:]:
        line.extend(bytes(counter))
    line.extend((maxModules * 2 + 1 - len(line)) * [0])
    dbset.extend(line)

header.extend(dbset)

# ---------------------------------------------------------------------------- #
H = open("./data/%s.gmt.bin" % (export_name, ), 'wb')
H.write(bytearray(header))
H.close()

H = open("./data/%s.offsets.json" % (export_name, ), 'w')
H.write(json.dumps({'genes': index, 'gse': [gse for gse in gseModulesList]}))
H.close()
