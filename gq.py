import os, sys, json, math
# import tracemalloc
# tracemalloc.start()

if len(sys.argv) < 3:
    print("Usage:\n  python3 gq.py [species] [genes list]")
    print("Example:\n  python3 gq.py mm \"15368 93692 104263 99929 15417 16770 15275 12576 66681\"\n")
    sys.exit(0)

f_data = './data/' + sys.argv[1] + '.bin'
f_meta = './data/' + sys.argv[1] + '.m.json'

for f in [f_data, f_meta]:
    if not os.path.isfile(f):
        print("File not found: " + f)
        sys.exit(1)

# --------------------------------------------------------------------------- #

genes = sys.argv[2].split(' ')
data = open(f_data, 'rb')
with open(f_meta, 'r') as e:
    meta = json.load(e)

def num(bin):
    return sum([256**(len(bin) - 1 - i) * v for i,v in enumerate(bin)])

glen = num(data.read(4)) + 13
mmax = num(data.read(1))
gmsx = num(data.read(4))
mcnt = num(data.read(4))

print([glen, mmax, gmsx, mcnt])

def reader():
    vec = []
    while True:
        point = ord(data.read(1))
        if point == 254:
            space = ord(data.read(1)) + 3
            vec.extend([255] * space)
        if point == 253:
            space = ord(data.read(1)) * 255 + ord(data.read(1)) + 3
            vec.extend([255] * space)
        if point <= 252 or point == 255:
            vec.append(point)
        if len(vec) == gmsx: return vec

tree = {}
genes_in_gse = {}
for gene in genes:
    data.seek(meta['genes'][gene] + 13, 0)
    modules = reader()

    for i, gse in enumerate(meta['gse']):
        if modules[i] == 255: continue
        pos = modules[i] + 1

        if gse not in tree:
            tree[gse] = {}
            genes_in_gse[gse] = 0

        genes_in_gse[gse] += 1
        K = '#' + str(modules[i])
        if K not in tree[gse]:
            data.seek(glen + (mmax * 2 + 1) * i + 3 + modules[i] * 2, 0)
            tree[gse][K] = [num(data.read(2)), 0]
        tree[gse][K][1] += 1

# current, peak = tracemalloc.get_traced_memory()
# print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")

# --------------------------------------------------------------------------- #
total = 100000
f = [0 for i in range(0, total)] # log factorials
for i in range(1, total):
    f[i] = f[i - 1] + math.log(i)

def pv(a, b, c, d):
    return math.exp(f[a + b] + f[c + d] + f[a + c] + f[b + d] - f[a + b + c + d] - f[a] - f[b] - f[c] - f[d])

def rightTail(a, b, c, d, psum = 0.0):
    while True:
        psum += pv(a, b, c, d)
        [a, b, c, d] = [a + 1, b - 1, c - 1, d + 1]
        if c <= 0 or b <= 0: break
    return psum

finded = []
for gse in tree:
    request = genes_in_gse[gse]
    total = 6000

    for m in tree[gse]:
        module_count, intersection = tree[gse][m]
        B, C = [module_count - intersection, request - intersection]
        pval = rightTail(intersection, B, C, total - intersection - B - C, 0)
        adj_pval = pval * mcnt
        # print([[intersection, B, C, total - intersection - B - C], adj_pval, pval])
        if adj_pval > 0.01: continue
        finded.append([gse, m, intersection, module_count, math.log10(pval), math.log10(adj_pval)])
        # print([gse, m, intersection, module_count, request, total, math.log10(rtail)])

finded.sort(key = lambda x: x[4])
for result in finded:
    print(result)

# current, peak = tracemalloc.get_traced_memory()
# print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
# tracemalloc.stop()

# python3 gq.py [species] [genes list]
#
