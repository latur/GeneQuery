import genequery
import os, sys, json, math, glob

# --------------------------------------------------------------------------- #
def log(val, color = '37;3'):
    print("\033[" + color + "m" + str(val) + "\033[0m")

args = {'db': './default.gqdb', 's': 'hs:hs', 'g': ''}
for k, v in enumerate(sys.argv):
    if v[0] == '-' and v[1:] in args: args[v[1:]] = sys.argv[k+1]
genes = args['g'].replace(' ', ',').split(',')

# --------------------------------------------------------------------------- #
if len(genes) == 1:
    log('Usage:', '37;1')
    log('  python3 gq.py -g <genes> -s <sp>:<sp> [-db <path>]')
    log('\nArguments:', '37;1')
    log('  -g  <genes>    Gene names, separated by space or comma')
    log('  -s  <sp>:<sp>  Query species shortcode : Database species shortcode')
    log('  -db <path>     Qustom database path.  Default: default.gqdb')
    log('\nExamples:', '37;1')
    log('  python3 gq.py -g "15368 93692 ... 104263 99929"')
    log('  python3 gq.py -g 15368,93692,...,104263,99929')
    log('  python3 gq.py -s mm:rn -g 15368,93692,...,104263,99929')
    log('  python3 gq.py -s gp:gp -db custom.gqdb -g 93692,15368,...,104263')
    log('')
    sys.exit()

# --------------------------------------------------------------------------- #
species = {}
for bin in glob.glob("%s/*/gmt.bin" % args['db']):
    spc = os.path.basename(os.path.dirname(bin))
    with open(bin.replace('gmt.bin', 'offsets.json'), 'r') as e:
        species[spc] = {'gmt': bin, 'meta': json.load(e) }

xref = {}
with open("%s/xref.json" % args['db'], 'r') as e:
    xref = json.load(e)

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
        group = xref['groups'][ xref['links'][code_from][id] ]
        if code_to not in group: continue
        request_target.extend(group[code_to])
    return request_target

# --------------------------------------------------------------------------- #
code_from, code_to = args['s'].split(':')
genes = entrez(genes, code_from)
if code_from != code_to:
    genes = use_orthology(genes, code_from, code_to)

offsets = [] # All entrz -> offsets
for id in genes:
    if id not in species[code_to]['meta']['genes']: continue
    offsets.append(species[code_to]['meta']['genes'][id])

gse = genequery.run(offsets, dbname=species[code_to]['gmt'])
gse.sort(key = lambda x: x[4])

if len(gse) == 0:
    log('Not found', '31;1')
    sys.exit(1)

for item in gse:
    item[0] = species[code_to]['meta']['gse'][item[0]]
    print('\t'.join([str(e) for e in item]))
