# Usage:
# python3 pack.names.py [gse names file] [gsm names file] [database dir]
# python3 pack.names.py ../supplement/gse_names.tsv ../supplement/gsm_names.tsv default.gqdb

import sys, json, os, glob

if len(sys.argv) != 4 or not os.path.isfile(sys.argv[1]):
    print("Usage:    python3 %s [gse names file] [gsm names file] [database dir]" % sys.argv[0])
    print("Example:  python3 %s ../supplement/gse_names.tsv ../supplement/gsm_names.tsv default.gqdb" % sys.argv[0])
    sys.exit(1)

gsenames = sys.argv[1]
gsmnames = sys.argv[2]
database = sys.argv[3]


def content(file):
    with open(file, 'r') as handle:
        for line in handle:
            if line == "": continue
            yield line.replace('\n', '').split('\t')[0:2]

def filler(obj, file, export):
    for code, name in content(file):
        obj[code] = name

    titled = {name: obj[name] for name in obj if obj[name] != ''}
    untitled = [name for name in obj if obj[name] == '']

    print("\nFiller: %s" % file)
    print("Titles saved:     %i" % len(titled))
    print("Titles not found: %i" % len(untitled))
    if len(untitled) > 0: print(' '.join(untitled))

    H = open("%s/%s" % (database, export), 'w')
    H.write(json.dumps(titled))
    H.close()

# --------------------------------------------------------------------------- #

gsm = {}
for tsv in glob.glob("%s/*/eigens/*.tsv" % database):
    with open(tsv, 'r') as e:
        for code in next(e).replace('\n', '').split('\t')[1:]:
            gsm[code] = 'NONE'

filler(gsm, gsmnames, 'gsm_names.json')

# --------------------------------------------------------------------------- #

gse = {}
for js in glob.glob("%s/*/offsets.json" % database):
    with open(js, 'r') as e:
        offsets = json.load(e)
        for name in offsets['gse']: gse[name] = ''

filler(gse, gsenames, 'gse_names.json')

# --------------------------------------------------------------------------- #
