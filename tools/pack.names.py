# Usage:
# python3 pack.names.py [names file] [database dir]
# python3 pack.xref.py default.gqdb ../supplement/gse_names.tsv

import sys, json, os, glob

if len(sys.argv) != 3 or not os.path.isfile(sys.argv[1]):
    print("Usage:    python3 %s [names file] [database dir]" % sys.argv[0])
    print("Example:  python3 %s ../supplement/gse_names.tsv default.gqdb" % sys.argv[0])
    sys.exit(1)

gsenames = sys.argv[1]
database = sys.argv[2]

gse = {}
for js in glob.glob("%s/*/offsets.json" % database):
    with open(js, 'r') as e:
        offsets = json.load(e)
        for name in offsets['gse']: gse[name] = ''

with open(gsenames, 'r') as handle:
    for line in handle:
        if line == "": continue
        name, title = line.replace('\n', '').split('\t')[0:2]
        if name in gse: gse[name] = title

titles = {name: gse[name] for name in gse if gse[name] != ''}
untitled = [name for name in gse if gse[name] == '']

H = open(database + '/titles.json', 'w')
H.write(json.dumps(titles))
H.close()

if len(untitled) > 0: print(' '.join(untitled))
print("Titles not found: %i" % len(untitled))
print("Titles saved:     %i" % len(titles))
