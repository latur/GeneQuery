import genequery
import os, sys, json, math, glob

# --------------------------------------------------------------------------- #
def log(val, color = '37;3'):
    print("\033[" + color + "m" + str(val) + "\033[0m")

args = {'db': './default.gqdb', 'p': '9225', 'h': False}
for k, v in enumerate(sys.argv):
    if v[0] == '-' and v[1:] in args:
        args[v[1:]] = sys.argv[k+1] if (k+1) in sys.argv else ""

def js(name):
    f = "%s/%s.json" % (args['db'], name)
    if os.path.isfile(f):
        with open(f, 'r') as e:
            return json.load(e)
    return {}

# --------------------------------------------------------------------------- #
if args['h'] != False:
    log('Usage:', '37;1')
    log('  python3 gqserver.py [-p <port>] [-db <path>]')
    log('\nArguments:', '37;1')
    log('  -p  <port>  Server port.           Default: 9225')
    log('  -db <path>  Qustom database path.  Default: default.gqdb')
    log('\nExamples:', '37;1')
    log('  python3 gqserver.py')
    log('  python3 gqserver.py -db custom.gqdb -p 80')
    log('  python3 gqserver.py -h')
    log('')
    sys.exit()

# --------------------------------------------------------------------------- #
species = {}
for bin in glob.glob("%s/*/gmt.bin" % args['db']):
    spc = os.path.basename(os.path.dirname(bin))
    with open(bin.replace('gmt.bin', 'offsets.json'), 'r') as e:
        species[spc] = {'gmt': bin, 'meta': json.load(e), 'heatmap': {} }
        for hm in glob.glob("%s/*/eigens/*.tsv" % args['db']):
            species[spc]['heatmap'][os.path.basename(hm).split('.')[0]] = hm

xref = js('xref')
gse_names = js('gse_names')
gsm_names = js('gsm_names')


def offsets(obj, db): # Entrez -> Offsets in .bin file
    data = {}
    for name in obj:
        for id in obj[name]:
            if id not in species[db]['meta']['genes']: continue
            data[id] = species[db]['meta']['genes'][id]
    return data

def convertor(names, query, db):
    names = [name if name[0:3] != 'ENS' else name.split('.')[0] for name in names]
    entrez = {name: [name] for name in names}

    if 'xrefs' in xref:
        for name in entrez:
            if name in xref['xrefs'][query]:
                entrez[name] = xref['xrefs'][query][name]

    data = offsets(entrez, query)

    orthology = {}
    if 'links' in xref and query != db:
        for name in entrez:
            for ID in entrez[name]:
                orthology[ID] = []
                if ID not in xref['links'][query]: continue
                orthology[ID] = xref['groups'][ xref['links'][query][ID] ][db]
        data = offsets(orthology, db)

    return {'orthology': orthology, 'entrez': entrez, 'offsets': data}


# --------------------------------------------------------------------------- #
from flask import Flask, jsonify, request, render_template
app = Flask(__name__, static_url_path='', static_folder='web')
loaded = {}

@app.route('/app/<query>:<db>', methods=['POST'])
def query(query, db):
    if db not in species or query not in species:
        return jsonify({'Error': 'DB not found'})

    if db not in loaded:
        loaded[db] = genequery.load(species[db]['gmt'])

    genes = convertor(request.get_json(), query, db)
    req = [genes['offsets'][k] for k in genes['offsets']]
    gse = genequery.run(req, index=loaded[db])
    gse.sort(key = lambda x: x[4])

    for item in gse:
        item[0] = species[db]['meta']['gse'][item[0]]
        item.append(item[0] if item[0] not in gse_names else gse_names[item[0]])

    return jsonify({'gse': gse, 'genes': genes})


@app.route('/heatmap/<db>/<name>', methods=['POST'])
def heatmap(db, name):
    if db not in species:
        return jsonify({'Error': 'DB not found'})

    if name not in species[db]['heatmap']:
        return jsonify({'Error': 'Heatmap not found'})

    with open(species[db]['heatmap'][name], 'r') as e:
        matrix = [line.replace('\n', '').split('\t') for line in e]
        matrix[0][1:] = [id if id not in gsm_names else gsm_names[id] for id in matrix[0][1:]]
        return jsonify(matrix)


@app.route('/genes/<db>/<name>', methods=['POST'])
def overlap(db, name):
    if db not in species:
        return jsonify({'Error': 'DB not found'})

    if db not in loaded:
        loaded[db] = genequery.load(species[db]['gmt'])

    genes = species[db]['meta']['genes']
    rev = {genes[k]:k for k in genes}
    all = [genes[k] for k in genes]
    gse = species[db]['meta']['gse'].index(name)
    filtred = genequery.filter(all, loaded[db], gse=gse)

    modules = {}
    for m, offset in filtred:
        if m not in modules: modules[m] = []
        modules[m].append(rev[offset])

    return jsonify(modules)


@app.route('/')
def root():
    return render_template('index.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/<path>/<any>')
def all(path, any):
    return root()


if __name__ == "__main__":
    app.run(debug=False, host='0.0.0.0', port=int(args['p']))


# pip3 install flask
# python3 gqserver.py
