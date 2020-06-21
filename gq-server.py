import tracemalloc
tracemalloc.start()

import os, sys, json, glob
import genequery
from flask import Flask, jsonify, request

app = Flask(__name__, static_url_path='', static_folder='web')

# --------------------------------------------------------------------------- #

with open('./data/xref.json', 'r') as e:
    xref = json.load(e)

species = {}
for bin in glob.glob('./data/*.bin'):
    code = os.path.basename(bin).replace('.bin', '')
    species[code] = {'id': genequery.load(bin), 'info': {}}
    with open(bin.replace('.bin', '.m.json'), 'r') as e:
        species[code]['info'] = json.load(e)

# --------------------------------------------------------------------------- #
def query(s_from, s_to, genes):
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

    gse = genequery.run(offsets, index=species[s_to]['id'])
    gse.sort(key = lambda x: x[4])
    for item in gse:
        item[0] = species[s_to]['info']['gse'][item[0]]

    return gse

# --------------------------------------------------------------------------- #

@app.route('/run/<s_from>:<s_to>', methods=['POST'])
def execute(s_from, s_to):
    if 'genes' not in request.form: return ""
    gse = query(s_from, s_to, request.form['genes'].split(','))
    return jsonify({ 'gse': gse })


@app.route('/')
def root():
    return app.send_static_file('index.html')

@app.route('/<path>/<any>')
def all(path, any):
    return root()


if __name__ == "__main__":
    app.run(host='0.0.0.0', port=9110)

# pip3 install flask
# python3 gq-server.py
