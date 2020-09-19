let qs = (e) => document.querySelector(e);

let palette = [
  [49,   54, 149],
  [69,  117, 180],
  [116, 173, 209],
  [171, 217, 233],
  [224, 243, 248],
  [254, 224, 144],
  [253, 174, 97 ],
  [244, 109, 67 ],
  [215,  48, 39 ],
  [165,   0, 38 ]
];

function scale(num, a, b, x, y) {
  return x + (y - x) * (num - a)/(b - a);
}

function color(num, a, b) {
  let rgb = '0,0,1';
  let step = (b - a) / (palette.length - 1);

  if (num < a) num = a;
  if (num > b) num = b;

  for (let i = 0; i < palette.length - 1; i++) {
    const l = a + i * step;
    const r = a + (i + 1) * step;
    if (num < l || num > r) continue;

    rgb = [0,1,2].map((c) => {
      return scale(num, l, r, palette[i][c], palette[i + 1][c]);
    }).join(',');
  }

  return 'rgb(' + rgb + ')';
}

function Request(url, data, after) {
  let xhr = new XMLHttpRequest();
  xhr.open('POST', url, true);
  xhr.setRequestHeader("Content-type", "application/json");
  xhr.onload = function () {
    after(JSON.parse(this.responseText));
  };
  xhr.send(JSON.stringify(data));
}

function Heatmap(tsv, current) {
  return '<div class="tsv">' + tsv.map((col, i) => {
    return '<div class="tsv-col c'+i+' '+(current == (i+1) ? 'current' : '')+'">'+ col.map((v, k) => {
      if (i == 0 || k == 0) return ('<div class="name row">'+v+'</div>');
      let disp = (parseFloat(v)).toFixed(2).replace('0.', '.');
      return ('<div class="value row" style="background: '+color(parseFloat(v), -0.6, 0.6)+'" data-v="'+v+'">'+disp+'</div>');
    }).join('') + '</div>';
  }).join('') + '</div>';
}

function Modal(content) {
  qs('[data-modals]').innerHTML = Template('modal', {'content': content});
  qs('body').classList.add('modal');
  setTimeout(function () {
    qs('body').classList.add('modal-animation');
  }, 10);
}

function ModalClose() {
  qs('body').classList.remove('modal');
  qs('body').classList.remove('modal-animation');
  let m = qs('[data-modals] .modal-wrapper');
  m.parentElement.removeChild(m);
}

function ShowHeatmap(url, m) {
  Request('/heatmap/' + url, {}, function(res){
    Modal(Heatmap(res, m));
    qs('[data-modal-window]').style.width = res.length * 22 + 80 + 40 + 'px';
  });
}

function ShowOverlap(db, name, module) {
  Request('/genes/' + [db, name].join('/'), {}, function(res){
    if (!res[module]) return ;
    let common = {};
    res[module].map((id) => { if (window.last_request[id]) common[id] = true; });

    Modal(Template('overlap', {
      'name': name,
      'module': module,
      'common': Object.keys(common).join(' '),
      'full': res[module].filter((id) => !common[id]).join(' ')
    }));
  });
}

function DownloadGMT(db, name) {
  Request('/genes/' + [db, name].join('/'), {}, function(res){
    let tsv = '';
    for (let module in res) tsv += module + '\t' + res[module].join(' ') + '\n';

    let e = document.createElement('a');
    e.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(tsv));
    e.setAttribute('download', name + '.gmt.tsv');

    e.style.display = 'none';
    document.body.appendChild(e);
    e.click();
    document.body.removeChild(e);
  });
}

function FullGenesList() {
  qs('[data-genes-lists]').classList.toggle('show-all');
}

qs('[data-query]').addEventListener('click', function(){
  let req = [ qs('[name="query"]:checked').value, qs('[name="db"]:checked').value ];
  let genes = qs('[name="genes"]').value;
  genes = genes.replace(/,|\t|\n/g, ' ').replace(/\s\s+/g, ' ').split(' ');
  genes = genes.filter(function(e){ return e != ""; });

  Request('/app/' + req.join(':'), genes, function(res){
    window.last_request = res.genes.offsets;
    let gse = res.gse.map(function(row){
      return {
        'code': row[0],
        'name': row[6],
        'module': row[1],
        'genes': row[2],
        'total': row[3],
        'GSE': row[0].split('_')[0],
        'log10_apv': row[5].toFixed(3),
        'db': req[1]
      };
    });
    console.log(gse);
    console.log(res);
    if (gse.length > 0) {
      qs('[data-results]').innerHTML = Template('table', {'data': gse});
    } else {
      qs('[data-results]').innerHTML = Template('empty-result', {});
    }
  });
}, false);

document.addEventListener('keyup', function(e){
  if (e.keyCode === 27) ModalClose();
});
