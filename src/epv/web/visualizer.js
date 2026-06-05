(function (global) {
  'use strict';

  /*!
   * Project: Enrichment Profile Visualizer
   * Author: FEI YUAN <fei.yuan@bcm.edu>
   * License: MIT
   * (c) 2025 FEI YUAN
   */

  Dropzone.options.dropzone = {
    acceptedFiles: ".csv, .tsv, .json, .gz, text/csv, text/tab-separated-values, application/json, application/gzip",
    url: "#",
    autoProcessQueue: false,
    maxFiles: 1,
    uploadMultiple: false,
  };

  const DATA = {
    rows: [],
    columns: [],
    libraries: [],
    countColumns: [],
    scoreColumns: [],
    smilesColumns: [],
    hits: new Map(),
    tops: new Map(),
    cards: {},
    uniques: [],
    duplicates: {},
    x: null,
    y: null,
    vs: null,
    library: null,
    hitsTable: null,
    topHitsTable: null
  }
  let R = {
    root: document,
    els: null,
    io: null,
    config: null,
    utilities: null,
    ...DATA
  };

  const q = id => R.root.getElementById(id);
  const readConfigForm = () => {
    return {
      fontSize: Number(q('fontSize')?.value || 12),

      render: {
        SMILES: !!q('renderSmiles')?.checked,
        c1_smiles: !!q('renderBB1')?.checked,
        c2_smiles: !!q('renderBB2')?.checked,
        c3_smiles: !!q('renderBB3')?.checked,
      },

      structure: {
        width:  Number(q('structWidth')?.value || 250),
        height: Number(q('structHeight')?.value || 100),
      },

      nTopHits: Number(q('nTopHits')?.value || 5),
      zscoreCutoff: Number(q('zscoreCutoff')?.value || 1.0),

      colors: {
        mono: q('colorMono')?.value || '#0000ad',
        di:   q('colorDi')?.value   || '#00ff33',
        tri:  q('colorTri')?.value  || '#ffcc00',
      },
    };
  }

  R.els = {
    btnHitsModal         : q('btnHitsModal'),
    btnTopHitsModal      : q('btnTopHitsModal'),
    btnSaveSession       : q('btnSaveSession'),
    switchers            : q('switchers'),
    selectors            : q('selectors'),
    btnLibrary           : q('btnLibrary'),
    librarySel           : q('librarySel'),
    uploadPanel          : q('uploadPanel'),
    fileInput            : q('fileInput'),
    dz                   : q('dropzone'),
    chartPanel           : q('chartPanel'),
    configModal          : q('configModal'),
    encodingModal        : q('encodingModal'),
    hhModal              : q('hhModal'),
    uploadModal          : q('uploadModal'),
    encodingTable        : q('encodingTable'),
    encodingSMILES       : q('encodingSMILES'),
    hhTable              : q('hhTable'),
    hhSMILES             : q('hhSMILES'),
    topHitsModal         : q('topHitsModal'),
    hitsModal            : q('hitsModal'),
    hitsTable            : q('hitsTable'),
    numHits              : q('numHits'),
    numTopHits           : q('numTopHits'),
    topHitsTable         : q('topHitsTable'),
  }

  R.io = (() => {
    const isGzip = (u8) => u8 && u8.length >= 2 && u8[0] === 0x1f && u8[1] === 0x8b;
    const inferDelimiter = (name = "") => {
      if (/\.tsv(\.gz)?$/i.test(name)) return "\t";
      if (/\.csv(\.gz)?$/i.test(name)) return ",";
      return undefined;
    };
    const papaParse = (input, download=false) => new Promise((resolve, reject) => {
      const delimiter = inferDelimiter(input);
      Papa.parse(input, {
        delimiter,
        download: download,
        header: true,
        dynamicTyping: true,
        skipEmptyLines: true,
        worker: true,
        complete: ({ data }) => {
          if (data && data.length) {
            resolve(data.map((row, i) => {
              return { ...row, index: i }
            }));
          } else {
            reject(new Error("The text/file is empty or contains no valid data."));
          }
        },
        error: reject
      });
    });
    const saveSession = () => {
      const getCircularReplacer = () => {
        const seen = new WeakSet();
        return (key, value) => {
          if (typeof value === "object" && value !== null) {
            if (seen.has(value)) {
              return;
            }
            seen.add(value);
          }
          return value;
        };
      };
      const serializableTops = Object.fromEntries(
        Array.from(R.tops.entries()).map(([vsKey, rowSet]) => [vsKey, Array.from(rowSet)])
      );
      const serializableHits = Object.fromEntries(
          Array.from(R.hits.entries()).map(([vsKey, rowSet]) => [vsKey, Array.from(rowSet)])
      );

      const data = {
        config: R.config,
        rows: R.rows,
        columns: R.columns,
        libraries: R.libraries,
        countColumns: R.countColumns,
        scoreColumns: R.scoreColumns,
        smilesColumns: R.smilesColumns,
        hits: serializableHits,
        tops: serializableTops,
        uniques: R.uniques,
        duplicates: R.duplicates,
        x: R.x,
        y: R.y,
        vs: R.vs,
        library: R.library,
        cards: R.cards
      }
      const jsonString = JSON.stringify(data, getCircularReplacer(), 2);
      const compressedData = pako.gzip(jsonString, { to: 'string' });
      const blob = new Blob([compressedData], { type: 'application/gzip' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = `VEP.session.${R.utilities.getCompactTimestamp()}.json.gz`;
      link.style.display = 'none';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      URL.revokeObjectURL(url);
    }
    const load = async (input) => {
        if (!input) throw new Error("No file or URL provided.");

        let rawInput = input;
        let fileName = (typeof input === 'object' && input.name) ? input.name : (typeof input === 'string' ? input : null);

        if (typeof input === "string") {
            let res = await fetch(input, { mode: "cors" });
            if (!res.ok) throw new Error(`Fetch failed: ${res.status} ${res.statusText}`);

            const u8 = new Uint8Array(await res.arrayBuffer());
            if (isGzip(u8)) {
                if (!window.pako) throw new Error("pako library required for gzip decompression.");
                rawInput = window.pako.ungzip(u8, { to: "string" });
            } else {
                rawInput = new TextDecoder().decode(u8);
            }
            fileName = fileName || input;
        } else if (input instanceof Blob) {
            const u8 = new Uint8Array(await input.arrayBuffer());
            if (isGzip(u8)) {
                if (!window.pako) throw new Error("pako library required for gzip decompression.");
                rawInput = window.pako.ungzip(u8, { to: "string" });
            } else {
                rawInput = input; // Pass the original Blob for text/CSV parsing
            }
        } else {
            throw new Error("loadAndProcessFile expects a URL string or a File/Blob.");
        }

        if (fileName && (fileName.endsWith('.json') || fileName.endsWith('.json.gz'))) {
            let jsonString;
            if (typeof rawInput !== 'string') {
                jsonString = await new Response(rawInput).text();
            } else {
                jsonString = rawInput;
            }

            try {
                return JSON.parse(jsonString);
            } catch (e) {
                throw new Error(`Failed to parse session file: ${e.message}`);
            }
        } else {
          return await papaParse(rawInput, false);
        }
    };
    return { load, saveSession };
  })();

  R.utilities = (() => {
    const findColumns = (columns, prefix = '', suffix = '', case_sensitive = false) => {
      const normPrefix = case_sensitive ? prefix : prefix.toLowerCase();
      const normSuffix = case_sensitive ? suffix : suffix.toLowerCase();

      const cs = [];
      for (const column of columns) {
        const normColumn = case_sensitive ? column : column.toLowerCase();

        const okPrefix = prefix ? normColumn.startsWith(normPrefix) : true;
        const okSuffix = suffix ? normColumn.endsWith(normSuffix) : true;

        if (okPrefix && okSuffix) {
          cs.push(column);
        }
      }

      if (cs.length === 0) {
        throw Error(`No column was found with prefix "${prefix}" and suffix "${suffix}"`);
      }
      return cs;
    };
    const keyForRow = (row) => {
      if (Object.prototype.hasOwnProperty.call(row, 'compound') && row.compound != null && row.compound !== '') {
        return  String(row.compound);
      }
      return [String(row.library ?? ''), ...R.smilesColumns.map(c => String(row?.[c] ?? ''))]
          .join('|').replace(/\s+/g, '')  // Remove all whitespace
          .replace(/[^a-zA-Z0-9_-]/g, '') // Remove all non-safe characters
          .replace(/^[0-9_-]/, '')       // Remove leading numbers and underscores

    };
    const getSMILES = (row) => Object.entries(R.config.render).filter(([, v]) => v).map(([k]) => row?.[k] ?? '');
    const ClientXY = (id) => {
      const gd = R.els.chartPanel;
      const data = gd._fullData || gd.data || [];
      const rect = gd.getBoundingClientRect();
      let clientX = null;
      let clientY = null;

      for (let i = 0; i < data.length; i++) {
        const tr = data[i];
        const ids = tr?.ids;
        if (!ids) continue;

        const j = ids.indexOf(id);
        if (j >= 0) {
          const x = tr.x[j];
          const y = tr.y[j];

          const xa = gd._fullLayout[(tr.xaxis || 'x') + 'axis'];
          const ya = gd._fullLayout[(tr.yaxis || 'y') + 'axis'];

          const xPx = xa._offset + xa.l2p(xa.d2l(x));
          const yPx = ya._offset + ya.l2p(ya.d2l(y));

          clientX = rect.left + xPx;
          clientY = rect.top + yPx;
          break
        }
      }
      return {x1: clientX, y1: clientY};
    };
    const Draggable = (el, handle) => {
      let startX, startY, startLeft, startTop;

      const onDown = (e) => {
        if (e.button !== undefined && e.button !== 0) return;
        handle.setPointerCapture?.(e.pointerId);
        const r = el.getBoundingClientRect();
        startX = e.clientX; startY = e.clientY;
        startLeft = r.left;  startTop  = r.top;
        window.addEventListener('pointermove', onMove);
        window.addEventListener('pointerup', onUp, { once: true });
        e.preventDefault();
      };

      const onMove = (e) => {
        const dx = e.clientX - startX;
        const dy = e.clientY - startY;
        let left = startLeft + dx, top = startTop + dy;
        el.style.left = left + 'px';
        el.style.top  = top  + 'px';
        el.dispatchEvent(new CustomEvent('dragmove', { detail: { left, top } }));
      };

      const onUp = (e) => {
        handle.releasePointerCapture?.(e.pointerId);
        window.removeEventListener('pointermove', onMove);
      };

      handle.addEventListener('pointerdown', onDown);
      handle.style.cursor = 'grab';
      el.style.touchAction = 'none';
    };
    const saveCardPosition = (card) => {
      const {left, top} = card.state;
      const key = card.getAttribute('id');
      R.cards[R.vs] ??= {};
      R.cards[R.vs][key] = {left: left, top: top};
    };
    const updateConnector = (card) => {
      let x1, y1, x2, y2, left, top;

      const styleLeft = parseFloat(card.style.left);
      const styleTop = parseFloat(card.style.top);

      if (Number.isFinite(styleLeft) && Number.isFinite(styleTop)) {
        left = styleLeft;
        top = styleTop;
      } else if (card?.state) {
        ({left, top} = card.state);
      } else {
        const rect = card.getBoundingClientRect();
        left = rect.left;
        top = rect.top;
      }

      const rect = card.getBoundingClientRect();
      ({ x1, y1 } = ClientXY(card.getAttribute('id')));

      x2 = rect.left + rect.width / 2;
      y2 = rect.bottom;

      card.line.setAttribute('x1', x1);
      card.line.setAttribute('y1', y1);
      card.line.setAttribute('x2', x2);
      card.line.setAttribute('y2', y2);

      card.state = {left: left, top: top, x1: x1, y1: y1, x2: x2, y2: y2}
      saveCardPosition(card);
    };
    const removeConnector = (card) => {
      if (!card.line) return;
      card.removeEventListener('dragmove', card.updateFunction);
      card.line.remove();
      card.line = null;
      card.updateFunction = null;
    };
    const attachConnector = (card) => {
        let ov = document.getElementById('plot-overlay');
        if (!ov) {
          ov = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
          ov.id = 'plot-overlay';
          Object.assign(ov.style, {
            position: 'fixed', left: 0, top: 0, width: '100vw', height: '100vh',
            pointerEvents: 'none', zIndex: 1049
          });
          document.body.appendChild(ov);
        }

        const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        line.setAttribute('stroke', '#afafaf');
        line.setAttribute('stroke-width', '1');
        line.setAttribute('stroke-linecap', 'round');
        ov.appendChild(line);

        card.line = line;

        const updateFunction = () => updateConnector(card);
        card.updateFunction = updateFunction;

        requestAnimationFrame(updateFunction);
        card.addEventListener('dragmove', updateFunction);
        return card;
    };
    const showCompound = (id, restyle=true, addToTop=true) => {
      if (restyle) restylePoints(null, [id], 'add');
      const row = R.uniques.filter(r => r.key === id)[0];
      assembleCompoundCard(row);
      if (addToTop) {
        R.tops.get(R.vs).add(row);
        buildTopHitsTable();
      }
    };
    const showCompounds = () => {
      const topsArray = Array.from(R.tops.get(R.vs) || []);
      const hitsArray = Array.from(R.hits.get(R.vs) || []);
      const rows = [...(topsArray ?? []), ...(hitsArray ?? [])].filter(row => row.library === R.library);
      if (rows.length > 0) {
        restylePoints(rows, null, 'replace');
        for (const row of rows) {
          showCompound(row.key, false, false);
        }
      }
    };
    const removeCompounds = () => {
      const cards = document.querySelectorAll('.compound-card');
      if (cards.length > 0) {
        for (const card of cards) {
          removeCompound(card.getAttribute('id'), false)
        }
      }
    };
    const removeCompound = (id, removeFromTop=true) => {
      const card = q(id);
      if (card) {
        removeConnector(card);
        card.remove();
      }
      restylePoints(null, [id], 'remove');

      if (removeFromTop && R.tops.get(R.vs)) {
          const topSet = R.tops.get(R.vs);
          const hitsSet = R.hits.get(R.vs);

          for (const row of topSet) {
              if (row.key === id) {
                  topSet.delete(row);
                  break;
              }
          }
          if (hitsSet) {
            for (const row of hitsSet) {
              if (row.key === id) {
                    topSet.delete(row);
                    break;
                }
            }
          }

          buildTopHitsTable();
          buildHitsTable();
      }
    };
    const restylePoints = (rows=null, ids=null, mode='add') => {
      const uids = ids || rows?.map(row => row.key);
      if (uids) {
        const want = new Set(uids);
        const data = R.els.chartPanel.data;
        const updaters = []
        for (let i = 0; i < data.length; i++) {
          const tr = data[i];
          const ids = tr?.ids;
          if (!ids) continue;

          const matches = [];
          for (let j = 0; j < ids.length; j++) if (want.has(ids[j])) matches.push(j);
          let selected = Array.isArray(tr.selectedpoints) ? tr.selectedpoints.slice() : [];

          if (mode === 'add') {
            selected = Array.from(new Set(selected.concat(matches)));
          } else if (mode === 'remove') {
            const rm = new Set(matches);
            selected = selected.filter(j => !rm.has(j));
          } else {
            selected = matches;
          }

          const update = {
            selectedpoints: [selected],
            "selected.marker.size": 40,
            "unselected.marker.opacity": 0.5
          };
          updaters.push(Plotly.restyle(R.els.chartPanel, update, [i]));
        }
        return Promise.all(updaters);
      }
      return Promise.resolve();
    };
    const assembleCompoundCard = (row, visible=true) => {
      const smiles = getSMILES(row);
      const trs = smiles.map(s => `<tr><td colspan="2"><div data-smiles="${s}" data-height="125"></div></td></tr>`);
      const text = assembleHoverText(row);
      const parts = text.split('<br>');
      const title = row.compound
      for (let i=1; i < (parts.length-2); i++) {
        let [k, v] = parts[i].split(': ');
        trs.push(`<tr><td class="text-start">${k}</td><td class="text-end">${v}</td></tr>`);
      }
      const hhs = row.nHH > 0 ? `<div data-action="hh" data-key="${row.key}" role="button" tabindex="0" class="text-primary">${row.nHH}</div>` : `${row.nHH}`;
      trs.push(`<tr><td class="text-start">nHH</td><td class="text-end">${hhs}</td></tr>`);

      const key = row.key;
      let card = document.createElement('div');
      const width = R.config.structure.width <= 240 ? 250 : R.config.structure.width + 10;
      card.setAttribute('id', key)
      card.className = 'card compound-card';
      card.style.position = 'fixed';
      card.style.width = `${width}px`;
      card.style.zIndex = '1050';
      card.style.pointerEvents = 'auto';

      const header = document.createElement('div');
      header.className = 'card-header d-flex justify-content-center align-items-center fw-bold gap-2';
      header.innerHTML = title;
      card.appendChild(header);

      const body = document.createElement('div')
      body.innerHTML = `<table class="w-100">${trs.join('')}</table>`;
      body.style.padding = '5px';
      card.appendChild(body);

      const footer = document.createElement('div');
      footer.className = 'card-footer bg-white d-flex align-items-center p-1';
      footer.innerHTML = `<i class="bi bi-bag me-3" data-action="bag" data-key=${key} role="button" tabindex="0" title="Add to candidate"></i> ` +
        `<i class="bi bi-copy text-success me-3" data-action="copy" data-key=${key} role="button" tabindex="0" title="Copy"></i> ` +
        `<i class="bi bi-envelope text-primary me-3" data-action="email" data-key=${key} role="button" tabindex="0" title="E-mail"></i> `;

      if (row.copies > 1) {
        footer.innerHTML += `<span class="text-success" data-action="encoding" data-key=${key} role="button" tabIndex="0" title="Encodings">E: ${row.copies}</span> `;
      }

      footer.innerHTML += `<i class="bi bi-x-circle ms-auto text-danger" data-action="close" data-key=${key} role="button" tabindex="0" title="Close"></i>`;
      card.appendChild(footer);

      card.style.display = 'block';
      card.style.fontSize = '0.8rem';
      card.style.visibility = visible ? 'visible' : 'hidden';

      const card_position_vs = R.cards[R.vs] ?? {};
      const position = card_position_vs[key] ?? {};

      let left = position.left;
      let top = position.top;

      if (left === undefined || top === undefined) {
          const panel = R.els.chartPanel.getBoundingClientRect();
          const { x1, y1 } = ClientXY(row.key);
          const START_BOUND = panel.left + 65;
          const STOP_BOUND = panel.right - 10;
          const HALF_WIDTH = width / 2;

          left ??= x1 - HALF_WIDTH;
          if (left < START_BOUND) {
              left = START_BOUND;
          } else if (left + width > STOP_BOUND) {
              left = STOP_BOUND - width;
          }
          top ??= panel.top + 30;
      }

      R.cards[R.vs] ??= {}, R.cards[R.vs][key] = {left: left, top: top};

      card.style.left = String(Math.round(left)) + 'px';
      card.style.top  = String(Math.round(top)) + 'px';
      document.body.appendChild(card);

      card = attachConnector(card)
      Draggable(card, header);

      return card
    };
    const assembleKV = (k, v, tabulate=false) => { return tabulate ? `<tr><td>${k}</td><td>${v}</td></tr>` : `${k}: ${v}`};
    const assembleCountScore = (row, tabulate=false) => {
      let text = [];
      text.push(assembleKV(`<b>${R.x.replace('zscore_', '')}`, `${row[R.x.replace('zscore_', 'count_')]} (${row[R.x].toFixed(2)})</b>`, tabulate));
      text.push(assembleKV(`<b>${R.y.replace('zscore_', '')}`, `${row[R.y.replace('zscore_', 'count_')]} (${row[R.y].toFixed(2)})</b>`, tabulate));
      const scores = R.scoreColumns.filter(c => (c !== R.x && c !== R.y));
      for (const c of scores) text.push(assembleKV(`${c.replace('zscore_', '')}`, `${row[c.replace('zscore_', 'count_')]} (${row[c].toFixed(2)})`, tabulate));
      return text
    };
    const assembleHoverText = (row) => {
      let text = [`<b>${row.compound}</b>`];
      text.push(...assembleCountScore(row));
      text.push(`nHH: ${row.nHH}`)
      text.push(`Encodings: ${row.copies}`)
      return text.join('<br>')
    };
    const assemblePlainText = (row) => {
      const bb = { c1_smiles: 'BB1', c2_smiles: 'BB2', c3_smiles: 'BB3', SMILES: 'SMILES'}
      let smiles = Object.entries(R.config.render).map(([k, v]) => `${bb[k]}: ${row?.[k] ?? ''}`);
      let ss = [`Compound: ${row.compound}`, ...smiles]
      const parts = assembleHoverText(row).split('<br>');
      const tags = ['<b>', '</b>'];
      for (let i=1; i < parts.length; i++) {
        let s = parts[i];
        for (const tag of tags) {
          s = s.replace(tag, '');
        }
        ss.push(s)
      }
      return ss.join('\n')
    };
    const alignModebarWithLegend = () => {
      const mb  = R.els.chartPanel.querySelector('.modebar');
      const leg = R.els.chartPanel.querySelector('.legend');
      if (!mb || !leg) return;

      const gbox = R.els.chartPanel.getBoundingClientRect();
      const lbox = leg.getBoundingClientRect();
      const top  = Math.max(0, Math.round(lbox.top - gbox.top)) + 3;
      mb.style.top = top + 'px';
    };
    const buildColumns = () => {
      const smiles = {title: 'Structure', columns: [], hozAlign: 'center'}
      for (const c of R.smilesColumns) {
        smiles.columns.push({
          title: c.includes('_') ? c.split('_')[0].replace('c', 'BB') : c,
          field: c,
          width: R.config.structure.width,
          formatter: (cell) => {
            return `<div data-smiles="${cell.getValue()}" data-height="125"></div>`;
          }
        });
      }

      const metrics = { title: 'Count & z-score', columns: [] };
      for (const countCol of R.countColumns) {
          const sample = countCol.replace('count_', '');
          const scoreCol = `zscore_${sample}`;
          metrics.columns.push({
              title: sample,
              field: countCol,
              variableHeight: true,
              formatter: (cell) => {
                  const row = cell.getRow().getData();
                  const count = row[countCol];
                  const zscore = row[scoreCol];
                  return `
                      <div style="line-height:1.2;padding:2px 0;">
                          <div><strong>${count ?? ''}</strong></div>
                          <div class="text-secondary">(${zscore ?? ''})</div>
                      </div>
                  `;
              }
          });
      }

      const EnHH = {
        title: 'Encodings',
        columns: [
          {
            title: 'nHH',
            field: 'copies',
            formatter: (cell) => {
                  const row = cell.getRow().getData();
                  const encodings = row['copies'] ? row['copies'] : 1;
                  const hh = (row.history_hits ?? "").toString().split(',').filter(Boolean).length;
                  return `
                      <div style="line-height:1.2;padding:2px 0;">
                          <div><strong>${encodings}</strong></div>
                          <div class="text-secondary">(${hh})</div>
                      </div>
                  `;
              }
          }
        ]
      }

      return [{title: 'Library', field: 'library'}, smiles, metrics, EnHH]
    };
    const tabularize = (el, data, columns, layout='fitColumns', rowHeight=125) => {
      el.replaceChildren();
      return  new Tabulator(el, {
        data: data,
        columns: columns,
        index: 'key',
        layout: layout,
        height: '100%',
        rowHeight: rowHeight,
        nestedFieldSeparator: "->",
        columnDefaults: { hozAlign: "center",  vertAlign: "middle", headerHozAlign: "center", headerSort: false },
      });
    };
    const updateHitsCount = (n, btn, span) => {
      if (n > 0) {
        btn.classList.remove('disabled');
        span.innerText = `${n}`;
      } else {
        btn.classList.add('disabled');
        span.innerText = '';
      }
    };
    const getCompactTimestamp = () => {
      const now = new Date();
      const pad = (num) => String(num).padStart(2, '0');
      const mm = pad(now.getMonth() + 1);
      const dd = pad(now.getDate());
      const yyyy = now.getFullYear();
      const hh = pad(now.getHours());
      const min = pad(now.getMinutes());
      const ss = pad(now.getSeconds());
      return `${mm}${dd}${yyyy}${hh}${min}${ss}`;
    };
    const hideUploadModal = () => {
      const modalElement = R.els.uploadModal;
      const instance = bootstrap.Modal.getInstance(modalElement);
      if (instance) {
          instance.hide();
      } else {
          modalElement.close();
      }
    };
    const showUploadPanel = () => {
      R.els.switchers.classList.add('d-none');
      R.els.selectors.classList.add('d-none');
      R.els.chartPanel.classList.add('d-none');
      hideUploadModal();
      R.els.uploadPanel.classList.remove('d-none');
      Object.assign(R, DATA);
    }
    const bbCards = (row) => {
      let smiles = getSMILES(row);
      if (Array.isArray(smiles) && smiles.length > 1) {
          smiles = [...smiles.slice(1), smiles[0]];
      }

      let cards = [];
      const n = smiles.length;
      for (let i = 0; i < n; i++) {
        const s = smiles[i];
        let label = '';
        if (i < n) {
          label = `BB${i+1}`;
        }
        if (i === n - 1) {
          label = 'SMILES';
        }

        const h = `
          <div class="col">
              <div class="card p-2 h-100">
                <div class="text-muted text-center small">${label}</div>
                <div class="flex-grow-1 d-flex align-items-center justify-content-center" 
                     data-smiles="${s}" data-height="200"></div>
              </div>
          </div>
        `;
        cards.push(h)
      }
      return cards;
    }

    return { findColumns, assembleCompoundCard, getSMILES, assemblePlainText,
             assembleKV, assembleCountScore, assembleHoverText, alignModebarWithLegend,
             buildColumns, keyForRow, tabularize, updateHitsCount, showUploadPanel,
             removeCompounds, removeCompound, showCompound, showCompounds, getCompactTimestamp, bbCards
           };
  })();

  const buildSelector = () => {
    const buildOptions = (selector, options, tag, selected, btn) => {
      if (!selector) return;
      selector.innerHTML = '';
      options.forEach(option => {
        const li = document.createElement('li');
        const link = document.createElement('a');
        link.href = '#';
        link.textContent = `${tag}: ${option}`;
        link.setAttribute('data-value', option);
        link.setAttribute('class', option === selected ? 'dropdown-item active' : 'dropdown-item');
        li.appendChild(link);
        selector.appendChild(li);
      })
      btn.textContent = `${tag}: ${selected}`;
      R[tag.toLowerCase()] = selected;
    };

    let x;
    let y;
    const scoreColumns = R.scoreColumns;

    if (R.x && R.y) {
      x = R.x;
      y = R.y;
    } else if (scoreColumns.includes('zscore_NTC')) {
      x = scoreColumns.find(element => element !== 'zscore_NTC');
      y = 'zscore_NTC';
    } else {
      x = scoreColumns[0];
      y = scoreColumns[1];
    }
    R.x = x;
    R.y = y;

    buildOptions(R.els.librarySel, ['All', ...R.libraries], 'Library', R.library || 'All', R.els.btnLibrary)
    R.els.selectors.classList.remove('d-none')
    R.vs = `${x.replace('zscore_', '')}.vs.${y.replace('zscore_', '')}`
  };

  /**
   * Dynamically builds and displays an options menu floating next to the clicked axis title
   * @param {HTMLElement} anchorElement - The axis title DOM element that was clicked
   * @param {'X'|'Y'} axisTag - Which axis we are updating
   */
  function showAxisContextMenu(anchorElement, axisTag) {
    // Remove any existing dynamic menu to prevent duplicates
    document.getElementById('plotly-axis-context-menu')?.remove();

    // Create a floating menu container (using Bootstrap dropdown styles)
    const menu = document.createElement('ul');
    menu.id = 'plotly-axis-context-menu';
    menu.className = 'dropdown-menu show p-2 shadow-sm';
    menu.style.position = 'absolute';
    menu.style.zIndex = '1050';
    menu.style.margin = '0';

    const currentSelected = axisTag === 'X' ? R.x : R.y;

    // Build options list inside the floating menu
    R.scoreColumns.forEach(option => {
      const li = document.createElement('li');
      const link = document.createElement('a');
      link.href = '#';
      link.textContent = option;
      link.className = option === currentSelected ? 'dropdown-item active' : 'dropdown-item';

      link.onclick = (e) => {
        e.preventDefault();

        // Update global state
        if (axisTag === 'X') R.x = option;
        if (axisTag === 'Y') R.y = option;

        // Re-trigger selector math & chart redraw
        buildSelector();
        analyzeData();
        buildHitsTable();
        buildTopHitsTable();
        renderChart();

        menu.remove(); // Close menu after choice
      };

      li.appendChild(link);
      menu.appendChild(li);
    });

    document.body.appendChild(menu);

    // Get dimensions of both the clicked title and our new menu
    const rect = anchorElement.getBoundingClientRect();
    const menuRect = menu.getBoundingClientRect();

    let leftPosition = 0;
    let topPosition = 0;

    if (axisTag === 'X') {
      // 💡 ABOVE THE TITLE:
      // Align horizontally with the center of the title, push vertically above it by its own height
      leftPosition = rect.left + (rect.width / 2) - (menuRect.width / 2);
      topPosition = rect.top - menuRect.height - 8; // 8px padding buffer
    } else if (axisTag === 'Y') {
      // 💡 RIGHT NEXT TO THE TITLE:
      // Place to the right of the title element, center it vertically against the title text
      leftPosition = rect.right + 8; // 8px padding buffer
      topPosition = rect.top + (rect.height / 2) - (menuRect.height / 2);
    }

    // Fallback check: Ensure the menu doesn't fly off the left or top edge of the screen
    leftPosition = Math.max(8, leftPosition);
    topPosition = Math.max(8, topPosition);

    // Apply calculated pixel coordinates (factoring in page scroll)
    menu.style.left = `${leftPosition + window.scrollX}px`;
    menu.style.top = `${topPosition + window.scrollY}px`;

    // Close menu if user clicks anywhere else on the screen
    setTimeout(() => {
      const closeMenu = (e) => {
        if (!menu.contains(e.target)) {
          menu.remove();
          document.removeEventListener('click', closeMenu);
        }
      };
      document.addEventListener('click', closeMenu);
    }, 0);
  }

  const handleChartEvent = (gd) => {
    if (R.library !== 'All') {
      R.utilities.alignModebarWithLegend();

      function bindAxisTitleClicks() {
        const xTitles = gd.querySelectorAll('.g-xtitle, .g-x2title');
        const yTitles = gd.querySelectorAll('.g-ytitle, .g-y2title');

        xTitles.forEach(el => {
          el.onclick = function(e) {
            e.stopPropagation();
            showAxisContextMenu(this, 'X');
          };
        });

        yTitles.forEach(el => {
          el.onclick = function(e) {
            e.stopPropagation();
            showAxisContextMenu(this, 'Y');
          };
        });
      }

      gd.removeAllListeners?.('plotly_click');
      gd.removeAllListeners?.('plotly_afterplot');

      if (typeof gd.on === 'function') {
          bindAxisTitleClicks();

          gd.on('plotly_click', ev => {
            const id = ev.points[0].id;
            const card = document.querySelector(`#${CSS.escape(id)}.card`)
            card ? R.utilities.removeCompound(id) : R.utilities.showCompound(id);
          });

      // Ensure that if the plot re-renders (zoom, pan, data changes), we re-bind
      gd.on('plotly_afterplot', () => {
        bindAxisTitleClicks();
      });

      }
    } else {
      gd.on('plotly_clickannotation', function(event) {
        const library = event.annotation.text;
        R.els.btnLibrary.textContent = `Library: ${library}`;
        R.els.btnLibrary.dataset.value = library;
        R['library'] = library;
        R.vs = `${R.x}.vs.${R.y}`

        analyzeData();
        buildHitsTable();
        buildTopHitsTable();
        renderChart();
      });
    }
  };

  const renderChart = () => {
    const x = R.x;
    const y = R.y;
    const library = R.library || R.els.btnLibrary?.textContent.split(': ')[1];
    const cfg = R.config;

    const colorForAxis = (ax) => {
      if ([0, 1, 2].includes(ax)) return cfg.colors.mono;
      if ([3, 4, 5].includes(ax)) return cfg.colors.di;
      return cfg.colors.tri;
    };
    const makeTrace = (name, rows, size=30) => {
      if (!Array.isArray(rows) || rows.length === 0) return null;
      return {
        ids: name.includes('sython') ? rows.map(r => r.key) : [],
        name: name,
        type: 'scattergl',
        mode: 'markers',
        marker: { color: rows.map(r => colorForAxis(r?.axis ?? 0)), size: size},
        x: rows.map(r => r?.[x]),
        y: rows.map(r => r?.[y]),
        text: name.includes('sython') ? rows.map(r => R.utilities.assembleHoverText(r)) : '',
        hoverinfo: 'text',
        hovertemplate: `%{text}<extra></extra>`,
        showlegend: name.includes('sython')
      };
    };

    let layout = {
      margin: {l: 60, r: 5, t: 5, b: 70},
      hovermode: 'closest',
      legend: {orientation: 'h', x: 0, xanchor: 'left', y: 0.99, yanchor: 'top'},
      font: { family: cfg.font, size: cfg.fontSize }
    };
    let config = {responsive: true, displayModeBar: true, displaylogo: false}
    let traces = [];

    if (library === 'All') {
      const libraries = R.libraries.filter(v => v && v !== 'All');
      const n = libraries.length;
      const columns = Math.ceil(Math.sqrt(n));
      const rows = Math.ceil(n / columns);
      const lastRowCount = n - (rows - 1) * columns;
      const xsAll = R.uniques.map(r => r[x]);
      const ysAll = R.uniques.map(r => r[y]);

      const xMin = xsAll.length ? Math.min(0, ...xsAll) : 0;
      const xMax = xsAll.length ? Math.max(0, ...xsAll) : 1;
      const xOffset = (xMax - xMin) * 0.05
      const yMin = ysAll.length ? Math.min(0, ...ysAll) : 0;
      const yMax = ysAll.length ? Math.max(0, ...ysAll) : 1;
      const yOffset = (yMax - yMin) * 0.05

      const axisDefault = {
        mirror: true, nticks: 3, zeroline: false, showgrid: false, linecolor: '#939393',
        linewidth: 1, ticks: 'outside', ticklen: 1, tickwidth: 1, automargin: true,
      };

      config.staticPlot = true;
      layout.margin.l = 90;
      layout.grid = {rows, columns: columns, pattern: 'independent', xgap: 0.03, ygap: 0.05, roworder: 'top to bottom'};
      layout.annotations = [];

      for (let i = 0; i < n; i++) {
        const library = libraries[i];
        const subset = R.uniques.filter(r => r.library === library);
        const maximum = Math.max(...subset.map(r => Number(r[x])));
        const axisKey = k => k.replace(/^([xy])(.*)$/, '$1axis$2');
        const trace = makeTrace(library, subset, 15);

        const ax = i + 1;
        trace.xaxis = (ax === 1) ? 'x' : `x${ax}`;
        trace.yaxis = (ax === 1) ? 'y' : `y${ax}`;
        traces.push(trace);

        const rIdx = Math.floor(i / columns);
        const cIdx = i % columns;

        const isLastRow = (rIdx === rows - 1);
        const isSecondLastRow = (rIdx === rows - 2);
        const lastRowIncomplete = (lastRowCount > 0 && lastRowCount < columns);
        const showXticks = isLastRow || (lastRowIncomplete && isSecondLastRow && (cIdx >= lastRowCount));
        const showYticks = (cIdx === 0);

        layout[axisKey(trace.xaxis)] = Object.assign(layout[axisKey(trace.xaxis)] || {}, {
          showticklabels: !!showXticks, range: [xMin - xOffset, xMax + xOffset], ...axisDefault
        });
        layout[axisKey(trace.yaxis)] = Object.assign(layout[axisKey(trace.yaxis)] || {}, {
          showticklabels: !!showYticks, range: [yMin - yOffset, yMax + yOffset], ...axisDefault
        });

        // layout.title = {text: library, x: 0.5, xanchor: 'center', y: 0.98, yanchor: 'top',
        // font: {size: R.config.fontSize, family: R.config.fontFamily, color: maximum >= 1 ? '#ff0000' : '#000000'}}

        layout.annotations.push({
            text: String(library),
            xref: `${trace.xaxis} domain`,
            yref: `${trace.yaxis} domain`,
            x: 0.5,
            y: 0.98,
            xanchor: 'center',
            yanchor: 'top',
            showarrow: false,
            captureevents: true,
            font: {
                size: R.config.fontSize,
                family: R.config.fontFamily,
                color: maximum >= 1 ? '#ff0000' : '#000000'
            }
        });
      }

      layout.annotations.push(
        {
          text: `${x.replace('zscore_', '')} (z-score)`,
          xref: 'paper', yref: 'paper',
          x: 0.5, y: -0.04, showarrow: false, font: {size: R.config.fontSize}
        },
        {
          text: `${y.replace('zscore_', '')} (z-score)`,
          xref: 'paper', yref: 'paper',
          x: -0.05, y: 0.5, textangle: -90, showarrow: false, font: {size: R.config.fontSize}
        }
      );
    } else {
      const rows = R.uniques.filter(r => (r.library ?? 'Unknown') === library);
      traces.push(makeTrace('Mono-sython', rows.filter(r => ([0, 1, 2].includes(r.axis)))))
      traces.push(makeTrace('Di-sython', rows.filter(r => ([3, 4, 5].includes(r.axis)))))
      traces.push(makeTrace('Tri-sython', rows.filter(r => (r.axis === 6))))
      const axisDefault = {zeroline: false, showgrid: false, mirror: true, linecolor: '#939393', linewidth: 1, nticks: 5}

      const xs = rows.map(r => r?.[x]).filter(Number.isFinite);
      const ys = rows.map(r => r?.[y]).filter(Number.isFinite);

      if (xs.length > 0 && ys.length > 0) {
        let lo = Math.min(Math.min(...xs), Math.min(...ys));
        let hi = Math.max(Math.max(...xs), Math.max(...ys));
        const span = Math.max(hi - lo, 1);
        const pad = span * 0.01;
        lo -= pad;
        hi += pad;

        const range = [lo, hi];
        layout.xaxis = {title: {text: `${x.replace('zscore_', '')} (z-score)`}, range: range, ...axisDefault};
        layout.yaxis = {title: {text: `${y.replace('zscore_', '')} (z-score)`}, range: range, ...axisDefault};

        layout.shapes = [{
            type: 'line', xref: 'x', yref: 'y', x0: lo, y0: lo, x1: hi, y1: hi,
            line: { dash: 'dot', width: 1, color: '#d7d7d7' }, layer: 'below',
            _tag: 'diag_y_eq_x'
        }];
      } else {
        layout.xaxis = {title: {text: `${x.replace('zscore_', '')} (z-score)`}, ...axisDefault};
        layout.yaxis = {title: {text: `${y.replace('zscore_', '')} (z-score)`}, ...axisDefault};
      }
    }

    Plotly.react(R.els.chartPanel, traces.filter(r => r), layout, config)
      .then(handleChartEvent)
      .catch(err => console.error('Failed to make chart: ', err));
  };

  const buildHitsTable = () => {
    const data = Array.from(R.hits.get(R.vs) || []);
    let table = R.hitsTable;
    if (table) {
      table.replaceData(data);
    } else {
      const deleteCol = {
      title: "", width: 46, hozAlign: "center", headerSort: false,
      titleFormatter: () => `<i class="bi bi-trash text-danger" aria-label="Delete row"></i>`,
      formatter: () => `<button type="button" class="btn btn-sm btn-outline-danger" 
                          title="Delete row" data-action="del"><i class="bi bi-trash"></i></button>`,
      cellClick: (e, cell) => {
        const btn = e.target.closest('button[data-action="del"]');
        const table = R.hitsTable;
        if (!btn) return;
        e.preventDefault();
        e.stopPropagation();
        const row = cell.getRow();
        const data = row.getData();
        const hitsSet = R.hits.get(R.vs);

        if (hitsSet) {
            for (const hitRow of hitsSet) {
                if (hitRow.key === data.key) {
                    hitsSet.delete(hitRow);
                    break;
                }
            }
        }

        row.delete();
        R.utilities.updateHitsCount(table, R.els.btnHitsModal, R.els.numHits);

        data.hits = false;
        R.topHitsTable.updateData([data]);
        const key = row.getData().key;
        const icons = document.querySelectorAll(`[data-action="bag"][data-key="${key}"]`);
        Array.from(icons).forEach(icon => {
          icon.classList.replace('bi-bag-fill', 'bi-bag');
          icon.classList.remove('text-danger');
        })
      },
    };
      const columns = [deleteCol, ...R.utilities.buildColumns()];
      table = R.utilities.tabularize(R.els.hitsTable, data, columns);
    }
    R.hitsTable = table;
    R.utilities.updateHitsCount(data.length, R.els.btnHitsModal, R.els.numHits)
  };

  const buildTopHitsTable = () => {
    const hits = R.hitsTable?.getData?.() ?? [];
    const keys = hits.map(hit => hit.key);
    const tops = Array.from(R.tops.get(R.vs) || []).map(top => ({...top, hits: keys.includes(top.key)}));
    let table = R.topHitsTable;
    if (table) {
      table.replaceData(tops)
    } else {
      const deleteColumn = {
        title: "", width: 46, hozAlign: "center", headerSort: false,
        titleFormatter: () => `<i class="bi bi-trash text-danger" aria-label="Delete row"></i>`,
        formatter: () => `<button type="button" class="btn btn-sm btn-outline-danger" 
                            title="Delete row" data-action="del">
                            <i class="bi bi-trash"></i>
                          </button>`,
        cellClick: (e, cell) => {
          const btn = e.target.closest('button[data-action="del"]');
          if (!btn) return;
          e.preventDefault();
          e.stopPropagation();
          const row = cell.getRow();
          row.delete();
          R.utilities.updateHitsCount(table.getDataCount(), R.els.btnTopHitsModal, R.els.numTopHits);
          const key = cell.getRow().getData().key
          R.utilities.removeCompound(key);
        },
      };
      const hitsColumn = { title: "Hits", field: "hits", formatter:"tickCross",
        accessorClipboard: v => (v ? 1 : 0), accessorDownload:  v => (v ? 1 : 0),
        cellClick: (e, cell) => {
          const next = !Boolean(cell.getValue());
          cell.setValue(next, true);
          const data = cell.getRow().getData();
          const hitsSet = R.hits.get(R.vs) ?? R.hits.set(R.vs, new Set()).get(R.vs);

          const table = R.hitsTable;
          if (next) {
            table.updateOrAddData([data]).then( () => {
              hitsSet.add(data);
              R.utilities.updateHitsCount(R.hitsTable.getDataCount(), R.els.btnHitsModal, R.els.numHits);
            })
          } else {
            table.deleteRow(data.key);
            for (const row of hitsSet) {
                if (row.key === data.key) {
                    hitsSet.delete(row);
                    break;
                }
            }
            R.utilities.updateHitsCount(R.hitsTable.getDataCount(), R.els.btnHitsModal, R.els.numHits);
          }
        }
      }
      const columns = [deleteColumn, hitsColumn, ...R.utilities.buildColumns()]
      table = R.utilities.tabularize(R.els.topHitsTable, tops, columns);
    }
    R.topHitsTable = table;
    R.utilities.updateHitsCount(tops.length, R.els.btnTopHitsModal, R.els.numTopHits);
  };

  const analyzeData = () => {
    const norm = v => {
      const n = Number(v);
      return Number.isFinite(n) ? n : -Infinity;
    };

    if (R.uniques.length === 0) {
      const groups = new Map();               // key -> indices[]
      const groupSizeByKey = new Map();       // key -> count
      const bestByKey = new Map();            // key -> { idx, score }

      R.rows.forEach((row) => {
        const k = String(row?.key ?? "");
        if (!groups.has(k)) groups.set(k, []);
        groups.get(k).push(row.index);

        groupSizeByKey.set(k, (groupSizeByKey.get(k) ?? 0) + 1);

        const s = norm(row?.[R.x]);
        const prev = bestByKey.get(k);
        if (!prev || s > prev.score /* tie: keep first seen */) {
          bestByKey.set(k, { idx: row.index, score: s });
        }
      });

      R.uniques = [];
      for (const [k, info] of bestByKey.entries()) {
        const copies = groupSizeByKey.get(k) ?? 0;
        R.uniques.push({ ...R.rows[info.idx], copies });
      }

      R.duplicates = Object.fromEntries([...groups.entries()]
        .filter(([, idxs]) => idxs.length >= 2).map(([k, idxs]) => [k, idxs.slice()]));
    }

    const topsSet = new Set();
    const topsByLibrary = {};

    for (const row of R.uniques) (topsByLibrary[row.library] ??= []).push(row);

    for (const v of Object.values(topsByLibrary)) {
      const top = v.filter(a => (a.axis === 6) && (norm(a?.[R.x]) >= R.config.zscoreCutoff))
        .sort((a, b) => norm(b?.[R.x]) - norm(a?.[R.x])).slice(0, R.config.nTopHits);
      for (const t of top) {
        topsSet.add(t);
      }
    }

    if (!R.tops.has(R.vs)) {
      R.tops.set(R.vs, topsSet);
    } else {
      R.tops.set(R.vs, new Set([...R.tops.get(R.vs), ...topsSet]));
    }

    R.hits.set(R.vs, R.hits.get(R.vs) ?? new Set());
  };

  const bindEvents = () => {
    R.els.dz = Dropzone.forElement("#dropzone");
    R.els.dz.on('addedfile', async (file) => {
      if (!file) return;

      const dz = R.els.dz;
      try {
        const data = await R.io.load(file)
        dz.removeAllFiles(true);

        if (file.name.endsWith('.json') || file.name.endsWith('.json.gz')) {
          Object.assign(R, data);
          R.uniques = R.uniques.filter(Boolean);
          R.hits = new Map(Object.entries(R.hits).map(([vsKey, rowArray]) => [vsKey, new Set(rowArray)]));
          R.tops = new Map(Object.entries(R.tops).map(([vsKey, rowArray]) => [vsKey, new Set(rowArray)]));
          initializePage(R.rows);
        } else {
            initializePage(data);
        }
      } catch (error) {
        R.els.dz.removeAllFiles(true);
        R.onError?.(error) || (() => {
          R.els.chartPanel.innerHTML = `<div class="text-danger text-center fs-2 fw-bold pt-5">
            Data Load Error: ${error?.message || String(error)}</div>`;
        })();
      }
    });

    const bindDropdown = (menu, btn, label) => {
      menu.addEventListener('click', (e) => {
        const link = e.target.closest('a.dropdown-item');
        if (!link) return;
        e.preventDefault();

        const value = link.dataset.value || link.textContent.trim().split(': ')[1];
        btn.textContent = `${label}: ${value}`;
        btn.dataset.value = value;

        menu.querySelectorAll('.dropdown-item.active').forEach(a => a.classList.remove('active'));
        link.classList.add('active');

        R[label.toLowerCase()] = value;
        R.vs = `${R.x}.vs.${R.y}`

        try { bootstrap.Dropdown.getOrCreateInstance(btn).hide(); } catch (_) {}

        analyzeData();
        buildHitsTable();
        buildTopHitsTable();
        R.utilities.removeCompounds();
        renderChart();
        if (R.library !== 'All') R.utilities.showCompounds();
      });
    }

    bindDropdown(R.els.librarySel, R.els.btnLibrary, 'Library')

    // R.els.btnSaveConfig?.addEventListener('click', () => {
    //   const cfg = readConfigForm();
    //   bootstrap.Modal.getInstance(R.els.configModal)?.hide();
    //   GlobalConfig(cfg);
    //   renderChart();
    // });
    //
    // R.els.configModal.addEventListener('show.bs.modal', () => {
    //   if (R.config) populateConfigForm(R.config);
    //   bootstrap.Modal.getInstance(R.els.configModal)?.show();
    // });

    handleChartEvent(R.els.chartPanel);

    const clamp = (v, min, max) => Math.max(min, Math.min(max, v));
    const showToastAbove = (target, { message = 'Done!', delay = 3000, offsetY = 8} = {}) => {
      const toastEl = document.createElement('div');
      toastEl.className = 'toast shadow';
      toastEl.setAttribute('role', 'status');
      toastEl.setAttribute('aria-live', 'polite');
      toastEl.setAttribute('aria-atomic', 'true');
      toastEl.style.position = 'absolute';
      toastEl.style.zIndex = 2000;
      toastEl.style.visibility = 'hidden';
      toastEl.innerHTML = `<div class="toast-body py-2 px-3">${message}</div>`;
      document.body.appendChild(toastEl);

      const rect = target.getBoundingClientRect();
      const { scrollX, scrollY, innerWidth } = window;
      const tw = toastEl.offsetWidth;
      const th = toastEl.offsetHeight;
      let left = rect.left + rect.width / 2 - tw / 2 + scrollX;
      left = clamp(left, 8 + scrollX, innerWidth - tw - 8 + scrollX);
      const top = rect.top + scrollY - th - offsetY;

      toastEl.style.left = `${left}px`;
      toastEl.style.top  = `${Math.max(8 + scrollY, top)}px`;
      toastEl.style.visibility = 'visible';

      const toast = bootstrap.Toast.getOrCreateInstance(toastEl, {autohide: true, delay});
      toastEl.addEventListener('hidden.bs.toast', () => toastEl.remove());
      toast.show();
    };

    document.addEventListener('click', (e) => {
      const btn = e.target.closest('[data-action]');
      if (!btn) return;

      const action = btn.dataset.action;
      const key = btn.dataset.key;
      const actionsNeedKey = ['encoding', 'hh', 'bag', 'copy', 'email', 'close'];
      if (actionsNeedKey.includes(action) && !key) {
        console.error(`Action "${action}" requires a key, but none was found.`);
        return;
      }

      const row = R.uniques.find(r => r.key === key);
      const hitsSet = R.hits.get(R.vs) ?? R.hits.set(R.vs, new Set()).get(R.vs);

      switch (action) {
        case 'encoding': {
          const cards = R.utilities.bbCards(row);
          R.els.encodingSMILES.innerHTML = cards.join('\n');

          const ids = R.duplicates[key];
          const rows = R.rows.filter(r => ids.includes(r.index));
          const excludes = ['Library', 'Structure', 'Encodings']
          const columns = R.utilities.buildColumns().filter(x => !excludes.includes(x.title));
          R.utilities.tabularize(R.els.encodingTable, rows, columns, 'fitColumns', 40);

          const title = document.getElementById('encodingModalTitle');
          title.textContent = `${row.compound} Encoding Details (n=${rows.length})`;

          const modal = bootstrap.Modal.getOrCreateInstance(R.els.encodingModal, {
            backdrop: true,
            keyboard: true,
            focus: true
          });
          modal.show();
          break;
        }
        case 'hh': {
          const title = document.getElementById('hhModalTitle');
          title.textContent = `${row.compound} History Hits (n=${row.nHH})`;

          const cards = R.utilities.bbCards(row);

          R.els.hhSMILES.innerHTML = cards.join('\n');

          const rows = row.history_hits.split(';').map((group, index) => {
            const [ids, value] = group.split(':');
            const [selection, target] = ids.split('_');
            return { key: index, selection: selection, target: target, count: parseInt(value, 10)};
          });
          const columns = [{title: 'Selection', field: 'selection'}, {title: 'Target', field: 'target'}, {title: 'Count', field: 'count'}];
          R.utilities.tabularize(R.els.hhTable, rows, columns, 'fitColumns', 20);
          const modal = bootstrap.Modal.getOrCreateInstance(R.els.hhModal, {
            backdrop: true,
            keyboard: true,
            focus: true
          });
          modal.show();
          break;
        }
        case 'bag': {
          if (btn.classList.contains('bi-bag-fill')) {
            R.hitsTable.deleteRow(key);
            for (const hitRow of hitsSet) {
                if (hitRow.key === key) {
                    hitsSet.delete(hitRow);
                    break;
                }
            }
            btn.classList.replace('bi-bag-fill', 'bi-bag');
            btn.classList.remove('text-danger');
          } else {
            R.hitsTable.updateOrAddData([row]);
            hitsSet.add(row);
            btn.classList.remove('bi-bag');
            btn.classList.add('bi-bag-fill');
            btn.classList.add('text-danger');
          }
          R.utilities.updateHitsCount(R.hitsTable.getDataCount(), R.els.btnHitsModal, R.els.numHits);
          break;
        }
        case 'copy': {
          const text = R.utilities.assemblePlainText(row);
          navigator.clipboard?.writeText(text).catch(console.warn);
          showToastAbove(btn, { message: 'Compound info copied!'});
          break;
        }
        case 'email':
          const text = R.utilities.assemblePlainText(row);
          const subject = encodeURIComponent('Compound info');
          const body = encodeURIComponent(text);
          window.location.href = `mailto:?subject=${subject}&body=${body}`;
          break;
        case 'close':
          R.utilities.removeCompound(key)
          break;
        case 'downloadTopHits':
          R.topHitsTable.download('csv', `${R.vs}.top.hits.csv`)
          break;
        case 'downloadHits':
          R.hitsTable.download('csv', `${R.vs}.candidate.hits.csv`)
          break;
        case 'saveSession':
          R.io.saveSession();
          break;
        case 'uploadAnyway':
          R.utilities.removeCompounds();
          R.utilities.showUploadPanel();
          break;
        case 'saveAndUpload':
          R.utilities.removeCompounds();
          R.io.saveSession();
          R.utilities.showUploadPanel();
          break;
        case 'applyConfig':
          bootstrap.Modal.getInstance(R.els.configModal)?.hide();
          const cfg = readConfigForm();
          R.config = cfg
          analyzeData();
          buildHitsTable();
          buildTopHitsTable();
          R.utilities.removeCompounds();
          renderChart();
          if (R.library !== 'All') R.utilities.showCompounds();
          break;
        case 'postHits':
          const url = `${post_api}`;
          const tableData = R.hitsTable.getData();
          fetch(url, {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
              "X-CSRFToken": `${csrftoken}`
            },
            body: JSON.stringify({ rows: tableData })
          })
          .then(response => response.json())
          .then(data => {
              console.log(`data: ${data}`)
              const msgSpan = document.getElementById('hits-post-message');
              if (msgSpan) {
                  msgSpan.textContent = `Successfully posted ${data.hits} hits.`;
              }
          })
          .catch(err => {
              console.error('Error posting hits:', err);
              const msgSpan = document.getElementById('hits-post-message');
              if (msgSpan) {
                  msgSpan.textContent = `Failed to post hits.`;
              }
          });

          const hm = bootstrap.Modal.getOrCreateInstance(R.els.hitsModal, {
            backdrop: true,
            keyboard: true,
            focus: true
          });
          setTimeout(() => {
            hm.hide();
          }, 5000);
          break;
      }
    });
  };

  const initializePage = (rows) => {
    R.els.uploadPanel.classList.add('d-none');
    R.els.switchers.classList.remove('d-none');
    R.els.selectors.classList.remove('d-none');
    R.els.chartPanel.classList.remove('d-none');

    R.columns = Object.keys(rows[0]);
    if (!R.utilities) {
        throw new Error("R.utilities not yet initialized. Cannot process columns.");
    }
    R.countColumns = R.utilities.findColumns(R.columns, 'count_')
    R.scoreColumns = R.utilities.findColumns(R.columns, 'zscore_')
    R.smilesColumns = R.utilities.findColumns(R.columns, '', 'smiles')
    R.libraries = [...new Set((rows ?? []).map(r => r.library))];

    if (!R.rows || R.rows.length === 0) {
      R.rows = rows.map(row => {
          row.key = R.utilities.keyForRow(row);
          row.compound = row.compound || `V${row.index}`;
          row.nHH = (row.history_hits ?? "").toString().split(';').filter(Boolean).length;
          return row;
      });
    }

    buildSelector();
    analyzeData();
    buildHitsTable();
    buildTopHitsTable();
    renderChart();
  };

  global.Visualizer = {
    async init(input=null) {
      R.config = readConfigForm();
      bindEvents();

      if (typeof post_api === 'undefined' || !post_api) {
        document.getElementById('hitsModalFooter').style.display = 'none';
      }

      if (input) {
        const rows = await R.io.load(input);
        initializePage(rows);
      } else {
        R.els.uploadPanel.classList.remove('d-none')
      }
    },
  };
})(window);
