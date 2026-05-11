
# Enrichment Profile Visualizer (EPV) Usage

    Author: Dr. Fei Yuan (fei.yuan@bcm.edu)
    Version: 1.0
    Description: A lightweight, browser‑only app for visualizing enrichment profiles, 
                 exploring hits and rendered chemical structures.

---

## Overview
The **Enrichment Profile Visualizer (EPV)** is a single‑page web app for interactive inspection 
of enrichment data (e.g., DEL or HTS outputs). It supports:
- Drag‑and‑drop loading of CSV/TSV/JSON (optionally gzipped) file.
- Interactive Plotly scatter plot.
- Tabulated hit summaries via Tabulator.
- SMILES and building‑block rendering via RDKit.js.

All functionality runs client‑side in the browser, no backend required.

---

## Typical Workflow
1. Open `index.html` (or serve locally).
2. Load dataset (via drag-and-drop or **Upload** button).
3. Choose X/Library/Y fields using the selectors.
4. Zoom, hover, and click points for detail.
5. Review **Candidate Hits** and **Top Hits**.
6. Tune font sizes, structure rendering, z-score cutoff, and maximum number of hits per library in **Config**.
7. Open **Encoding** to inspect details of building blocks encoded by multiple barcodes.

---

## Data Format
Your table must include a header row and at least two numeric columns for plotting. Common/expected fields:

| Column | Description | Example |
|---|---|---|
| `library` | Library/category id | `qDOS1` |
| `zscore_x`, `zscore_y`, `...`| Values for X/Y axes (must start with zscore_) | `2.31` |
| `count_x`, `count_y`, `...` | Values paired with z-scores (must start with count_) | `10` |
| `c1_smiles`, `c2_smiles`, `c3_smiles` | Building‑block SMILES | `CCN(CC)C(=O)C1=CC=CC=C1` |
| `SMILES` | Full compound SMILES | `CCN(CC)C(=O)...` |

CSV/TSV or GZIP‑compressed versions are supported.

---

## UI Guide
### Toolbar (top)
- **Help ( ? )** — Open help modal
- **Candidate Hits (bag)** — Open the table of candidate hits
- **Top Hits (table)** — Open the table of top hits
- **Upload (cloud)** — Open the unsaved session warning/upload modal to load new data
- **Save (download)** — Save the current visualizer state as a JSON file
- **Config (gear)** — Open configuration modal

### Selectors (after data loads)
- **X**: choose the x‑axis field (numeric)
- **Library**: choose the grouping field
- **Y**: choose the y‑axis field (numeric)

### Chart Panel
Interactive Plotly scatter with zoom/pan/reset. Click points to inspect structures/details.

### Modals
- **Help Modal** — Overview
- **Candidate Hits** — All selected/candidate hits (Tabulator)
- **Top Hits** — AI‑flagged or top‑N per library
- **Config** — Visualization and rendering settings
- **Encoding** — Composition view (SMILES, BB1–BB3)
- **Upload** — Prompts to save the current session before loading new data

---

## Configuration (⚙️)
- **Font size**: for chart labels, axes, and compound card text
- **Structure size**: width × height (px)
- **Structure rendering toggles**: SMILES, BB1, BB2, BB3
- **Top hits per library**: integer (1–20)
- **z score cutoff**: threshold for determining 'hits'
- **Maximum number** of top hits per library: integer (1–20)
- **Colors**: mono/di/tri‑sython pickers
Click **Apply** to save/apply changes.

---

## Dependencies (CDN)
Loaded automatically by `index.html`:
- Bootstrap 5.3.3 + Bootstrap Icons 1.13.1
- Dropzone.js (drag‑and‑drop uploads)
- pako (gzip)
- PapaParse (CSV/TSV parsing)
- Plotly.js 3.1.0
- Tabulator 6.3.1
- RDKit.js (minimal build)
- marked (usage generating)

---

## Styling Notes
`style.css` sets:
- Chart container height to `calc(100vh - 65px)`
- Wide modals (95vw) and tall content (75vh) for tables
- Flat/transparent Tabulator theme with clear hover/selection states

Modify to taste.

---

## Troubleshooting
| Symptom | Cause | Fix |
|---|---|---|
| Upload unresponsive | Browser security/sandbox | Use Chrome/Edge or serve locally (`python -m http.server`) |
| Structures not rendering | RDKit not initialized | Reload after network is stable |
| Blank axes | No numeric columns | Include appropriate numeric fields |
| Empty hits | Filters too strict | Adjust **Top hits per library** in Config |

---

## Extensibility
- Change the default dataset URL in `Visualizer.init({ url: ... })`
- Extend chart logic and hit detection in `visualizer.js`
- Adjust SMILES rendering in `smiles-render.js`
- Style via `style.css`

---

## License
This project is licensed under the PolyForm Noncommercial License 1.0.0.
You are free to use, copy, modify, and distribute this software for noncommercial purposes — 
including academic research, internal analysis, and educational use — provided that this license 
notice is included in all copies or substantial portions of the software.
