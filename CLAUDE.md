# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Tahoe-100M Gene Explorer — a GTEx-style search and visualization tool for a ~100M single-cell dataset. Users search a gene and see its expression as a **violin plot** and a **cell line × perturbation heatmap**.

The core design: precompute per-group aggregates offline (the pipeline), serve them at millisecond latency (the API), and render them in the browser (the frontend).

## Commands

### API (Python / FastAPI)

```bash
pip install -r api/requirements.txt

# Build the sample store (~12k cells bundled in the repo)
python pipeline/build_aggregates.py --sample

# Serve the API (reads pipeline/out by default)
uvicorn api.main:app --reload --port 8000

# Point at a specific store directory
STORE_DIR=/path/to/pipeline/out uvicorn api.main:app --port 8000

# Run tests (requires pipeline/out to exist; skips cleanly if absent)
pytest api/tests
```

### Frontend (React / Vite)

```bash
cd web
npm install
npm run dev       # dev server on :5173, proxies /api → localhost:8000
npm run build     # tsc then vite build
npm run lint      # oxlint
```

### Pipeline: larger stores

```bash
# Cancer panel (~122 genes) across 256 real shards (~10 min, 8 cores)
python pipeline/build_aggregates.py --gene-set cancer --shards 256 --workers 8

# Inline gene list
python pipeline/build_aggregates.py --genes TP53,EGFR,KRAS --shards 256
```

## Architecture

Three independent layers joined by a JSON contract (defined in `project.md` §4):

```
Tahoe-100M shards → pipeline/ (offline) → api/ (FastAPI) → web/ (React)
```

**`pipeline/build_aggregates.py`** — two execution engines:
- `single`: one DuckDB query with `quantile_cont` (exact deciles, doesn't scale to full dataset due to OOM)
- `parallel` (default for `--full`): multiprocessing map-reduce; each worker emits a mergeable fixed-bin histogram of expression values; a DuckDB reduce phase sums partials; deciles reconstructed with NumPy. Single data pass, bounded memory, checkpointable.

Output layout (written to `pipeline/out/`, git-ignored):
```
gene_index.parquet                    # symbol → name, for autocomplete
aggregates/gene_symbol=<SYM>/*.parquet  # one Hive partition per gene
run_manifest.json                     # provenance
```

**`api/store.py`** — `Store` class loads the gene index fully into memory at startup (so search is a pure in-memory scan). Gene partitions are read on demand via an in-process DuckDB connection. Supports **multi-root stores**: `[pipeline/out_demo, pipeline/out]` — first root to provide a gene's partition wins. The demo store (`out_demo`) provides richer pre-built data for 4 curated genes (KRAS, TP53, EGFR, BRCA1).

**`api/main.py`** — FastAPI app. Two contract endpoints:
- `GET /api/search?q=<prefix>` — autocomplete, prefix matches rank first
- `GET /api/gene/{symbol}` — heatmap matrix + violin array per group

Gene payloads are `@lru_cache`'d (immutable for a given store build).

**`web/src/api.ts`** — typed fetch wrappers. `VITE_API_URL` env var overrides the default `http://localhost:8000`. Demo genes (`KRAS`, `TP53`, `EGFR`, `BRCA1`) are loaded from static JSON files in `web/public/` so the frontend works without the API running.

**`web/src/App.tsx`** — single-page app: `SearchBar` autocompletes against `/api/search`, selecting a gene fetches `/api/gene/{symbol}` and renders `Heatmap` and `ViolinPlot` (both using Plotly).

## JSON contract (section 4 of project.md)

The shape returned by `GET /api/gene/{symbol}`:

```json
{
  "gene": "TP53",
  "ensembl_id": "ENSG...",
  "n_groups": 42,
  "heatmap": {
    "cell_lines": ["A549", ...],
    "perturbations": ["DMSO", ...],
    "mean": [[0.1, null, ...], ...],
    "cell_line_organs": {"A549": "lung", ...}
  },
  "violin": [
    {
      "cell_line": "A549", "organ": "lung",
      "perturbation": "DMSO",
      "n": 1000, "deciles": [0.0, 0.1, ..., 2.5],
      "pct_expressing": 0.85, "mean": 0.6
    }
  ]
}
```

`deciles` is always 11 values: q0, q10, …, q100. `null` in the heatmap matrix means the group was not observed.
