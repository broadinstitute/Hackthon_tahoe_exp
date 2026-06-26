# Tahoe-100M Gene Explorer

A GTEx-style search & visualization tool for the **Tahoe-100M** single-cell
dataset. Search for a gene and see its expression rendered as a **violin plot**
and a **cell line × perturbation heatmap** — the equivalent of a single
[GTEx gene page](https://gtexportal.org/home/gene/TAGLN).

The dataset is ~100M single cells stored row-oriented (h5ad / Parquet shards),
so naive per-gene scans over every cell are far too slow for an interactive
site. The core idea is to **precompute per-group aggregates offline** and serve
them at millisecond latency.

## Architecture

Three layers, each owned end-to-end by one team member and joined by a frozen
JSON contract (see [`project.md`](project.md) §4):

```
 Tahoe-100M shards                                          Browser
        │                                                      ▲
        ▼                                                      │ JSON
 ┌──────────────┐  Parquet/   ┌──────────────┐   HTTP   ┌──────┴─────┐
 │ 1. Precompute │── DuckDB ──▶│ 2. Query API │ ────────▶│ 3. Frontend │
 │   pipeline    │   store     │  (FastAPI)   │          │  (planned)  │
 └──────────────┘             └──────────────┘          └────────────┘
   offline, batch               ms-latency               search + plots
```

1. **Precompute pipeline** (`pipeline/`) — streams the raw expression shards
   once and computes, for every `gene × cell_line × perturbation` group:
   `n_cells`, `pct_expressing`, `mean`, and 11 deciles (the violin silhouette).
   Output is Parquet partitioned by gene plus a gene index for autocomplete.
2. **Query API** (`api/`) — FastAPI service that reads the store and serves
   `GET /api/search?q=` (autocomplete) and `GET /api/gene/{symbol}` (heatmap +
   violin payload) at ms latency, with an in-memory gene index and LRU caching.
3. **Frontend** (`web/`, planned) — gene search box + gene page rendering the
   violin and heatmap (Plotly).

See [`project.md`](project.md) for the full design brainstorm: the precompute
trade-offs, the JSON contract, scale estimates, and deferred nice-to-haves.

## Repository layout

| Path                         | What                                                              |
| ---------------------------- | ---------------------------------------------------------------- |
| `pipeline/`                  | Precompute pipeline (`build_aggregates.py`) + tiny sample data    |
| `pipeline/sample_data/`      | Downsampled 12k-cell slice for dev without the full ~328 GB set   |
| `api/`                       | FastAPI Query API serving the precompute store                   |
| `project.md`                 | Design brainstorm, architecture, and the frozen JSON contract    |
| `pipeline/out/`              | Built store (git-ignored, regenerated from raw data)             |

## Quick start

Build the bundled sample store and serve it locally:

```bash
pip install -r api/requirements.txt

# 1. Build the precompute store from the bundled sample (-> pipeline/out/)
python pipeline/build_aggregates.py --sample

# 2. Serve the Query API (reads pipeline/out by default)
uvicorn api.main:app --reload --port 8000
```

Then open the interactive API docs at <http://localhost:8000/docs>.

### A bigger demo store: the curated cancer panel

For a richer store than the tiny sample, build the curated ~122-gene cancer
panel ([`pipeline/gene_sets/cancer.txt`](pipeline/gene_sets/cancer.txt)) over a
slice of the real shards. This uses the cheap targeted ("no-unnest") path, so
even 256 shards (~7M cells) finishes in ~10 min on an 8-core box:

```bash
# All 122 cancer genes across 256 real shards (-> pipeline/out/)
python pipeline/build_aggregates.py --gene-set cancer --shards 256 --workers 8
```

`--gene-set <name>` loads `pipeline/gene_sets/<name>.txt`; `--genes-file PATH`
takes any symbol list, and `--genes TP53,EGFR,...` an inline list. Add `--full`
to run all 3,388 shards (~2 h for the panel). For **all ~62k genes**, drop the
gene restriction entirely — that engages the heavier `unnest` map-reduce.

For the full dataset, `pipeline/build_aggregates.py` provides a parallel
map-reduce engine (`--full`, `--workers`); see its module docstring and
[`pipeline/`](pipeline/) for details. Per-component docs live in
[`api/README.md`](api/README.md) and
[`pipeline/sample_data/SAMPLE_README.md`](pipeline/sample_data/SAMPLE_README.md).

## Tests

```bash
pytest api/tests        # runs against pipeline/out; skips cleanly if absent
```

## License

MIT © Broad Institute — see [`LICENSE`](LICENSE).
