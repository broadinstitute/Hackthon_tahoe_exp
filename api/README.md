# Gene Explorer Query API (Member B)

FastAPI service that serves gene search + per-gene expression aggregates to the
frontend, reading the precompute store built by the pipeline (Member A). This
is the "Query API" layer in `project.md` §2, implementing the §4 contract.

## Endpoints

| Method & path           | Returns                                                        |
| ----------------------- | -------------------------------------------------------------- |
| `GET /api/search?q=`    | `{"matches": [{"symbol", "name"}, ...]}` — autocomplete        |
| `GET /api/gene/{symbol}`| `{gene, ensembl_id, heatmap{cell_lines, perturbations, mean}, violin[...]}` |
| `GET /api/health`       | liveness + which store/build is being served                   |

`/api/gene/{symbol}` payload (§4 contract, plus a few pass-through extras):

```jsonc
{
  "gene": "TAGLN",
  "ensembl_id": "ENSG00000149591",
  "n_groups": 30,
  "heatmap": {
    "cell_lines":    ["BT-474", "C-33 A", ...],          // matrix rows
    "perturbations": ["Acetohexamide@5.0uM", ...],        // matrix columns
    "mean":          [[null, 0.98, ...], ...]             // mean[row][col], null if group unobserved
  },
  "violin": [
    { "cell_line": "BT-474", "perturbation": "Filgotinib@5.0uM",
      "n": 1, "deciles": [...11 values...], "pct_expressing": 1.0, "mean": 1.36 }
    // one entry per observed cell_line × perturbation group
  ]
}
```

Unknown gene → `404`. Empty `q` → `422` (validation).

## Run

From the repo root:

```bash
pip install -r api/requirements.txt

# 1. Have a built store. Either unpack a teammate's tarball to pipeline/out,
#    or build the bundled sample:
python pipeline/build_aggregates.py --sample      # -> pipeline/out/

# 2. Serve (defaults to reading ../pipeline/out):
uvicorn api.main:app --reload --port 8000
# or point at any store:
STORE_DIR=/path/to/pipeline/out uvicorn api.main:app --port 8000
```

Interactive docs: <http://localhost:8000/docs>.

## How it reads the store

- `gene_index.parquet` is loaded fully into memory once at startup → `/api/search`
  is a pure in-memory scan (no disk per keystroke), prefix matches ranked first.
- `/api/gene/{symbol}` reads only that gene's Hive partition
  (`aggregates/gene_symbol=<SYM>/*.parquet`) via DuckDB — tens to low-thousands
  of rows — then pivots into the heatmap matrix + violin list. Results are
  LRU-cached (the store is immutable for a given build).

## Test

```bash
pytest api/tests        # runs against pipeline/out; skips cleanly if absent
```
