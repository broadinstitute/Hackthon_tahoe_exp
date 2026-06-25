# Tahoe-100M Gene Explorer — Project Brainstorm

> A GTEx-style search & visualization website for the Tahoe-100M single-cell
> dataset, served from the Manifold platform.
> Reference UX: https://gtexportal.org/home/gene/TAGLN

**Status:** brainstorm / design draft · **Date:** 2026-06-25 · **Team:** 3

---

## 1. Goal

Build a search engine over the Tahoe-100M dataset. A user searches for a
**gene**, and the site renders, for that gene:

- a **violin plot** of its expression distribution, and
- a **heatmap** of mean expression across **cell line × perturbation**.

The data lives as large **h5ad single-cell** files (~100M cells). The core
engineering challenge is a **fast query backend** — h5ad is cell/row-oriented,
so naive per-gene queries over 100M cells are far too slow for an interactive
site.

### MVP (hackathon target)

**Single gene → both plots.** Search one gene symbol, render the violin plot
and the cell-line × perturbation heatmap. This is the equivalent of a single
GTEx gene page and is fully demoable.

Out of scope for MVP (nice-to-haves, see §8): multi-gene comparison,
cross-cell-line filtering, dose-response views, export, auth.

---

## 2. Architecture

Three layers, which map cleanly onto our three team members:

```
 h5ad (100M cells)                                            Browser
        │                                                        ▲
        ▼                                                        │ JSON
 ┌───────────────┐   Parquet/    ┌───────────────┐   HTTP   ┌────┴──────┐
 │ 1. Precompute  │── DuckDB ───▶ │ 2. Query API   │ ───────▶ │ 3. Front-  │
 │    pipeline    │    store      │    (FastAPI)   │          │    end     │
 └───────────────┘               └───────────────┘          └───────────┘
   offline, batch                   ms-latency                search + plots
```

1. **Precompute pipeline** (offline, batch): stream the h5ad once, compute
   per-group summary statistics, write to a fast columnar store.
2. **Query API**: read the store, serve gene queries + search autocomplete at
   millisecond latency.
3. **Frontend**: gene search box + gene page rendering the two plots.

---

## 3. The key decision: precompute aggregates

We **do not** query raw cells at request time. Instead, an offline job
precomputes, for each **gene × cell-line × perturbation** group:

| Field            | Use                                                |
| ---------------- | -------------------------------------------------- |
| `n_cells`        | group size / confidence                            |
| `mean`           | heatmap color (raw or z-scored)                    |
| `pct_expressing` | fraction of cells with expression > 0              |
| `deciles`        | 11 values (q0, q10 … q90, q100) → **violin shape** |

**Storage:** Parquet **partitioned by gene** (or a DuckDB table indexed on
gene). A single-gene query then touches only that gene's partition
(~50 cell lines × ~1,100 perturbations of rows) and returns in milliseconds.

**Trade-off to be aware of:** because we store aggregates (not raw cells), the
violin is reconstructed from **precomputed deciles**, not a true KDE over raw
values. Deciles give a faithful violin silhouette without storing 100M raw
values. If we later want exact KDE violins, we add a subsampled-cells path
(see §8) — this is the "hybrid" option we deferred.

**Scale sanity check:** ~50 cell lines × ~1,100 perturbations ≈ 55k groups per
gene; with ~20–30k genes that's ~1.5B aggregate rows total. Large but routine
for Parquet/DuckDB when partitioned by gene; per-query cost stays tiny.

---

## 4. The contract (freeze this on Day 1)

The single most important step for parallel work: **agree the JSON schema
first**, then everyone builds against a mock of it.

```jsonc
// GET /api/gene/{symbol}
{
  "gene": "TAGLN",
  "heatmap": {
    "cell_lines":    ["A549", "MCF7", ...],          // rows
    "perturbations": ["DMSO", "Drug_A@1uM", ...],     // columns
    "mean":          [[0.12, 0.98, ...], ...]         // [cell_line][perturbation]
  },
  "violin": [
    {
      "cell_line": "A549",
      "perturbation": "DMSO",
      "n": 1234,
      "deciles": [0.0, 0.1, 0.2, ...],                // 11 values q0..q100
      "pct_expressing": 0.81
    }
    // ... one entry per cell_line × perturbation group
  ]
}
```

```jsonc
// GET /api/search?q=TAG   → autocomplete
{ "matches": [ { "symbol": "TAGLN", "name": "transgelin" }, ... ] }
```

Once frozen, the frontend and pipeline no longer block on each other.

---

## 5. Divide & conquer — 3 members

Each member **owns one layer end-to-end** and exposes a well-defined interface.
The contract in §4 is the seam between them.

### Member A — Data / Precompute pipeline (`pipeline/`)
- Owns: streaming the h5ad → per-group aggregates; the storage format
  (DuckDB/Parquet); the gene-symbol index used for search.
- **Day-1 deliverable:** a small **sample dataset** (a handful of genes, real
  schema) so B and C are unblocked immediately.
- Done when: the full store builds reproducibly from raw h5ad via one command.

### Member B — Query API (`api/`)
- Owns: FastAPI service reading the store; `GET /api/gene/{symbol}`,
  `GET /api/search?q=`; response shaping to the §4 contract; caching.
- **Day-1 deliverable:** a **mock server** returning the §4 schema from static
  fixtures, so C can build the UI before the real data exists.
- Done when: real endpoints serve A's store at ms latency.

### Member C — Frontend (`web/`)
- Owns: search box + autocomplete, gene page, violin + heatmap rendering
  (Plotly), loading/empty/error states, layout.
- Builds against B's mock until the real API is live.
- Done when: search a gene → both plots render from the live API.

> **Why this split works:** each layer has one purpose, a clear interface, and
> can be developed and tested independently. The only shared artifact is the
> §4 contract and the sample dataset — both delivered Day 1.

---

## 6. Tech stack (one open decision)

Backend is **Python** regardless (the h5ad pipeline is Python/scanpy world).
Two options for serving + frontend:

| Option                            | Pros                                                 | Cons                                |
| --------------------------------- | ---------------------------------------------------- | ----------------------------------- |
| **FastAPI + React/Plotly.js** ⭐  | Cleanest 3-way parallel split; GTEx-level polish     | More setup; two languages           |
| **Plotly Dash (Python all-in)**   | Fastest to demo; one language; great for scientists  | Less custom UX; front/back coupled  |

**Recommendation:** FastAPI + React for the clean parallel split and polish.
**Fallback if time is tight:** Dash, all-in Python.

> ⛳ **OPEN DECISION — pick before kickoff.** This changes Member C's stack only;
> A and B are unaffected.

---

## 7. Repo conventions (3 people on `main`)

- **Branch per layer:** `data/<topic>`, `api/<topic>`, `web/<topic>`.
  Short-lived; rebase/merge often.
- **PRs:** one reviewer before merge to `main`. Keep `main` green.
- **Top-level dirs:** `pipeline/`, `api/`, `web/` so members rarely edit the
  same files → minimal merge conflicts.
- **Shared contract** (§4) lives in one place (e.g. `api/schema.json` or a
  shared `CONTRACT.md`); changes to it require a heads-up to all three.

---

## 8. Future / deferred (YAGNI for MVP)

- Exact KDE violins via a subsampled-raw-cells path (the "hybrid" store).
- Multi-gene and multi-perturbation comparison views.
- Dose-response curves; cross-cell-line filtering and ranking.
- PNG/CSV export, shareable URLs, gene metadata panel.
- Manifold platform integration / deployment specifics.

---

## 9. Suggested timeline (compress to your hackathon length)

| Phase     | A (pipeline)                       | B (API)                         | C (frontend)                         |
| --------- | ---------------------------------- | ------------------------------- | ------------------------------------ |
| **Day 1** | Sample dataset + freeze §4 schema  | Mock server from fixtures       | Page scaffold + search box           |
| **Day 2** | Full aggregate build on real h5ad  | Real `/gene` + `/search`        | Violin + heatmap vs mock             |
| **Day 3** | Tune storage/partitioning, index   | Wire to A's store; caching      | Wire to live API; polish; demo       |

---

## 10. Open questions to confirm at kickoff

1. Heatmap axes: **cell line (rows) × perturbation (columns)** — confirm
   orientation and whether color is raw mean or z-score.
2. Violins from **precomputed deciles** (no raw cells stored) — acceptable for
   MVP?
3. Stack choice from §6 (FastAPI+React vs Dash).
4. How the Manifold platform hosts/serves the final site (deployment target).
