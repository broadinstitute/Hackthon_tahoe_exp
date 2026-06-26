#!/usr/bin/env python3
"""Tahoe-100M precompute pipeline (Member A).

Streams the Tahoe-100M per-cell expression Parquet shards once and computes,
for every (gene x cell_line x perturbation) group, the summary statistics the
Gene Explorer needs:

    n_cells          group size (all cells, incl. non-expressing)
    pct_expressing   fraction of cells with count > 0
    mean             mean normalized expression (incl. implicit zeros)
    deciles          11 values q0,q10,..,q100 -> violin silhouette

Output is Parquet partitioned by gene symbol (one small partition per gene),
plus a gene index for search autocomplete. See pipeline/README.md and the
output schema documented there; this maps directly onto the API contract in
project.md section 4.

Two engines
-----------
``single``    The original single-DuckDB-query path: one query over all shards
              using ``quantile_cont`` (exact deciles). Simple, but ``quantile_cont``
              is a *holistic* aggregate that materialises every value per group,
              so it does not scale to the full 314 GB / ~95M-cell dataset, and
              ``--gene-batches`` only trades that OOM for an N-pass re-scan.

``parallel``  A multiprocessing map-reduce (``--workers`` flag). Each worker
              reads a chunk of shards once, emits a *mergeable* fixed-bin
              histogram of the expressing-cell distribution plus running sums,
              and writes those partials to disk. A final DuckDB reduce sums the
              partials (out-of-core, bounded memory) and the deciles are
              reconstructed from the merged histogram with vectorised NumPy.
              Single pass over the data, bounded per-worker memory, all cores,
              and checkpointable partials. This is the path for ``--full``.

The histogram makes deciles approximate at the bin resolution (``--bins``,
default 64) rather than exact -- squarely within the "faithful violin
silhouette, not an exact KDE" trade-off documented in project.md section 3.
``mean``, ``pct_expressing`` and the cell counts remain exact in both engines.

The raw data
------------
Each row of expression_data/*.parquet is one cell, stored sparsely:
  genes        list<int64>   token ids of expressed genes (+ a leading special
                             token id 1 with sentinel value -2, which we drop)
  expressions  list<float>   raw UMI counts, parallel to `genes`
  cell_line_id string        Cellosaurus id (CVCL_*), mapped to a cell name
  drug         string        perturbation compound
  sample       string        joins sample_metadata for the concentration/dose
A gene appears in a cell's list iff its count >= 1, so pct_expressing is
n_expressing / n_cells, and mean/deciles fold in the (n_cells - n_expressing)
implicit zeros.

Normalization: per-cell counts-per-10k then log1p -- the single-cell standard,
making expression comparable across cells of differing sequencing depth.

Estimating full-run cost
-------------------------
Run the parallel engine on a handful of real shards and read the per-phase
timing + the projected full-run estimate it logs, e.g.

    python pipeline/build_aggregates.py --shards 32 --shard-stride 100 --workers 8

`--shard-stride` spreads the sampled shards across the dataset (each shard
covers only a few drugs) so the per-shard group cardinality -- and therefore
the throughput -- is representative of the full run.
"""
from __future__ import annotations

import argparse
import ast
import json
import math
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds

DEFAULT_DATA_DIR = "/home/jovyan/organization/raw/public-datasets/tahoe_100m"
SAMPLE_DATA_DIR = Path(__file__).parent / "sample_data"
N_SHARDS_TOTAL = 3388
SPECIAL_TOKEN_MIN = 3      # real genes have token_id >= 3; 1 is a CLS sentinel
NORM_TARGET = 10_000.0     # counts-per-10k before log1p
N_DECILES = 11             # q0, q10, ..., q100
DECILE_PROBS = np.linspace(0.0, 1.0, N_DECILES)
# Histogram (parallel engine) covers log1p(counts-per-10k). The theoretical max
# is a cell expressing a single gene: ln(1 + NORM_TARGET).
HIST_LO = 0.0
HIST_HI = float(np.log1p(NORM_TARGET))
DEFAULT_BINS = 64
# One Hive partition per gene symbol; the source has ~62k genes, far above
# pyarrow's default max_partitions of 1024. A single batch (e.g. the sample,
# one shard, no gene-batching) can touch thousands of genes at once.
MAX_PARTITIONS = 100_000


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def system_gb() -> float:
    """Total physical RAM in GiB (best-effort; falls back to a safe guess)."""
    try:
        return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 2**30
    except (ValueError, OSError, AttributeError):
        return 16.0


# --------------------------------------------------------------------------- #
# Lookup tables (small)                                                       #
# --------------------------------------------------------------------------- #
def parse_conc(literal: str | None):
    """Parse a drugname_drugconc cell, e.g. "[('Infigratinib', 0.05, 'uM')]".

    Returns (concentration: float|None, unit: str|None).
    """
    if not literal or not isinstance(literal, str):
        return None, None
    try:
        parsed = ast.literal_eval(literal)
    except (ValueError, SyntaxError):
        return None, None
    if not parsed:
        return None, None
    first = parsed[0]
    if not isinstance(first, (list, tuple)) or len(first) < 3:
        return None, None
    _, conc, unit = first[0], first[1], first[2]
    try:
        conc = float(conc)
    except (TypeError, ValueError):
        conc = None
    return conc, (unit if isinstance(unit, str) else None)


def build_lookups(con: duckdb.DuckDBPyConnection, data_dir: Path) -> None:
    """Create gene_map, cell_map and sample_conc tables in the connection."""
    gene_pq = data_dir / "gene_metadata" / "gene_metadata.parquet"
    cell_pq = data_dir / "cell_line_metadata" / "cell_line_metadata.parquet"
    sample_pq = data_dir / "sample_metadata" / "sample_metadata.parquet"

    con.execute(
        """
        CREATE OR REPLACE TABLE gene_map AS
        SELECT token_id, gene_symbol, ensembl_id
        FROM read_parquet(?)
        WHERE token_id >= ?
        """,
        [str(gene_pq), SPECIAL_TOKEN_MIN],
    )
    # cell_line_metadata has one row per cell; collapse to one row per line.
    con.execute(
        """
        CREATE OR REPLACE TABLE cell_map AS
        SELECT Cell_ID_Cellosaur AS cvcl,
               any_value(cell_name) AS cell_name,
               any_value(Organ)     AS organ
        FROM read_parquet(?)
        WHERE Cell_ID_Cellosaur IS NOT NULL
        GROUP BY Cell_ID_Cellosaur
        """,
        [str(cell_pq)],
    )
    # Parse the python-literal concentration in pandas, then register.
    sm = pd.read_parquet(sample_pq, columns=["sample", "drugname_drugconc"])
    concs = sm["drugname_drugconc"].map(parse_conc)
    sm["conc"] = [c for c, _ in concs]
    sm["unit"] = [u for _, u in concs]
    sample_conc = sm[["sample", "conc", "unit"]].copy()
    con.register("sample_conc_df", sample_conc)
    con.execute("CREATE OR REPLACE TABLE sample_conc AS SELECT * FROM sample_conc_df")
    con.unregister("sample_conc_df")

    n_genes = con.sql("SELECT count(*) FROM gene_map").fetchone()[0]
    n_cells = con.sql("SELECT count(*) FROM cell_map").fetchone()[0]
    n_samp = con.sql("SELECT count(*) FROM sample_conc WHERE conc IS NOT NULL").fetchone()[0]
    log(f"lookups: {n_genes} genes, {n_cells} cell lines, {n_samp} samples with dose")


def persist_lookups(con: duckdb.DuckDBPyConnection, lookups_dir: Path) -> dict:
    """Write the in-connection lookup tables to Parquet for worker processes.

    Workers run in separate processes and cannot see this connection's tables,
    so the parallel engine hands them small Parquet snapshots instead.
    """
    lookups_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "gene_map": lookups_dir / "gene_map.parquet",
        "cell_map": lookups_dir / "cell_map.parquet",
        "sample_conc": lookups_dir / "sample_conc.parquet",
    }
    for tbl, path in paths.items():
        con.execute(f"COPY {tbl} TO '{path}' (FORMAT PARQUET)")
    return {k: str(v) for k, v in paths.items()}


# --------------------------------------------------------------------------- #
# Single-query engine (exact quantiles; does not scale to --full)             #
# --------------------------------------------------------------------------- #
def aggregate_query(shard_files: list[str], token_lo: int | None, token_hi: int | None,
                    grid_probs: list[float]) -> tuple[str, list]:
    """Build the DuckDB aggregation SQL + params for the single-query engine.

    Groups every expressed (token, cell_line, perturbation) and returns, per
    group: n_expr, sum_norm, n_total (group size incl. non-expressing) and a
    fine quantile grid over the expressing-cell normalized values.
    """
    token_filter = ""
    params: list = [shard_files]
    if token_lo is not None:
        token_filter = "AND token >= ? AND token < ?"
    sql = f"""
    WITH base AS (
        SELECT cell_line_id, drug, sample, genes, expressions,
               list_sum(list_filter(expressions, x -> x > 0)) AS cell_total
        FROM read_parquet(?)
    ),
    basef AS (SELECT * FROM base WHERE cell_total > 0),
    basec AS (
        SELECT b.cell_line_id, b.drug, b.genes, b.expressions, b.cell_total,
               sc.conc, sc.unit,
               CASE WHEN sc.conc IS NULL THEN b.drug
                    ELSE b.drug || '@' || CAST(sc.conc AS VARCHAR) || COALESCE(sc.unit, '')
               END AS pert
        FROM basef b
        LEFT JOIN sample_conc sc USING (sample)
    ),
    sizes AS (
        SELECT cell_line_id, pert, count(*) AS n_total
        FROM basec GROUP BY cell_line_id, pert
    ),
    expl AS (
        SELECT unnest(genes) AS token, unnest(expressions) AS expr,
               cell_line_id, drug, conc, unit, pert, cell_total
        FROM basec
    ),
    exf AS (
        SELECT token, cell_line_id, drug, conc, unit, pert,
               ln(1 + expr / cell_total * {NORM_TARGET}) AS norm
        FROM expl
        WHERE expr > 0 AND token >= {SPECIAL_TOKEN_MIN} {token_filter}
    ),
    agg AS (
        SELECT token, cell_line_id, drug, conc, unit, pert,
               count(*)            AS n_expr,
               sum(norm)           AS sum_norm,
               quantile_cont(norm, {grid_probs}) AS qgrid
        FROM exf
        GROUP BY token, cell_line_id, drug, conc, unit, pert
    )
    SELECT a.token, a.cell_line_id, a.drug, a.conc, a.unit, a.pert,
           a.n_expr, a.sum_norm, a.qgrid,
           s.n_total,
           g.gene_symbol, g.ensembl_id,
           COALESCE(c.cell_name, a.cell_line_id) AS cell_line, c.organ
    FROM agg a
    JOIN sizes s USING (cell_line_id, pert)
    JOIN gene_map g ON g.token_id = a.token
    LEFT JOIN cell_map c ON c.cvcl = a.cell_line_id
    """
    if token_lo is not None:
        # token filter params sit inside exf, before the trailing SELECT params.
        # read_parquet(?) is the only other placeholder, already first.
        params.append(token_lo)
        params.append(token_hi)
    return sql, params


# --------------------------------------------------------------------------- #
# Parallel map-reduce engine (mergeable histogram; scales to --full)          #
# --------------------------------------------------------------------------- #
def basec_with() -> str:
    """A WITH clause exposing per-cell rows + resolved perturbation as `basec`.

    Used as a streaming CTE (NOT a materialised table) so DuckDB pipelines the
    scan -> unnest -> aggregate and prunes the heavy genes/expressions columns
    when a downstream query (e.g. sizes) doesn't read them. Reads parquet once
    via a single `?` placeholder.
    """
    return f"""
    WITH base AS (
        SELECT cell_line_id, drug, sample, genes, expressions,
               list_sum(list_filter(expressions, x -> x > 0)) AS cell_total
        FROM read_parquet(?)
    ),
    basef AS (SELECT * FROM base WHERE cell_total > 0),
    basec AS (
        SELECT b.cell_line_id, b.drug, b.genes, b.expressions, b.cell_total,
               sc.conc, sc.unit,
               CASE WHEN sc.conc IS NULL THEN b.drug
                    ELSE b.drug || '@' || CAST(sc.conc AS VARCHAR) || COALESCE(sc.unit, '')
               END AS pert
        FROM basef b
        LEFT JOIN sample_conc sc USING (sample)
    )"""


def partial_hist_query(n_bins: int, binw: float) -> str:
    """Per-(group, bin) histogram + running sums over the expressing cells.

    Mergeable: across shards, the same (token, cell_line, pert, bin) key sums
    its `cnt` and `sum_norm`. `floor(norm/binw)` is clamped into [0, n_bins-1].
    Streams from the `basec` CTE -- no full materialisation of the sparse lists.
    """
    return basec_with() + f""",
    expl AS (
        SELECT unnest(genes) AS token, unnest(expressions) AS expr,
               cell_line_id, drug, conc, unit, pert, cell_total
        FROM basec
    ),
    exf AS (
        SELECT token, cell_line_id, drug, conc, unit, pert,
               ln(1 + expr / cell_total * {NORM_TARGET}) AS norm
        FROM expl
        WHERE expr > 0 AND token >= {SPECIAL_TOKEN_MIN}
    )
    SELECT token, cell_line_id, drug, conc, unit, pert,
           least({n_bins - 1}, floor(norm / {binw!r}))::INT AS bin,
           count(*)  AS cnt,
           sum(norm) AS sum_norm
    FROM exf
    GROUP BY token, cell_line_id, drug, conc, unit, pert, bin"""


def partial_sizes_query() -> str:
    """Per-(cell_line, pert) group size. Reads only the light columns (the
    optimiser prunes genes/expressions from the `basec` CTE)."""
    return basec_with() + """
    SELECT cell_line_id, pert, count(*) AS n_total
    FROM basec GROUP BY cell_line_id, pert"""


# --- targeted few-genes path (no unnest) ----------------------------------- #
# For a handful of genes we never explode the full per-cell gene list. Instead
# we pull each target gene's count directly with list_position/list_extract
# (NULL when the gene is absent in that cell -> treated as not expressing). This
# turns ~3000 unnested rows per cell into K cheap lookups, so even a full-shard
# scan is light and memory-trivial. The per-cell `basec` then holds no lists, so
# it is safe to MATERIALISE once per shard and reference K+1 times (sizes + one
# branch per gene) without re-reading the parquet.
def targeted_basec_query(tokens: list[int]) -> str:
    """CREATE-TABLE body: per-cell row with one expression column per token."""
    ext = ",\n".join(
        f"list_extract(expressions, list_position(genes, {t})) AS e{i}"
        for i, t in enumerate(tokens)
    )
    return f"""
    WITH base AS (
        SELECT cell_line_id, drug, sample,
               list_sum(list_filter(expressions, x -> x > 0)) AS cell_total,
               {ext}
        FROM read_parquet(?)
    ),
    basef AS (SELECT * FROM base WHERE cell_total > 0)
    SELECT b.* EXCLUDE (sample), sc.conc, sc.unit,
           CASE WHEN sc.conc IS NULL THEN b.drug
                ELSE b.drug || '@' || CAST(sc.conc AS VARCHAR) || COALESCE(sc.unit, '')
           END AS pert
    FROM basef b
    LEFT JOIN sample_conc sc USING (sample)"""


def targeted_hist_query(tokens: list[int], n_bins: int, binw: float) -> str:
    """Histogram partial from a materialised targeted `basec` (no parquet read)."""
    unions = "\nUNION ALL\n".join(
        f"SELECT cell_line_id, drug, conc, unit, pert, {t} AS token, "
        f"e{i} AS expr, cell_total FROM basec WHERE e{i} > 0"
        for i, t in enumerate(tokens)
    )
    return f"""
    WITH vals AS (
        {unions}
    ),
    exf AS (
        SELECT token, cell_line_id, drug, conc, unit, pert,
               ln(1 + expr / cell_total * {NORM_TARGET}) AS norm
        FROM vals
    )
    SELECT token, cell_line_id, drug, conc, unit, pert,
           least({n_bins - 1}, floor(norm / {binw!r}))::INT AS bin,
           count(*)  AS cnt,
           sum(norm) AS sum_norm
    FROM exf
    GROUP BY token, cell_line_id, drug, conc, unit, pert, bin"""


def map_task(args: dict) -> dict:
    """Worker: aggregate one chunk of shards into partial Parquet files.

    Runs in a separate process. Streams each shard once (twice for the cheap
    sizes pass) and writes the histogram + group-size partials. Insertion-order
    preservation is off and spill is enabled, so memory stays bounded even when
    a chunk's expressing rows far exceed the worker's memory_limit.
    """
    chunk = args["chunk"]
    idx = args["idx"]
    con = duckdb.connect()
    con.execute(f"PRAGMA threads={args['threads']}")
    con.execute("PRAGMA preserve_insertion_order=false")
    if args["memory"]:
        con.execute(f"PRAGMA memory_limit='{args['memory']}'")
    if args["temp"]:
        con.execute(f"PRAGMA temp_directory='{args['temp']}'")

    con.execute(
        "CREATE TEMP TABLE sample_conc AS SELECT * FROM read_parquet(?)",
        [args["sample_conc_pq"]],
    )

    hist_pq = f"{args['partials_dir']}/hist-{idx:05d}.parquet"
    sizes_pq = f"{args['partials_dir']}/sizes-{idx:05d}.parquet"
    tokens = args.get("tokens")
    if tokens:
        # Targeted few-genes path: materialise the (list-free) per-cell table
        # once, then derive both partials from it without re-reading parquet.
        con.execute(
            f"CREATE TEMP TABLE basec AS {targeted_basec_query(tokens)}", [chunk]
        )
        con.execute(
            f"COPY ({targeted_hist_query(tokens, args['n_bins'], args['binw'])}) "
            f"TO '{hist_pq}' (FORMAT PARQUET)"
        )
        con.execute(
            "COPY (SELECT cell_line_id, pert, count(*) AS n_total "
            "FROM basec GROUP BY cell_line_id, pert) "
            f"TO '{sizes_pq}' (FORMAT PARQUET)"
        )
    else:
        con.execute(
            f"COPY ({partial_hist_query(args['n_bins'], args['binw'])}) "
            f"TO '{hist_pq}' (FORMAT PARQUET)",
            [chunk],
        )
        con.execute(
            f"COPY ({partial_sizes_query()}) TO '{sizes_pq}' (FORMAT PARQUET)",
            [chunk],
        )
    n_cells = con.sql(f"SELECT sum(n_total) FROM read_parquet('{sizes_pq}')").fetchone()[0]
    con.close()
    return {"idx": idx, "n_shards": len(chunk), "n_cells": int(n_cells or 0)}


def reduce_query(hist_glob: str) -> str:
    """Merge partial histograms into one row per group, with bins+counts lists."""
    return """
    WITH merged AS (
        SELECT token, cell_line_id, drug, conc, unit, pert, bin,
               sum(cnt)      AS cnt,
               sum(sum_norm) AS sum_norm
        FROM read_parquet(?)
        GROUP BY token, cell_line_id, drug, conc, unit, pert, bin
    ),
    grp AS (
        SELECT token, cell_line_id, drug, conc, unit, pert,
               sum(cnt)      AS n_expr,
               sum(sum_norm) AS sum_norm,
               list(bin)     AS bins,
               list(cnt)     AS cnts
        FROM merged
        GROUP BY token, cell_line_id, drug, conc, unit, pert
    )
    SELECT g.token, g.cell_line_id, g.drug, g.conc, g.unit, g.pert,
           g.n_expr, g.sum_norm, g.bins, g.cnts,
           s.n_total,
           gm.gene_symbol, gm.ensembl_id,
           COALESCE(cm.cell_name, g.cell_line_id) AS cell_line, cm.organ
    FROM grp g
    JOIN sizes s USING (cell_line_id, pert)
    JOIN gene_map gm ON gm.token_id = g.token
    LEFT JOIN cell_map cm ON cm.cvcl = g.cell_line_id
    ORDER BY gene_symbol
    """


# --------------------------------------------------------------------------- #
# Decile reconstruction (vectorised)                                          #
# --------------------------------------------------------------------------- #
def fold_zeros_vectorized(qgrid: np.ndarray, grid_probs: np.ndarray,
                          pct: np.ndarray) -> np.ndarray:
    """Single engine: fold the zero mass into expressing-cell quantile grids.

    `qgrid` is (G, len(grid_probs)) -- values of each group's expressing
    distribution at `grid_probs`. With zero fraction z = 1 - pct, the full
    quantile at prob p is 0 for p <= z, else the expressing quantile at
    (p - z)/pct. Vectorised over all G groups and all 11 deciles at once.
    """
    g = qgrid.shape[0]
    out = np.zeros((g, N_DECILES), dtype=np.float64)
    z = 1.0 - pct
    safe_pct = np.where(pct > 0, pct, 1.0)
    rows = np.arange(g)
    for j, p in enumerate(DECILE_PROBS):
        qc = np.clip((p - z) / safe_pct, 0.0, 1.0)
        # Interpolate qc on the shared, monotonic grid_probs (per-row y-values).
        idx = np.clip(np.searchsorted(grid_probs, qc, side="right"), 1, len(grid_probs) - 1)
        x0, x1 = grid_probs[idx - 1], grid_probs[idx]
        y0, y1 = qgrid[rows, idx - 1], qgrid[rows, idx]
        denom = x1 - x0
        t = np.where(denom > 0, (qc - x0) / np.where(denom > 0, denom, 1.0), 0.0)
        val = y0 + t * (y1 - y0)
        val[(p <= z) | (pct <= 0)] = 0.0
        out[:, j] = val
    return out


def hist_deciles(cum: np.ndarray, n_expr: np.ndarray, pct: np.ndarray,
                 edges: np.ndarray) -> np.ndarray:
    """Parallel engine: full-distribution deciles from merged histograms.

    `cum` is (G, n_bins) cumulative expressing-cell counts; `edges` the
    (n_bins+1) bin edges. Cumulative fraction reaches 1.0 at the last occupied
    bin, so the inverse-CDF lookup yields that bin's upper edge for q100. Zeros
    are folded in exactly as in `fold_zeros_vectorized`.
    """
    g, n_bins = cum.shape
    out = np.zeros((g, N_DECILES), dtype=np.float64)
    z = 1.0 - pct
    safe_pct = np.where(pct > 0, pct, 1.0)
    safe_expr = np.where(n_expr > 0, n_expr, 1.0)
    frac = cum / safe_expr[:, None]      # fraction at each bin's upper edge
    upper = edges[1:]                     # value at each bin's upper edge
    rows = np.arange(g)
    for j, p in enumerate(DECILE_PROBS):
        qc = np.clip((p - z) / safe_pct, 0.0, 1.0)
        ge = frac >= qc[:, None]
        idx = ge.argmax(axis=1)           # first bin whose cum fraction >= qc
        idx[~ge.any(axis=1)] = n_bins - 1
        prev_frac = np.where(idx > 0, frac[rows, np.maximum(idx - 1, 0)], 0.0)
        prev_x = np.where(idx > 0, upper[np.maximum(idx - 1, 0)], edges[0])
        cur_frac = frac[rows, idx]
        denom = cur_frac - prev_frac
        t = np.where(denom > 0, (qc - prev_frac) / np.where(denom > 0, denom, 1.0), 0.0)
        val = prev_x + t * (upper[idx] - prev_x)
        val[(p <= z) | (pct <= 0)] = 0.0
        out[:, j] = val
    return out


# --------------------------------------------------------------------------- #
# Output assembly (shared by both engines)                                    #
# --------------------------------------------------------------------------- #
def assemble_output(df: pd.DataFrame, pct: np.ndarray, mean: np.ndarray,
                    deciles: np.ndarray) -> pa.Table:
    """Build the API-contract output table (one row per group)."""
    out = pd.DataFrame(
        {
            "gene_symbol": df["gene_symbol"].astype(str),
            "ensembl_id": df["ensembl_id"].astype(str),
            "token_id": df["token"].astype("int64"),
            "cell_line": df["cell_line"].astype(str),
            "organ": df["organ"].astype("object"),
            "perturbation": df["pert"].astype(str),
            "drug": df["drug"].astype(str),
            "concentration": df["conc"].astype("float64"),
            "conc_unit": df["unit"].astype("object"),
            "n_cells": df["n_total"].astype("int64"),
            "n_expressing": df["n_expr"].astype("int64"),
            "pct_expressing": np.round(pct, 6),
            "mean": np.round(mean, 6),
        }
    )
    table = pa.Table.from_pandas(out, preserve_index=False)
    return table.append_column(
        "deciles", pa.array(deciles.tolist(), type=pa.list_(pa.float64()))
    )


def finalize_batch(df: pd.DataFrame, grid_probs: np.ndarray) -> pa.Table:
    """Single engine: raw aggregate rows -> output schema (vectorised)."""
    n_total = df["n_total"].to_numpy(dtype=np.float64)
    n_expr = df["n_expr"].to_numpy(dtype=np.float64)
    pct = np.where(n_total > 0, n_expr / n_total, 0.0)
    mean = np.where(n_total > 0, df["sum_norm"].to_numpy(dtype=np.float64) / n_total, 0.0)
    qgrid = np.asarray(df["qgrid"].tolist(), dtype=np.float64)
    deciles = fold_zeros_vectorized(qgrid, grid_probs, pct)
    return assemble_output(df, pct, mean, deciles)


def finalize_batch_hist(df: pd.DataFrame, n_bins: int, edges: np.ndarray) -> pa.Table:
    """Parallel engine: merged-histogram rows -> output schema (vectorised)."""
    n_total = df["n_total"].to_numpy(dtype=np.float64)
    n_expr = df["n_expr"].to_numpy(dtype=np.float64)
    pct = np.where(n_total > 0, n_expr / n_total, 0.0)
    mean = np.where(n_total > 0, df["sum_norm"].to_numpy(dtype=np.float64) / n_total, 0.0)

    # Scatter the ragged (bins, cnts) lists into a dense (G, n_bins) matrix.
    g = len(df)
    counts = np.zeros((g, n_bins), dtype=np.float64)
    lengths = df["bins"].map(len).to_numpy()
    if lengths.sum() > 0:
        row_idx = np.repeat(np.arange(g), lengths)
        bin_idx = np.concatenate([np.asarray(b, dtype=np.int64) for b in df["bins"]])
        cnt_val = np.concatenate([np.asarray(c, dtype=np.float64) for c in df["cnts"]])
        counts[row_idx, bin_idx] = cnt_val
    cum = np.cumsum(counts, axis=1)
    deciles = hist_deciles(cum, n_expr, pct, edges)
    return assemble_output(df, pct, mean, deciles)


# --------------------------------------------------------------------------- #
# Shard / token helpers                                                       #
# --------------------------------------------------------------------------- #
def select_shards(data_dir: Path, n_shards: int | None, full: bool,
                  stride: int = 1) -> list[str]:
    expr_dir = data_dir / "expression_data"
    files = sorted(str(p) for p in expr_dir.glob("train-*.parquet"))
    if not files:
        sys.exit(f"No expression shards found in {expr_dir}")
    if full or n_shards is None or (stride <= 1 and n_shards >= len(files)):
        return files
    if stride > 1:
        files = files[::stride]
    return files[:n_shards]


def chunk_shards(shards: list[str], size: int) -> list[list[str]]:
    return [shards[i:i + size] for i in range(0, len(shards), size)]


def token_buckets(con: duckdb.DuckDBPyConnection, n_buckets: int):
    """Yield (lo, hi) half-open token ranges covering all gene tokens."""
    lo_tok, hi_tok = con.sql("SELECT min(token_id), max(token_id) FROM gene_map").fetchone()
    if n_buckets <= 1:
        yield (None, None)
        return
    edges = np.linspace(lo_tok, hi_tok + 1, n_buckets + 1).astype(int)
    for i in range(n_buckets):
        yield (int(edges[i]), int(edges[i + 1]))


# --------------------------------------------------------------------------- #
# Engine drivers                                                              #
# --------------------------------------------------------------------------- #
def run_single(con, shards, args, grid_probs, grid_list, agg_dir):
    """Original single-query engine. Exact deciles; bounded by memory only via
    --gene-batches (each pass re-scans all shards)."""
    seen_symbols: set[str] = set()
    total_rows = 0
    for bi, (lo, hi) in enumerate(token_buckets(con, args.gene_batches)):
        label = "all tokens" if lo is None else f"tokens [{lo}, {hi})"
        log(f"pass {bi + 1}/{args.gene_batches}: aggregating {label} ...")
        sql, params = aggregate_query(shards, lo, hi, grid_list)
        reader = con.execute(sql, params).fetch_record_batch(rows_per_batch=200_000)
        wrote_any = False
        for rb in reader:
            df = rb.to_pandas()
            if df.empty:
                continue
            table = finalize_batch(df, grid_probs)
            seen_symbols.update(table.column("gene_symbol").to_pylist())
            total_rows += table.num_rows
            ds.write_dataset(
                table, base_dir=str(agg_dir), format="parquet",
                partitioning=["gene_symbol"], partitioning_flavor="hive",
                existing_data_behavior="overwrite_or_ignore",
                basename_template=f"part-b{bi}-{{i}}.parquet",
                max_partitions=MAX_PARTITIONS,
            )
            wrote_any = True
        if not wrote_any:
            log(f"  (no expressed rows in {label})")
    return seen_symbols, total_rows, {}


def run_parallel(con, shards, args, agg_dir, out_dir, tokens=None):
    """Map-reduce engine: parallel partial histograms -> out-of-core reduce.

    When `tokens` is given, workers use the targeted few-genes path (no unnest).
    """
    n_bins = args.bins
    edges = np.linspace(HIST_LO, HIST_HI, n_bins + 1)
    binw = HIST_HI / n_bins

    partials_dir = out_dir / "partials"
    lookups_dir = out_dir / "lookups"
    tmp_dir = out_dir / ".duckdb_tmp"
    for d in (partials_dir, lookups_dir, tmp_dir):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    lookup_paths = persist_lookups(con, lookups_dir)

    total_threads = os.cpu_count() or 2
    worker_threads = max(1, total_threads // args.workers)
    if args.worker_memory:
        worker_mem = args.worker_memory
    else:
        worker_mem = f"{max(2, int(0.7 * system_gb() / args.workers))}GB"

    chunks = chunk_shards(shards, args.shards_per_task)
    tasks = [
        {
            "chunk": chunk, "idx": i,
            "n_bins": n_bins, "binw": binw,
            "sample_conc_pq": lookup_paths["sample_conc"],
            "partials_dir": str(partials_dir),
            "threads": worker_threads, "memory": worker_mem,
            "temp": str(tmp_dir), "tokens": tokens,
        }
        for i, chunk in enumerate(chunks)
    ]
    mode = f"targeted {len(tokens)} tokens (no unnest)" if tokens else "all genes"
    log(f"MAP [{mode}]: {len(shards)} shards in {len(tasks)} tasks "
        f"({args.shards_per_task}/task) x {args.workers} workers "
        f"({worker_threads} duckdb-threads, {worker_mem}/worker)")

    t_map0 = time.time()
    done_shards = done_cells = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(map_task, t) for t in tasks]
        for k, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            done_shards += r["n_shards"]
            done_cells += r["n_cells"]
            if k % max(1, len(tasks) // 10) == 0 or k == len(tasks):
                log(f"  map {k}/{len(tasks)} tasks "
                    f"({done_shards} shards, {done_cells:,} cells)")
    t_map = time.time() - t_map0
    log(f"MAP done: {done_shards} shards, {done_cells:,} cells in {t_map:.1f}s "
        f"({done_shards / t_map:.2f} shards/s)")

    # -- reduce ----------------------------------------------------------- #
    reduce_mem = args.memory_limit or f"{max(4, int(0.7 * system_gb()))}GB"
    con.execute(f"PRAGMA memory_limit='{reduce_mem}'")
    con.execute(f"PRAGMA temp_directory='{tmp_dir}'")
    con.execute("PRAGMA preserve_insertion_order=false")
    if args.threads:
        con.execute(f"PRAGMA threads={args.threads}")
    con.execute(
        "CREATE OR REPLACE TABLE sizes AS "
        "SELECT cell_line_id, pert, sum(n_total) AS n_total "
        f"FROM read_parquet('{partials_dir}/sizes-*.parquet') "
        "GROUP BY cell_line_id, pert"
    )
    log(f"REDUCE: merging partials ({reduce_mem} limit) ...")
    t_red0 = time.time()
    hist_glob = f"{partials_dir}/hist-*.parquet"
    reader = con.execute(reduce_query(hist_glob), [hist_glob]).fetch_record_batch(
        rows_per_batch=200_000
    )
    seen_symbols: set[str] = set()
    total_rows = 0
    t_write_total = 0.0
    for rb in reader:
        df = rb.to_pandas()
        if df.empty:
            continue
        table = finalize_batch_hist(df, n_bins, edges)
        seen_symbols.update(table.column("gene_symbol").to_pylist())
        total_rows += table.num_rows
        tw0 = time.time()
        ds.write_dataset(
            table, base_dir=str(agg_dir), format="parquet",
            partitioning=["gene_symbol"], partitioning_flavor="hive",
            existing_data_behavior="overwrite_or_ignore",
            basename_template="part-{i}.parquet",
            max_partitions=MAX_PARTITIONS,
        )
        t_write_total += time.time() - tw0
    t_reduce = time.time() - t_red0 - t_write_total

    if not args.keep_partials:
        shutil.rmtree(partials_dir, ignore_errors=True)
        shutil.rmtree(lookups_dir, ignore_errors=True)
    shutil.rmtree(tmp_dir, ignore_errors=True)

    timings = {
        "map_seconds": round(t_map, 1),
        "reduce_seconds": round(t_reduce, 1),
        "write_seconds": round(t_write_total, 1),
        "shards_per_second": round(done_shards / t_map, 3) if t_map else None,
        "reduce_rows_per_second": round(total_rows / t_reduce, 0) if t_reduce else None,
        "n_cells": done_cells,
    }
    log(f"REDUCE done: {total_rows:,} group rows, "
        f"reduce {t_reduce:.1f}s + write {t_write_total:.1f}s")
    log_full_estimate(done_shards, done_cells, t_map, t_reduce, t_write_total,
                      total_rows, targeted=bool(tokens))
    return seen_symbols, total_rows, timings


# project.md section 3 scale sanity check: ~1.5B aggregate group rows total.
EST_FULL_GROUP_ROWS = 1.5e9


def log_full_estimate(n_shards, n_cells, t_map, t_reduce, t_write, total_rows,
                      targeted=False):
    """Project full-run cost from a partial run.

    Map is embarrassingly parallel and scales ~linearly with shards -- a
    reliable projection. Reduce+write cost tracks the number of *output group
    rows*, which SATURATES on real data (the same gene x cell_line x pert group
    recurs across shards) rather than growing with shard count. So we project
    reduce from measured rows/s against the dataset's expected total group-row
    count (project.md s3, ~1.5B), not by shard scaling -- and flag both the
    assumption and the fact that a 1-shard / bundled-sample run is structurally
    unrepresentative (near 1 cell per group => far more, far smaller groups).
    """
    if n_shards <= 0 or t_map <= 0 or n_shards >= N_SHARDS_TOTAL:
        return
    scale = N_SHARDS_TOTAL / n_shards
    proj_map = t_map * scale
    log("-" * 70)
    log(f"FULL-RUN ESTIMATE (from {n_shards} shards / {n_cells:,} cells "
        f"-> {N_SHARDS_TOTAL} shards):")
    log(f"  map           ~ {proj_map / 60:.1f} min   "
        f"(linear x{scale:.0f}; reliable)")
    if targeted:
        # Group count is bounded (~n_genes x cell_lines x perts) and already
        # near-saturated, so the reduce stays ~constant as shards grow.
        log(f"  reduce+write  ~ {(t_reduce + t_write):.1f}s (~constant; group "
            "count saturates for a few-genes build, doesn't grow with shards)")
        log(f"  TOTAL         ~ {proj_map / 60:.1f} min   (map-dominated)")
    else:
        rd_rows_per_s = total_rows / t_reduce if t_reduce > 0 else 0
        proj_rd = (EST_FULL_GROUP_ROWS / rd_rows_per_s) if rd_rows_per_s else float("nan")
        log(f"  reduce+write  ~ {(proj_rd + t_write * scale) / 60:.1f} min   "
            f"(at {rd_rows_per_s:,.0f} group-rows/s vs ~{EST_FULL_GROUP_ROWS:.1g} "
            "expected rows; depends on final group count, not shards)")
        log(f"  TOTAL (rough) ~ {(proj_map + proj_rd + t_write * scale) / 60:.1f} min")
        log("  NOTE: run on REAL shards with --shard-stride (NOT the bundled 1-shard")
        log("  sample) for a representative reduce estimate -- the sample has ~1 cell")
        log("  per group, inflating group count and reduce time by orders of magnitude.")
    log("-" * 70)


def resolve_gene_tokens(con: duckdb.DuckDBPyConnection, genes_arg: str) -> list[int]:
    """Map comma-separated gene symbols to their token_ids via gene_map.

    A symbol may resolve to more than one token; all are kept. Unknown symbols
    are reported and skipped. Exits if none resolve.
    """
    symbols = [s.strip() for s in genes_arg.split(",") if s.strip()]
    rows = con.execute(
        "SELECT gene_symbol, token_id FROM gene_map WHERE gene_symbol IN "
        f"({','.join(['?'] * len(symbols))}) ORDER BY gene_symbol, token_id",
        symbols,
    ).fetchall()
    found = {}
    for sym, tok in rows:
        found.setdefault(sym, []).append(int(tok))
    missing = [s for s in symbols if s not in found]
    if missing:
        log(f"WARNING: gene symbol(s) not in gene_metadata, skipped: {missing}")
    tokens = [t for toks in found.values() for t in toks]
    if not tokens:
        sys.exit(f"None of the requested genes resolved to a token: {symbols}")
    log(f"targeting {len(tokens)} token(s) for {len(found)} gene(s): "
        + ", ".join(f"{s}={found[s]}" for s in found))
    return tokens


# --------------------------------------------------------------------------- #
# Driver                                                                      #
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=None, type=Path,
                    help=f"Tahoe-100M dataset root (default: {DEFAULT_DATA_DIR})")
    ap.add_argument("--sample", action="store_true",
                    help=f"process the tiny bundled sample at {SAMPLE_DATA_DIR} "
                         "(all of its shards); a shortcut for "
                         "'--data-dir pipeline/sample_data --full'. Use to produce "
                         "example aggregates for API development.")
    ap.add_argument("--out", default=Path(__file__).parent / "out", type=Path,
                    help="output directory (default: pipeline/out)")
    ap.add_argument("--genes", default=None,
                    help="comma-separated gene symbols to restrict to, e.g. "
                         "'TP53,EGFR,KRAS,BRCA1'. Uses a targeted no-unnest scan "
                         "(parallel engine) that is cheap enough to run over all "
                         "shards (combine with --full) for a few-gene demo build.")
    ap.add_argument("--shards", type=int, default=4,
                    help="number of expression shards to process (default: 4, a "
                         "Day-1 sample). Ignored when --full is set.")
    ap.add_argument("--shard-stride", type=int, default=1,
                    help="pick every Nth shard before taking --shards of them, to "
                         "spread the sample across drugs for a representative "
                         "throughput estimate (default: 1). Ignored when --full.")
    ap.add_argument("--full", action="store_true",
                    help=f"process all {N_SHARDS_TOTAL} shards (~95M cells; slow, heavy)")
    ap.add_argument("--engine", choices=["auto", "single", "parallel"], default="auto",
                    help="aggregation engine (default: auto -> single for --sample, "
                         "else parallel). 'single' = exact quantiles, one query, "
                         "does not scale; 'parallel' = multiprocessing map-reduce "
                         "with histogram deciles, for --full.")
    ap.add_argument("--workers", type=int, default=0,
                    help="parallel engine: worker processes (default: 0 = auto = "
                         "min(8, cpus-2)).")
    ap.add_argument("--shards-per-task", type=int, default=4,
                    help="parallel engine: shards aggregated per worker task "
                         "(default: 4). Larger = fewer/larger partials, more "
                         "per-worker memory.")
    ap.add_argument("--bins", type=int, default=DEFAULT_BINS,
                    help=f"parallel engine: histogram bins over log1p(cp10k) in "
                         f"[0, {HIST_HI:.2f}] (default: {DEFAULT_BINS}). Higher = "
                         "finer deciles, larger partials.")
    ap.add_argument("--worker-memory", default=None,
                    help="parallel engine: per-worker DuckDB memory_limit, e.g. "
                         "'4GB' (default: ~0.7*RAM/workers).")
    ap.add_argument("--keep-partials", action="store_true",
                    help="parallel engine: keep the intermediate partials/ and "
                         "lookups/ dirs (for debugging / checkpointing).")
    ap.add_argument("--grid-points", type=int, default=21,
                    help="single engine: quantile-grid resolution for the "
                         "expressing distribution (default: 21).")
    ap.add_argument("--gene-batches", type=int, default=1,
                    help="single engine: split aggregation into N token-range "
                         "passes to bound memory (default: 1). Each pass rescans "
                         "shards. (Unused by the parallel engine.)")
    ap.add_argument("--threads", type=int, default=None, help="DuckDB threads")
    ap.add_argument("--memory-limit", default=None,
                    help="DuckDB memory limit for the single engine / parallel "
                         "reduce, e.g. '16GB' (enables disk spill)")
    args = ap.parse_args()

    if args.sample and args.data_dir is not None:
        ap.error("--sample and --data-dir are mutually exclusive")
    if args.sample:
        data_dir = SAMPLE_DATA_DIR.resolve()
        full = True  # the sample is small; always use all of its shards
    else:
        data_dir = (args.data_dir or Path(DEFAULT_DATA_DIR)).resolve()
        full = args.full
    out_dir = args.out.resolve()
    agg_dir = out_dir / "aggregates"

    if args.workers <= 0:
        args.workers = min(8, max(1, (os.cpu_count() or 2) - 2))

    shards = select_shards(data_dir, args.shards, full, args.shard_stride)

    engine = args.engine
    if engine == "auto":
        engine = "single" if args.sample else "parallel"
    if args.genes and engine != "parallel":
        log("note: --genes forces the parallel engine (targeted no-unnest path)")
        engine = "parallel"
    if engine == "parallel" and args.workers == 1:
        log("note: --workers 1 with parallel engine (no concurrency)")

    grid_probs = np.linspace(0.0, 1.0, args.grid_points)
    grid_list = grid_probs.tolist()

    log(f"data_dir={data_dir}")
    log(f"out={out_dir}")
    if args.sample:
        mode = "bundled sample"
    elif full:
        mode = "FULL"
    else:
        mode = "sample"
    log(f"processing {len(shards)} / {N_SHARDS_TOTAL} shards "
        f"({mode}); engine={engine}"
        + (f"; gene_batches={args.gene_batches}" if engine == "single"
           else f"; bins={args.bins}"))

    if out_dir.exists():
        shutil.rmtree(out_dir)
    agg_dir.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    if args.threads and engine == "single":
        con.execute(f"PRAGMA threads={args.threads}")
    if args.memory_limit and engine == "single":
        con.execute(f"PRAGMA memory_limit='{args.memory_limit}'")
    if engine == "single":
        con.execute("PRAGMA enable_progress_bar")

    t0 = time.time()
    build_lookups(con, data_dir)

    tokens = resolve_gene_tokens(con, args.genes) if args.genes else None

    if engine == "single":
        seen_symbols, total_rows, timings = run_single(
            con, shards, args, grid_probs, grid_list, agg_dir)
    else:
        seen_symbols, total_rows, timings = run_parallel(
            con, shards, args, agg_dir, out_dir, tokens=tokens)

    write_gene_index(con, out_dir, seen_symbols)
    write_manifest(out_dir, data_dir, shards, args, full, total_rows,
                   len(seen_symbols), time.time() - t0, engine, timings)
    con.close()
    log(f"DONE: {total_rows} group rows across {len(seen_symbols)} genes "
        f"in {time.time() - t0:.1f}s -> {agg_dir}")


def write_gene_index(con: duckdb.DuckDBPyConnection, out_dir: Path,
                     seen_symbols: set[str]) -> None:
    """Gene index for /api/search autocomplete: every gene that has data."""
    gm = con.sql("SELECT gene_symbol, ensembl_id, token_id FROM gene_map").df()
    gm = gm[gm["gene_symbol"].isin(seen_symbols)].copy()
    gm = gm.sort_values("gene_symbol").reset_index(drop=True)
    # contract uses {symbol, name}; no description in source, so name := symbol.
    gm["symbol"] = gm["gene_symbol"]
    gm["name"] = gm["gene_symbol"]
    gm[["symbol", "name", "ensembl_id", "token_id"]].to_parquet(
        out_dir / "gene_index.parquet", index=False
    )
    log(f"gene_index: {len(gm)} genes -> {out_dir / 'gene_index.parquet'}")


def write_manifest(out_dir: Path, data_dir: Path, shards: list[str], args, full: bool,
                   total_rows: int, n_genes: int, elapsed: float,
                   engine: str, timings: dict) -> None:
    manifest = {
        "data_dir": str(data_dir),
        "sample": bool(args.sample),
        "engine": engine,
        "n_shards_processed": len(shards),
        "n_shards_total": N_SHARDS_TOTAL,
        "full": bool(full),
        "normalization": "log1p(counts_per_10k)",
        "perturbation": "drug@concentration+unit",
        "special_token_min": SPECIAL_TOKEN_MIN,
        "n_deciles": N_DECILES,
        "n_group_rows": total_rows,
        "n_genes": n_genes,
        "elapsed_seconds": round(elapsed, 1),
    }
    if args.genes:
        manifest["genes"] = [s.strip() for s in args.genes.split(",") if s.strip()]
        manifest["targeted"] = True
    if engine == "single":
        manifest["deciles_method"] = "quantile_cont (exact)"
        manifest["grid_points"] = args.grid_points
        manifest["gene_batches"] = args.gene_batches
    else:
        manifest["deciles_method"] = f"histogram ({args.bins} bins, approximate)"
        manifest["bins"] = args.bins
        manifest["workers"] = args.workers
        manifest["shards_per_task"] = args.shards_per_task
        manifest["timings"] = timings
    (out_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    log(f"manifest -> {out_dir / 'run_manifest.json'}")


if __name__ == "__main__":
    main()
