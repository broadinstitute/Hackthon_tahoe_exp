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

Deciles: computed by sampling the expressing-cell distribution on a fine
quantile grid in DuckDB (bounded memory, single streaming pass), then folding
in the zero mass in Python. This is a faithful violin silhouette, not an exact
KDE over raw cells (an explicit MVP trade-off, see project.md section 3).
"""
from __future__ import annotations

import argparse
import ast
import json
import shutil
import sys
import time
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
# One Hive partition per gene symbol; the source has ~62k genes, far above
# pyarrow's default max_partitions of 1024. A single batch (e.g. the sample,
# one shard, no gene-batching) can touch thousands of genes at once.
MAX_PARTITIONS = 100_000


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


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


# --------------------------------------------------------------------------- #
# Core aggregation                                                            #
# --------------------------------------------------------------------------- #
def aggregate_query(shard_files: list[str], token_lo: int | None, token_hi: int | None,
                    grid_probs: list[float]) -> tuple[str, list]:
    """Build the DuckDB aggregation SQL + params.

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
           COALESCE(c.cell_name, a.cell_line_id) AS cell_line,
           c.organ
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


def finalize_deciles(qgrid: np.ndarray, grid_probs: np.ndarray, pct: float) -> np.ndarray:
    """Fold the zero mass into the expressing-cell quantile grid.

    `qgrid` are values of the expressing distribution at `grid_probs`. With a
    zero fraction z = 1 - pct, the full (zero-inclusive) quantile at prob p is 0
    for p <= z, else the expressing quantile at (p - z) / pct.
    """
    z = 1.0 - pct
    out = np.zeros(N_DECILES, dtype=np.float64)
    for i, p in enumerate(DECILE_PROBS):
        if p <= z or pct <= 0:
            out[i] = 0.0
        else:
            p_adj = (p - z) / pct
            out[i] = float(np.interp(p_adj, grid_probs, qgrid))
    return out


def finalize_batch(df: pd.DataFrame, grid_probs: np.ndarray) -> pa.Table:
    """Turn raw aggregate rows into the output schema (one row per group)."""
    n_total = df["n_total"].to_numpy(dtype=np.float64)
    n_expr = df["n_expr"].to_numpy(dtype=np.float64)
    pct = np.where(n_total > 0, n_expr / n_total, 0.0)
    mean = np.where(n_total > 0, df["sum_norm"].to_numpy(dtype=np.float64) / n_total, 0.0)

    deciles = [
        finalize_deciles(np.asarray(g, dtype=np.float64), grid_probs, float(p)).tolist()
        for g, p in zip(df["qgrid"], pct)
    ]

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
    table = table.append_column(
        "deciles", pa.array(deciles, type=pa.list_(pa.float64()))
    )
    return table


# --------------------------------------------------------------------------- #
# Driver                                                                      #
# --------------------------------------------------------------------------- #
def select_shards(data_dir: Path, n_shards: int | None, full: bool) -> list[str]:
    expr_dir = data_dir / "expression_data"
    files = sorted(str(p) for p in expr_dir.glob("train-*.parquet"))
    if not files:
        sys.exit(f"No expression shards found in {expr_dir}")
    if full or n_shards is None or n_shards >= len(files):
        return files
    return files[:n_shards]


def token_buckets(con: duckdb.DuckDBPyConnection, n_buckets: int):
    """Yield (lo, hi) half-open token ranges covering all gene tokens."""
    lo_tok, hi_tok = con.sql("SELECT min(token_id), max(token_id) FROM gene_map").fetchone()
    if n_buckets <= 1:
        yield (None, None)
        return
    edges = np.linspace(lo_tok, hi_tok + 1, n_buckets + 1).astype(int)
    for i in range(n_buckets):
        yield (int(edges[i]), int(edges[i + 1]))


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
    ap.add_argument("--shards", type=int, default=4,
                    help="number of expression shards to process (default: 4, a "
                         "Day-1 sample). Ignored when --full is set.")
    ap.add_argument("--full", action="store_true",
                    help=f"process all {N_SHARDS_TOTAL} shards (~95M cells; slow, heavy)")
    ap.add_argument("--grid-points", type=int, default=21,
                    help="quantile-grid resolution for the expressing distribution "
                         "(default: 21). Higher = finer violin silhouette.")
    ap.add_argument("--gene-batches", type=int, default=1,
                    help="split the aggregation into N token-range passes to bound "
                         "memory in --full mode (default: 1). Each pass rescans shards.")
    ap.add_argument("--threads", type=int, default=None, help="DuckDB threads")
    ap.add_argument("--memory-limit", default=None,
                    help="DuckDB memory limit, e.g. '16GB' (enables disk spill)")
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

    shards = select_shards(data_dir, args.shards, full)
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
        f"({mode}); gene_batches={args.gene_batches}")

    if out_dir.exists():
        shutil.rmtree(out_dir)
    agg_dir.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    if args.threads:
        con.execute(f"PRAGMA threads={args.threads}")
    if args.memory_limit:
        con.execute(f"PRAGMA memory_limit='{args.memory_limit}'")
    con.execute("PRAGMA enable_progress_bar")

    t0 = time.time()
    build_lookups(con, data_dir)

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

    write_gene_index(con, out_dir, seen_symbols)
    write_manifest(out_dir, data_dir, shards, args, full, total_rows,
                   len(seen_symbols), time.time() - t0)
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
                   total_rows: int, n_genes: int, elapsed: float) -> None:
    manifest = {
        "data_dir": str(data_dir),
        "sample": bool(args.sample),
        "n_shards_processed": len(shards),
        "n_shards_total": N_SHARDS_TOTAL,
        "full": bool(full),
        "normalization": "log1p(counts_per_10k)",
        "perturbation": "drug@concentration+unit",
        "special_token_min": SPECIAL_TOKEN_MIN,
        "grid_points": args.grid_points,
        "n_deciles": N_DECILES,
        "n_group_rows": total_rows,
        "n_genes": n_genes,
        "elapsed_seconds": round(elapsed, 1),
    }
    (out_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    log(f"manifest -> {out_dir / 'run_manifest.json'}")


if __name__ == "__main__":
    main()
