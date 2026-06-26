#!/usr/bin/env python3
"""Build a tiny, self-contained Tahoe-100M sample dataset.

The real Tahoe-100M dataset is ~328 GB of per-cell expression Parquet (3388
shards) plus a 2.2 GB obs table and an 83 GB pseudobulk table. That is far too
large to copy to a laptop or a CI box. This script carves out a *tiny* slice
that mirrors the real directory layout exactly, so the precompute pipeline
(`build_aggregates.py`) — and anything else that reads the dataset — runs
against it unchanged, just with `--data-dir` pointed at the sample.

What it produces (under --out, default pipeline/sample_data/):

    expression_data/train-00000-of-00001.parquet   sampled cells (the big one)
    gene_metadata/gene_metadata.parquet             copied whole (small)
    cell_line_metadata/cell_line_metadata.parquet   copied whole (small)
    drug_metadata/drug_metadata.parquet             copied whole (small)
    sample_metadata/sample_metadata.parquet         copied whole (small)
    SAMPLE_README.md                                provenance + how to use

The metadata tables are lookups (gene token -> symbol, CVCL -> cell name,
sample -> dose); they are kept whole and tiny, so every reference in the
sampled cells still resolves. We deliberately skip obs_metadata/ and
pseudobulk_differential_expression/ — the pipeline does not read them.

Cells are drawn from several shards because each shard covers only a few drugs;
sampling across shards gives the cell-line x perturbation grid the Gene
Explorer needs. Sampling is seeded, so reruns are reproducible.
"""
from __future__ import annotations

import argparse
import shutil
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

DEFAULT_DATA_DIR = "/home/jovyan/organization/raw/public-datasets/tahoe_100m"

# Small lookup tables copied verbatim: (subdir, filename).
META_FILES = [
    ("gene_metadata", "gene_metadata.parquet"),
    ("cell_line_metadata", "cell_line_metadata.parquet"),
    ("drug_metadata", "drug_metadata.parquet"),
    ("sample_metadata", "sample_metadata.parquet"),
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def pick_shards(expr_dir: Path, n_shards: int, stride: int) -> list[Path]:
    """Spread picks across the shard range so we see diverse drugs."""
    files = sorted(expr_dir.glob("train-*.parquet"))
    if not files:
        raise SystemExit(f"No expression shards found in {expr_dir}")
    picked = files[::stride][:n_shards]
    if len(picked) < n_shards:
        picked = files[:n_shards]
    return picked


def sample_shard(path: Path, rows: int, rng: np.random.Generator) -> pa.Table:
    """Read one shard and return a random subset of `rows` cells."""
    table = pq.read_table(path)
    n = table.num_rows
    if rows >= n:
        return table
    idx = np.sort(rng.choice(n, size=rows, replace=False))
    return table.take(pa.array(idx))


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--data-dir", default=DEFAULT_DATA_DIR, type=Path,
                    help="real Tahoe-100M root (default: %(default)s)")
    ap.add_argument("--out", default=Path(__file__).parent / "sample_data", type=Path,
                    help="output directory for the sample (default: pipeline/sample_data)")
    ap.add_argument("--shards", type=int, default=24,
                    help="number of shards to draw cells from (default: 24)")
    ap.add_argument("--rows-per-shard", type=int, default=500,
                    help="cells sampled per shard (default: 500)")
    ap.add_argument("--stride", type=int, default=140,
                    help="shard stride for spreading picks across drugs "
                         "(default: 140 ~ 24 picks over 3376 shards)")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed (default: 0)")
    args = ap.parse_args()

    data_dir = args.data_dir.resolve()
    out_dir = args.out.resolve()
    expr_dir = data_dir / "expression_data"
    rng = np.random.default_rng(args.seed)

    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    # 1. Sample expression cells across shards.
    shards = pick_shards(expr_dir, args.shards, args.stride)
    log(f"sampling {args.rows_per_shard} cells from each of {len(shards)} shards")
    parts = []
    for i, sp in enumerate(shards):
        parts.append(sample_shard(sp, args.rows_per_shard, rng))
        if (i + 1) % 5 == 0 or i == len(shards) - 1:
            log(f"  {i + 1}/{len(shards)} shards read")
    cells = pa.concat_tables(parts)

    (out_dir / "expression_data").mkdir()
    expr_out = out_dir / "expression_data" / "train-00000-of-00001.parquet"
    pq.write_table(cells, expr_out, compression="zstd")
    expr_mb = expr_out.stat().st_size / 1e6

    # 2. Copy the small lookup tables verbatim.
    for subdir, fname in META_FILES:
        src = data_dir / subdir / fname
        if not src.exists():
            log(f"  WARN: missing {src}, skipping")
            continue
        (out_dir / subdir).mkdir(exist_ok=True)
        shutil.copy2(src, out_dir / subdir / fname)

    # 3. Summary stats for the README.
    cl = sorted(set(cells.column("cell_line_id").to_pylist()))
    drugs = sorted(set(cells.column("drug").to_pylist()))
    samples = sorted(set(cells.column("sample").to_pylist()))
    total_mb = sum(p.stat().st_size for p in out_dir.rglob("*")) / 1e6

    readme = f"""# Tahoe-100M — tiny sample

Generated by `pipeline/make_sample.py` from the full dataset at
`{data_dir}`.

This is a downsampled slice for development on machines without the full
~328 GB dataset. It mirrors the real directory layout, so the precompute
pipeline runs against it unchanged:

```bash
python pipeline/build_aggregates.py \\
    --data-dir pipeline/sample_data --full --out pipeline/out
```

(`--full` here just means "use all shards in this sample" — there is one.)

## Contents
- `expression_data/train-00000-of-00001.parquet` — {cells.num_rows:,} cells
- `gene_metadata/`, `cell_line_metadata/`, `drug_metadata/`,
  `sample_metadata/` — full lookup tables, copied verbatim (small), so every
  gene / cell-line / sample reference in the sampled cells resolves.

Skipped from the source (not read by the pipeline): `obs_metadata/` (2.2 GB),
`pseudobulk_differential_expression/` (83 GB).

## Sample coverage
- cells: {cells.num_rows:,}
- distinct cell lines: {len(cl)}
- distinct drugs: {len(drugs)}
- distinct samples (drug x dose): {len(samples)}
- shards sampled: {len(shards)} (seed={args.seed}, {args.rows_per_shard} cells each)
- total size on disk: {total_mb:.1f} MB

Because cells are spread thin across many (gene x cell_line x perturbation)
groups, per-group cell counts are small — this is a structural/dev sample, not
statistically representative of the full dataset.
"""
    (out_dir / "SAMPLE_README.md").write_text(readme)

    log(f"expression parquet: {expr_mb:.1f} MB, {cells.num_rows:,} cells")
    log(f"coverage: {len(cl)} cell lines, {len(drugs)} drugs, {len(samples)} samples")
    log(f"DONE: {total_mb:.1f} MB total -> {out_dir}")


if __name__ == "__main__":
    main()
