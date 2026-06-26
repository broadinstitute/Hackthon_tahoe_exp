"""Data access for the Gene Explorer API (Member B).

Reads the precompute store produced by `pipeline/build_aggregates.py`:

    <store>/gene_index.parquet                       symbol -> name lookup (search)
    <store>/aggregates/gene_symbol=<SYM>/*.parquet   one partition per gene
    <store>/run_manifest.json                        provenance (optional)

Each aggregate row is one (cell_line x perturbation) group for that gene, with
the columns the API contract (project.md section 4) needs: cell_line,
perturbation, n_cells, pct_expressing, mean, deciles (+ extras we pass through).

The store is read-only and immutable for a given build, so we:
  * load the gene index fully into memory once at startup (search is then a
    pure in-memory scan -- no disk per keystroke), and
  * read a single gene's partition on demand via DuckDB, which is tiny
    (tens to low-thousands of rows). Per-gene results are cached upstream.
"""
from __future__ import annotations

import json
import os
import urllib.parse
from pathlib import Path
from typing import Any

import duckdb

DEFAULT_STORE = Path(__file__).resolve().parent.parent / "pipeline" / "out"


class GeneNotFound(Exception):
    """Raised when a gene symbol has no partition in the store."""


class Store:
    def __init__(self, root: str | os.PathLike | None = None) -> None:
        self.root = Path(root or os.environ.get("STORE_DIR", DEFAULT_STORE)).resolve()
        self.agg_dir = self.root / "aggregates"
        self.index_path = self.root / "gene_index.parquet"
        if not self.index_path.exists():
            raise FileNotFoundError(
                f"gene_index.parquet not found under {self.root}. "
                "Point STORE_DIR at a built pipeline output, or run "
                "`python pipeline/build_aggregates.py --sample`."
            )
        # One connection, reused read-only across requests. DuckDB connections
        # are safe to share for concurrent reads.
        self.con = duckdb.connect(database=":memory:")
        self.manifest = self._load_manifest()
        self._symbol_to_dir: dict[str, Path] = {}
        self._index: list[dict[str, str]] = []
        self._load_index()

    # -- startup loaders ---------------------------------------------------- #
    def _load_manifest(self) -> dict[str, Any]:
        mpath = self.root / "run_manifest.json"
        if mpath.exists():
            return json.loads(mpath.read_text())
        return {}

    def _load_index(self) -> None:
        """Load the search index and map every gene symbol to its partition dir.

        Partition directory names are Hive-encoded (`gene_symbol=<url-quoted>`);
        we url-decode them so symbols with special characters resolve to the
        right directory without guessing the encoding.
        """
        rows = self.con.execute(
            "SELECT symbol, name FROM read_parquet(?) ORDER BY symbol",
            [str(self.index_path)],
        ).fetchall()
        self._index = [{"symbol": s, "name": n} for s, n in rows]

        if self.agg_dir.exists():
            prefix = "gene_symbol="
            for entry in self.agg_dir.iterdir():
                if entry.is_dir() and entry.name.startswith(prefix):
                    symbol = urllib.parse.unquote(entry.name[len(prefix):])
                    self._symbol_to_dir[symbol] = entry

    # -- queries ------------------------------------------------------------ #
    def search(self, q: str, limit: int = 20) -> list[dict[str, str]]:
        """Autocomplete: case-insensitive, prefix matches first then substring."""
        q = (q or "").strip().lower()
        if not q:
            return []
        prefix: list[dict[str, str]] = []
        contains: list[dict[str, str]] = []
        for row in self._index:
            sym = row["symbol"].lower()
            if sym.startswith(q):
                prefix.append(row)
            elif q in sym or q in row["name"].lower():
                contains.append(row)
            if len(prefix) >= limit:
                break
        return (prefix + contains)[:limit]

    def gene(self, symbol: str) -> dict[str, Any]:
        """Return the section-4 contract payload for one gene symbol.

        Raises GeneNotFound if the symbol has no partition.
        """
        part_dir = self._symbol_to_dir.get(symbol)
        if part_dir is None:
            raise GeneNotFound(symbol)

        glob = str(part_dir / "*.parquet")
        rows = self.con.execute(
            """
            SELECT cell_line, perturbation, n_cells, pct_expressing, mean,
                   deciles, drug, concentration, conc_unit, organ,
                   n_expressing, ensembl_id
            FROM read_parquet(?)
            ORDER BY cell_line, perturbation
            """,
            [glob],
        ).fetchall()

        cols = [d[0] for d in self.con.description]
        records = [dict(zip(cols, r)) for r in rows]

        cell_lines = sorted({r["cell_line"] for r in records})
        perturbations = sorted({r["perturbation"] for r in records})
        cl_idx = {c: i for i, c in enumerate(cell_lines)}
        pt_idx = {p: i for i, p in enumerate(perturbations)}

        # Heatmap: mean expression per cell_line (row) x perturbation (col);
        # None where a group was not observed.
        mean_matrix: list[list[float | None]] = [
            [None] * len(perturbations) for _ in cell_lines
        ]
        for r in records:
            mean_matrix[cl_idx[r["cell_line"]]][pt_idx[r["perturbation"]]] = r["mean"]

        violin = [
            {
                "cell_line": r["cell_line"],
                "perturbation": r["perturbation"],
                "n": r["n_cells"],
                "deciles": r["deciles"],
                "pct_expressing": r["pct_expressing"],
                "mean": r["mean"],
            }
            for r in records
        ]

        ensembl_id = records[0]["ensembl_id"] if records else None
        return {
            "gene": symbol,
            "ensembl_id": ensembl_id,
            "n_groups": len(records),
            "heatmap": {
                "cell_lines": cell_lines,
                "perturbations": perturbations,
                "mean": mean_matrix,
            },
            "violin": violin,
        }
