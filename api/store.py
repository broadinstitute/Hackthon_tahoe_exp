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

_PIPELINE = Path(__file__).resolve().parent.parent / "pipeline"
DEFAULT_STORE = _PIPELINE / "out"
DEMO_STORE   = _PIPELINE / "out_demo"


class GeneNotFound(Exception):
    """Raised when a gene symbol has no partition in any store."""


class Store:
    """Multi-root store: roots are checked in order; first hit wins per gene.

    Typical setup: [out_demo, out]  — the demo store provides richer data for
    4 curated genes; the full out store covers the remaining ~35k genes.
    """

    def __init__(self, roots: list[str | os.PathLike] | None = None) -> None:
        if roots is None:
            # Env override still works for the primary store; demo store is
            # auto-discovered next to it when present.
            primary = Path(os.environ.get("STORE_DIR", DEFAULT_STORE)).resolve()
            candidate_demo = primary.parent / "out_demo"
            roots = ([candidate_demo, primary] if candidate_demo.is_dir()
                     else [primary])

        self.roots = [Path(r).resolve() for r in roots]

        # Validate at least one root has a gene index.
        valid = [r for r in self.roots if (r / "gene_index.parquet").exists()]
        if not valid:
            raise FileNotFoundError(
                f"No gene_index.parquet found in any of {self.roots}. "
                "Run `python pipeline/build_aggregates.py --sample`."
            )

        # Expose the primary (last-fallback) root for the health endpoint.
        self.root = self.roots[-1]

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
        """Merge indexes from all roots; first root wins for duplicates.

        Partition directory names are Hive-encoded (`gene_symbol=<url-quoted>`);
        we url-decode them so symbols with special characters resolve correctly.
        """
        seen_symbols: set[str] = set()

        for root in self.roots:
            index_path = root / "gene_index.parquet"
            if not index_path.exists():
                continue
            rows = self.con.execute(
                "SELECT symbol, name FROM read_parquet(?) ORDER BY symbol",
                [str(index_path)],
            ).fetchall()
            for sym, name in rows:
                if sym not in seen_symbols:
                    self._index.append({"symbol": sym, "name": name})
                    seen_symbols.add(sym)

            agg_dir = root / "aggregates"
            if agg_dir.exists():
                prefix = "gene_symbol="
                for entry in agg_dir.iterdir():
                    if entry.is_dir() and entry.name.startswith(prefix):
                        symbol = urllib.parse.unquote(entry.name[len(prefix):])
                        # First root to provide a gene's partition wins.
                        if symbol not in self._symbol_to_dir:
                            self._symbol_to_dir[symbol] = entry

        self._index.sort(key=lambda r: r["symbol"])

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
                "organ": r["organ"],
                "perturbation": r["perturbation"],
                "n": r["n_cells"],
                "deciles": r["deciles"],
                "pct_expressing": r["pct_expressing"],
                "mean": r["mean"],
            }
            for r in records
        ]

        ensembl_id = records[0]["ensembl_id"] if records else None
        # organ lookup keyed by cell_line name, for the frontend to annotate axes
        cl_organ: dict[str, str | None] = {}
        for r in records:
            if r["cell_line"] not in cl_organ:
                cl_organ[r["cell_line"]] = r["organ"]
        return {
            "gene": symbol,
            "ensembl_id": ensembl_id,
            "n_groups": len(records),
            "heatmap": {
                "cell_lines": cell_lines,
                "perturbations": perturbations,
                "mean": mean_matrix,
                "cell_line_organs": cl_organ,
            },
            "violin": violin,
        }
