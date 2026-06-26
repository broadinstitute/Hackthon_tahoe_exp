"""Contract tests for the Gene Explorer API.

These run against the real sample store under pipeline/out (built by
`python pipeline/build_aggregates.py --sample`). They are skipped with a clear
message if the store is absent, so they never fail spuriously on a fresh clone.

    cd api && pytest            # or: pytest api/tests
"""
from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from api.main import app, store

client = TestClient(app)


def _a_gene_with_data() -> str:
    """A symbol that actually has a partition in the current store."""
    assert store._symbol_to_dir, "store has no gene partitions"
    return next(iter(sorted(store._symbol_to_dir)))


def test_health():
    r = client.get("/api/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"
    assert r.json()["n_genes_indexed"] > 0


def test_search_prefix_ranks_first():
    sym = _a_gene_with_data()
    r = client.get("/api/search", params={"q": sym[:3]})
    assert r.status_code == 200
    matches = r.json()["matches"]
    assert matches, "expected at least one match"
    # every returned row has the contract shape {symbol, name}
    assert all({"symbol", "name"} <= m.keys() for m in matches)
    # first result is a genuine prefix match
    assert matches[0]["symbol"].lower().startswith(sym[:3].lower())


def test_search_empty_query_rejected():
    # min_length=1 -> 422 from FastAPI validation
    assert client.get("/api/search", params={"q": ""}).status_code == 422


def test_gene_contract_shape():
    sym = _a_gene_with_data()
    r = client.get(f"/api/gene/{sym}")
    assert r.status_code == 200
    body = r.json()
    assert body["gene"] == sym
    hm = body["heatmap"]
    assert {"cell_lines", "perturbations", "mean"} <= hm.keys()
    # mean matrix dimensions match the axes
    assert len(hm["mean"]) == len(hm["cell_lines"])
    if hm["cell_lines"]:
        assert all(len(row) == len(hm["perturbations"]) for row in hm["mean"])
    # violin entries carry the section-4 fields
    assert body["violin"], "expected violin groups"
    v = body["violin"][0]
    assert {"cell_line", "perturbation", "n", "deciles", "pct_expressing"} <= v.keys()
    assert len(v["deciles"]) == 11


def test_gene_not_found():
    r = client.get("/api/gene/NOT_A_REAL_GENE_XYZ")
    assert r.status_code == 404
