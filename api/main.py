"""Gene Explorer Query API (Member B) -- FastAPI service.

Serves the two endpoints in the project contract (project.md section 4) from
the precompute store built by `pipeline/build_aggregates.py`:

    GET /api/search?q=<prefix>   autocomplete -> {"matches": [{symbol, name}, ...]}
    GET /api/gene/{symbol}       gene page    -> {gene, heatmap{...}, violin[...]}

Run (from repo root):
    uvicorn api.main:app --reload --port 8000
    # or point at a specific store:
    STORE_DIR=/path/to/pipeline/out uvicorn api.main:app --port 8000

Interactive docs at /docs once running.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from .store import GeneNotFound, Store

app = FastAPI(
    title="Tahoe-100M Gene Explorer API",
    version="0.1.0",
    description="Gene search + per-gene expression aggregates for the Gene Explorer.",
)

# Frontend (Member C) calls this from the browser during dev; allow any origin.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET"],
    allow_headers=["*"],
)

# Single store instance, loaded once at import/startup.
store = Store()


@lru_cache(maxsize=1024)
def _gene_cached(symbol: str) -> dict:
    """Per-gene payloads are immutable for a given store build, so cache them."""
    return store.gene(symbol)


@app.get("/api/health")
def health() -> dict:
    """Liveness + which stores are being served."""
    return {
        "status": "ok",
        "stores": [str(r) for r in store.roots],
        "n_genes_indexed": len(store._index),
        "manifest": store.manifest,
    }


@app.get("/api/search")
def search(
    q: str = Query(..., min_length=1, description="gene symbol prefix, e.g. 'TAG'"),
    limit: int = Query(20, ge=1, le=100),
) -> dict:
    """Autocomplete over gene symbols. Prefix matches rank first."""
    return {"matches": store.search(q, limit=limit)}


@app.get("/api/gene/{symbol}")
def gene(symbol: str) -> dict:
    """Full gene payload: heatmap matrix + per-group violin data (section 4)."""
    try:
        return _gene_cached(symbol)
    except GeneNotFound:
        raise HTTPException(status_code=404, detail=f"gene '{symbol}' not found")


# Serve the built React frontend. Must come last so /api/* routes take precedence.
_dist = Path(__file__).parent.parent / "web" / "dist"
if _dist.exists():
    app.mount("/", StaticFiles(directory=_dist, html=True), name="frontend")
