import { useState } from "react";
import SearchBar from "./components/SearchBar";
import Heatmap from "./components/Heatmap";
import ViolinPlot from "./components/ViolinPlot";
import { fetchGene, fetchDemoGene, DEMO_GENES, type DemoGene, type GeneResponse } from "./api";
import "./App.css";

type State = "idle" | "loading" | "ok" | "error";

const GENE_DESC: Record<DemoGene, string> = {
  KRAS: "oncogene · RAS signalling",
  TP53: "tumour suppressor",
  EGFR: "receptor tyrosine kinase",
  BRCA1: "DNA repair",
};

export default function App() {
  const [state, setState] = useState<State>("idle");
  const [data, setData] = useState<GeneResponse | null>(null);
  const [error, setError] = useState("");
  const [loadingGene, setLoadingGene] = useState<string | null>(null);

  async function handleSelect(symbol: string) {
    setState("loading");
    setData(null);
    setError("");
    setLoadingGene(symbol);
    try {
      const result = await fetchGene(symbol);
      setData(result);
      setState("ok");
    } catch (e: unknown) {
      setError(e instanceof Error ? e.message : "Unknown error");
      setState("error");
    } finally {
      setLoadingGene(null);
    }
  }

  async function handleDemo(gene: DemoGene) {
    setState("loading");
    setData(null);
    setError("");
    setLoadingGene(gene);
    try {
      const result = await fetchDemoGene(gene);
      setData(result);
      setState("ok");
    } catch (e: unknown) {
      setError(e instanceof Error ? e.message : "Unknown error");
      setState("error");
    } finally {
      setLoadingGene(null);
    }
  }

  return (
    <div className="app">
      <header className="app-header">
        <div className="header-inner">
          <div className="logo">
            <span className="logo-mark">T</span>
            <span className="logo-text">Tahoe-100M Gene Explorer</span>
          </div>
          <SearchBar onSelect={handleSelect} />
        </div>
      </header>

      <main className="app-main">
        {state === "idle" && (
          <div className="splash">
            <h1>Search a gene to explore its expression</h1>
            <p>
              Single-cell expression aggregates across 50+ cancer cell lines and 1,100+
              perturbations from the Tahoe-100M dataset.
            </p>

            <div className="demo-section">
              <p className="demo-label">Featured genes — click to explore</p>
              <div className="demo-grid">
                {DEMO_GENES.map(gene => (
                  <button
                    key={gene}
                    className="demo-card"
                    onClick={() => handleDemo(gene)}
                  >
                    <span className="demo-card-gene">{gene}</span>
                    <span className="demo-card-desc">{GENE_DESC[gene]}</span>
                  </button>
                ))}
              </div>
            </div>
          </div>
        )}

        {state === "loading" && (
          <div className="status-msg">
            <div className="spinner" />
            Loading {loadingGene ?? ""}…
          </div>
        )}

        {state === "error" && (
          <div className="status-msg error">{error}</div>
        )}

        {state === "ok" && data && (
          <div className="gene-page">
            <div className="gene-title">
              <h2>{data.gene}</h2>
              {data.ensembl_id && (
                <span className="ensembl">{data.ensembl_id}</span>
              )}
              <span className="n-groups">{data.n_groups.toLocaleString()} groups</span>
              <button className="back-btn" onClick={() => { setState("idle"); setData(null); }}>
                ← Back
              </button>
            </div>

            <div className="plot-grid">
              <section className="plot-card heatmap-card">
                <Heatmap data={data.heatmap} gene={data.gene} />
              </section>
              <section className="plot-card violin-card">
                <ViolinPlot violin={data.violin} gene={data.gene} />
              </section>
            </div>
          </div>
        )}
      </main>
    </div>
  );
}
