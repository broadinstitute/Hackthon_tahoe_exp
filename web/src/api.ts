const BASE = import.meta.env.VITE_API_URL ?? "http://localhost:8000";

export interface SearchMatch {
  symbol: string;
  name: string;
}

export interface HeatmapData {
  cell_lines: string[];
  perturbations: string[];
  mean: (number | null)[][];
  cell_line_organs: Record<string, string | null>;
}

export interface ViolinEntry {
  cell_line: string;
  organ: string | null;
  perturbation: string;
  n: number;
  deciles: number[];
  pct_expressing: number;
  mean: number;
}

export interface GeneResponse {
  gene: string;
  ensembl_id: string | null;
  n_groups: number;
  heatmap: HeatmapData;
  violin: ViolinEntry[];
}

export async function searchGenes(q: string): Promise<SearchMatch[]> {
  if (!q.trim()) return [];
  try {
    const r = await fetch(`${BASE}/api/search?q=${encodeURIComponent(q)}&limit=10`);
    if (!r.ok) return [];
    const data = await r.json();
    return data.matches as SearchMatch[];
  } catch {
    return [];
  }
}

export async function fetchGene(symbol: string): Promise<GeneResponse> {
  const r = await fetch(`${BASE}/api/gene/${encodeURIComponent(symbol)}`);
  if (r.status === 404) throw new Error(`Gene "${symbol}" not found`);
  if (!r.ok) throw new Error(`Server error: ${r.status}`);
  return r.json();
}

// Load a bundled static demo gene without needing the API running.
export const DEMO_GENES = ["KRAS", "TP53", "EGFR", "BRCA1"] as const;
export type DemoGene = (typeof DEMO_GENES)[number];

export async function fetchDemoGene(gene: DemoGene): Promise<GeneResponse> {
  const r = await fetch(`/${gene}.json`);
  if (!r.ok) throw new Error(`Demo data for ${gene} not found`);
  return r.json();
}
