import { useEffect, useRef, useState, useMemo } from "react";
import Plotly from "plotly.js-dist-min";
import { type ViolinEntry } from "../api";
import "./ViolinPlot.css";

interface Props {
  violin: ViolinEntry[];
  gene: string;
}

function synthesizePoints(deciles: number[], n = 120): number[] {
  const probs = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
  const dMin = deciles[0];
  const dMax = deciles[deciles.length - 1];
  const range = dMax - dMin;
  // 0.8% jitter breaks point-mass singularities (e.g. 90% cells at 0),
  // but clamped to [dMin, dMax] so no negative expression values appear.
  const jitter = range > 0 ? range * 0.008 : 0;

  let seed = 42;
  const rand = () => {
    seed = (seed * 1664525 + 1013904223) & 0xffffffff;
    return (seed >>> 0) / 0xffffffff;
  };

  return Array.from({ length: n }, (_, i) => {
    const p = n === 1 ? 0.5 : i / (n - 1);
    let val = dMax;
    for (let j = 0; j < probs.length - 1; j++) {
      if (p <= probs[j + 1]) {
        const span = probs[j + 1] - probs[j];
        const t = span === 0 ? 0 : (p - probs[j]) / span;
        val = deciles[j] + t * (deciles[j + 1] - deciles[j]);
        break;
      }
    }
    // Clamp to data range — no negative expression values from jitter
    const noisy = val + (rand() - 0.5) * jitter;
    return Math.max(dMin, Math.min(dMax, noisy));
  });
}

// Compute a reasonable KDE bandwidth from an array of decile sets.
// Uses the IQR of all non-zero decile midpoints across the dataset.
function estimateBandwidth(entries: ViolinEntry[]): number {
  const mids = entries.map(v => v.deciles[5]); // median decile
  mids.sort((a, b) => a - b);
  const q25 = mids[Math.floor(mids.length * 0.25)] ?? 0;
  const q75 = mids[Math.floor(mids.length * 0.75)] ?? 0;
  const iqr = q75 - q25;
  return Math.max(0.08, iqr * 0.5);
}

function buildGroupedTrace(
  categoryPoints: Array<{ cat: string; pts: number[] }>,
  traceName: string,
  isControl: boolean,
  bandwidth: number,
): Plotly.Data {
  const x: string[] = [];
  const y: number[] = [];
  for (const { cat, pts } of categoryPoints) {
    for (const pt of pts) { x.push(cat); y.push(pt); }
  }
  return {
    type: "violin" as const,
    orientation: "v",
    x,
    y,
    name: traceName,
    points: false,
    marker: { opacity: 0, size: 0 },
    box: { visible: false },
    meanline: { visible: true },
    bandwidth,
    spanmode: "hard",
    opacity: isControl ? 0.45 : 0.75,
    ...(isControl ? {
      fillcolor: "rgba(160,160,160,0.3)",
      line: { color: "rgba(100,100,100,0.55)" },
      meanline: { visible: true, color: "rgba(80,80,80,0.7)" },
    } : {}),
    hovertemplate: `<b>%{x}</b><br>${traceName.length > 30 ? traceName.slice(0, 28) + "…" : traceName}<extra></extra>`,
  };
}

type ViewMode = "by-cell-line" | "by-perturbation";

const MAX_CELL_LINES = 12;
const N_PRESETS = 5;
// Top-N perturbations to show per cell line in by-perturbation view
const MAX_PERTS_SHOWN = 40;
const CTRL_RE = /^DMSO/i;

function pickDefaultPert(perts: string[]): string {
  return perts.find(p => CTRL_RE.test(p)) ?? perts[0] ?? "";
}

function topPerturbations(violin: ViolinEntry[], n: number): string[] {
  const totals: Record<string, number> = {};
  for (const v of violin) totals[v.perturbation] = (totals[v.perturbation] ?? 0) + v.n;
  const ctrl = Object.keys(totals).find(p => CTRL_RE.test(p));
  const sorted = Object.entries(totals)
    .filter(([p]) => p !== ctrl)
    .sort((a, b) => b[1] - a[1])
    .map(([p]) => p);
  return ctrl ? [ctrl, ...sorted.slice(0, n - 1)] : sorted.slice(0, n);
}

export default function ViolinPlot({ violin, gene }: Props) {
  const ref = useRef<HTMLDivElement>(null);
  const [view, setView] = useState<ViewMode>("by-cell-line");
  const [groupByTissue, setGroupByTissue] = useState(false);
  const [selectedCL, setSelectedCL] = useState<string>("");
  const [selectedPert, setSelectedPert] = useState<string>("");
  const [expandedCL, setExpandedCL] = useState(false);
  const [searchText, setSearchText] = useState<string>("");
  const [dropdownOpen, setDropdownOpen] = useState(false);

  const cellLines = useMemo(() => [...new Set(violin.map(v => v.cell_line))].sort(), [violin]);
  const perturbations = useMemo(() => [...new Set(violin.map(v => v.perturbation))].sort(), [violin]);
  const presets = useMemo(() => topPerturbations(violin, N_PRESETS), [violin]);
  const ctrlPert = useMemo(() => perturbations.find(p => CTRL_RE.test(p)) ?? "", [perturbations]);

  // Reset selections when gene/dataset changes
  useEffect(() => {
    setSelectedCL("");
    setSelectedPert("");
  }, [violin]);

  useEffect(() => {
    if (cellLines.length && !selectedCL) setSelectedCL(cellLines[0]);
  }, [cellLines, selectedCL]);

  useEffect(() => {
    if (perturbations.length && !selectedPert) setSelectedPert(pickDefaultPert(perturbations));
  }, [perturbations, selectedPert]);

  const dropdownOptions = useMemo(() =>
    !searchText.trim() ? perturbations
      : perturbations.filter(p => p.toLowerCase().includes(searchText.toLowerCase())),
    [perturbations, searchText]
  );

  const { traces } = useMemo(() => {
    // ── By perturbation ──
    if (view === "by-perturbation") {
      const cellEntries = violin.filter(v => v.cell_line === selectedCL);
      // Limit to top-N perturbations by n_cells to keep chart readable
      const byN = [...cellEntries].sort((a, b) => b.n - a.n);
      const ctrl = byN.find(v => CTRL_RE.test(v.perturbation));
      const topEntries = byN.filter(v => !CTRL_RE.test(v.perturbation)).slice(0, MAX_PERTS_SHOWN);
      const shown = ctrl ? [ctrl, ...topEntries] : topEntries;

      const bw = estimateBandwidth(shown);
      const ctrlPts = ctrl ? synthesizePoints(ctrl.deciles) : [];

      const ctrlTrace = ctrlPts.length ? buildGroupedTrace(
        shown.map(v => ({ cat: v.perturbation, pts: ctrlPts })),
        ctrl!.perturbation, true, bw,
      ) : null;

      const pertTrace = buildGroupedTrace(
        shown.map(v => ({ cat: v.perturbation, pts: synthesizePoints(v.deciles) })),
        selectedCL, false, bw,
      );

      return { traces: [ctrlTrace, pertTrace].filter(Boolean) as Plotly.Data[], bw };
    }

    // ── By cell line (grouped by tissue) ──
    if (groupByTissue) {
      const selectedSubset = violin.filter(v => v.perturbation === selectedPert && v.organ);
      const ctrlSubset    = violin.filter(v => v.perturbation === ctrlPert && v.organ);
      const bw = estimateBandwidth([...selectedSubset, ...ctrlSubset]);

      const poolByOrgan = (entries: ViolinEntry[]) => {
        const map: Record<string, number[]> = {};
        // 50 pts per cell line keeps data size reasonable while still giving good KDE
        for (const v of entries) (map[v.organ!] ??= []).push(...synthesizePoints(v.deciles, 50));
        return Object.entries(map)
          .sort(([a], [b]) => a.localeCompare(b))
          .map(([cat, pts]) => ({ cat, pts }));
      };

      const out: Plotly.Data[] = [];
      if (ctrlPert && selectedPert !== ctrlPert && ctrlSubset.length) {
        out.push(buildGroupedTrace(poolByOrgan(ctrlSubset), "Control (DMSO)", true, bw));
      }
      if (selectedSubset.length) {
        out.push(buildGroupedTrace(poolByOrgan(selectedSubset), selectedPert, false, bw));
      }
      return { traces: out, bw };
    }

    // ── By cell line (individual) ──
    const selectedSubset = violin.filter(v => v.perturbation === selectedPert);
    const ctrlMap: Record<string, ViolinEntry> = {};
    for (const v of violin.filter(v => v.perturbation === ctrlPert)) ctrlMap[v.cell_line] = v;
    const bw = estimateBandwidth(selectedSubset);

    const label = (v: ViolinEntry) => v.organ ? `${v.cell_line} (${v.organ})` : v.cell_line;
    const out: Plotly.Data[] = [];

    if (ctrlPert && selectedPert !== ctrlPert) {
      const ctrlCats = selectedSubset
        .filter(v => ctrlMap[v.cell_line])
        .map(v => ({ cat: label(v), pts: synthesizePoints(ctrlMap[v.cell_line].deciles) }));
      if (ctrlCats.length) out.push(buildGroupedTrace(ctrlCats, "Control (DMSO)", true, bw));
    }
    out.push(buildGroupedTrace(
      selectedSubset.map(v => ({ cat: label(v), pts: synthesizePoints(v.deciles) })),
      selectedPert, false, bw,
    ));
    return { traces: out, bw };
  }, [view, violin, selectedCL, selectedPert, groupByTissue, ctrlPert]);

  const xLabel = view === "by-perturbation" ? "Perturbation"
    : groupByTissue ? "Tissue / Organ" : "Cell Line";
  const subtitle = view === "by-perturbation" ? selectedCL : selectedPert;

  useEffect(() => {
    if (!ref.current || !traces.length) return;
    Plotly.react(ref.current, traces, {
      title: { text: `${gene} — Expression by ${xLabel}<br><sup>${subtitle}</sup>`, font: { size: 14 } },
      violinmode: "group",
      showlegend: true,
      legend: { orientation: "h", y: -0.25 },
      yaxis: {
        title: { text: "Expression (log1p CPT10K)" },
        zeroline: true, zerolinecolor: "#cfd8dc",
        rangemode: "tozero",
        automargin: true,
      },
      xaxis: {
        type: "category",          // force categorical — prevents Plotly numeric-index fallback
        tickangle: -45,
        tickfont: { size: 9 },
        title: { text: xLabel },
        automargin: true,
      },
      height: 460,
      margin: { l: 60, r: 20, t: 70, b: 160 },
      plot_bgcolor: "#fafcff",
      paper_bgcolor: "#fff",
    } as Partial<Plotly.Layout>, { responsive: true, displayModeBar: false });
  }, [traces, gene, xLabel, subtitle]);

  function selectPert(p: string) {
    setSelectedPert(p);
    setSearchText("");
    setDropdownOpen(false);
  }

  const pertLimitNote = view === "by-perturbation"
    ? violin.filter(v => v.cell_line === selectedCL).length > MAX_PERTS_SHOWN + 1
      ? `Top ${MAX_PERTS_SHOWN} perturbations by cell count`
      : null
    : null;

  return (
    <div className="violin-wrap">
      <div className="controls-row">
        <div className="view-toggle">
          {(["by-cell-line", "by-perturbation"] as ViewMode[]).map(v => (
            <button key={v} className={view === v ? "vt-btn active" : "vt-btn"} onClick={() => setView(v)}>
              {v === "by-cell-line" ? "By cell line" : "By perturbation"}
            </button>
          ))}
        </div>
        {view === "by-cell-line" && (
          <label className="tissue-toggle">
            <input type="checkbox" checked={groupByTissue} onChange={e => setGroupByTissue(e.target.checked)} />
            Group by tissue
          </label>
        )}
        {pertLimitNote && <span className="pert-limit-note">{pertLimitNote}</span>}
      </div>

      {view === "by-cell-line" && (
        <div className="condition-selector">
          <span className="cl-label">Condition:</span>
          <div className="preset-chips">
            {presets.map(p => (
              <button key={p} className={p === selectedPert ? "preset-chip active" : "preset-chip"}
                onClick={() => selectPert(p)} title={p}>
                {CTRL_RE.test(p) ? "Control (DMSO)" : p.length > 18 ? p.slice(0, 16) + "…" : p}
              </button>
            ))}
          </div>
          <div className="pert-dropdown">
            <input className="pert-search" type="text" placeholder="Search all conditions…"
              value={dropdownOpen ? searchText : (presets.includes(selectedPert) ? "" : selectedPert)}
              onChange={e => { setSearchText(e.target.value); setDropdownOpen(true); }}
              onFocus={() => { setSearchText(""); setDropdownOpen(true); }}
              onBlur={() => setTimeout(() => setDropdownOpen(false), 160)}
            />
            {dropdownOpen && (
              <ul className="pert-list">
                {dropdownOptions.slice(0, 40).map(p => (
                  <li key={p} className={p === selectedPert ? "active" : ""} onMouseDown={() => selectPert(p)}>
                    {CTRL_RE.test(p) && <span className="pert-badge">Control</span>}{p}
                  </li>
                ))}
                {dropdownOptions.length > 40 && <li className="pert-overflow">{dropdownOptions.length - 40} more — keep typing</li>}
              </ul>
            )}
          </div>
        </div>
      )}

      {view === "by-perturbation" && (
        <div className="cl-tabs">
          <span className="cl-label">Cell line:</span>
          {cellLines.slice(0, expandedCL ? cellLines.length : MAX_CELL_LINES).map(cl => (
            <button key={cl} className={cl === selectedCL ? "cl-tab active" : "cl-tab"} onClick={() => setSelectedCL(cl)}>{cl}</button>
          ))}
          {cellLines.length > MAX_CELL_LINES && (
            <button className="cl-more" onClick={() => setExpandedCL(e => !e)}>
              {expandedCL ? "Show less" : `+${cellLines.length - MAX_CELL_LINES} more`}
            </button>
          )}
        </div>
      )}

      <div ref={ref} style={{ width: "100%" }} />
    </div>
  );
}
