import { useEffect, useRef } from "react";
import Plotly from "plotly.js-dist-min";
import { type HeatmapData } from "../api";

interface Props {
  data: HeatmapData;
  gene: string;
}

// Match only the DMSO vehicle control, not drugs that happen to contain "DMSO" in their name.
const CTRL_RE = /^DMSO/i;

function reorderPerts(perts: string[], mean: (number | null)[][]): {
  perts: string[];
  mean: (number | null)[][];
} {
  const ctrlIdx = perts.findIndex(p => CTRL_RE.test(p));
  if (ctrlIdx <= 0) return { perts, mean }; // already first or not found
  const newPerts = [perts[ctrlIdx], ...perts.filter((_, i) => i !== ctrlIdx)];
  const newMean = mean.map(row => [row[ctrlIdx], ...row.filter((_, i) => i !== ctrlIdx)]);
  return { perts: newPerts, mean: newMean };
}

export default function Heatmap({ data, gene }: Props) {
  const ref = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!ref.current) return;

    const { perts, mean } = reorderPerts(data.perturbations, data.mean);
    const nRows = data.cell_lines.length;

    // Percentile-based color range so variation is visible even when values cluster.
    const flat = mean.flat().filter((v): v is number => v !== null && !isNaN(v));
    flat.sort((a, b) => a - b);
    const p2  = flat[Math.floor(flat.length * 0.02)]  ?? 0;
    const p98 = flat[Math.floor(flat.length * 0.98)]  ?? flat[flat.length - 1] ?? 1;

    // y-axis labels: "CellLine (Organ)" when organ is available
    const yLabels = data.cell_lines.map(cl => {
      const organ = data.cell_line_organs?.[cl];
      return organ ? `${cl} (${organ})` : cl;
    });

    // "Control" annotation above the first column
    const hasCtrl = CTRL_RE.test(perts[0] ?? "");
    const annotations: object[] = hasCtrl ? [{
      xref: "x",
      yref: "paper",
      x: 0,
      y: 1.02,
      text: "▼ Control",
      showarrow: false,
      font: { size: 10, color: "#1565c0", weight: 700 },
      xanchor: "center",
      yanchor: "bottom",
    }] : [];

    // Professional sequential colorscale: pale cream → teal → deep navy.
    // Perceptually uniform and readable in print/projector.
    const colorscale: [number, string][] = [
      [0,    "#f7fbff"],
      [0.15, "#c6dbef"],
      [0.35, "#6baed6"],
      [0.60, "#2171b5"],
      [0.80, "#084594"],
      [1,    "#08306b"],
    ];

    const trace: Plotly.Data = {
      type: "heatmap",
      x: perts,
      y: yLabels,
      z: mean,
      colorscale,
      reversescale: false,
      zmin: p2,
      zmax: p98,
      colorbar: { title: { text: "Mean expr<br>(log1p CPT10K)", side: "right" }, thickness: 14 },
      hoverongaps: false,
    };

    const layout = {
      title: { text: `${gene} — Mean Expression Heatmap`, font: { size: 15 } },
      margin: { l: 140, r: 20, t: 60, b: 160 },
      xaxis: {
        tickangle: -45,
        tickfont: { size: 9 },
        title: { text: "Perturbation" },
        automargin: true,
      },
      yaxis: {
        tickfont: { size: 10 },
        title: { text: "Cell Line" },
        automargin: true,
      },
      annotations,
      // Height scales with rows; width is handled by the container (responsive).
      // More rows → taller; cap at 900 to stay scrollable rather than huge.
      height: Math.min(900, Math.max(320, nRows * 26 + 200)),
      plot_bgcolor: "#fafcff",
      paper_bgcolor: "#fff",
    };

    Plotly.react(ref.current, [trace], layout, { responsive: true, displayModeBar: false });
  }, [data, gene]);

  return <div ref={ref} style={{ width: "100%" }} />;
}
