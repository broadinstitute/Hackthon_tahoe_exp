import { useState, useEffect, useRef } from "react";
import { searchGenes, type SearchMatch } from "../api";
import "./SearchBar.css";

interface Props {
  onSelect: (symbol: string) => void;
}

export default function SearchBar({ onSelect }: Props) {
  const [query, setQuery] = useState("");
  const [suggestions, setSuggestions] = useState<SearchMatch[]>([]);
  const [open, setOpen] = useState(false);
  const [highlighted, setHighlighted] = useState(-1);
  const debounce = useRef<ReturnType<typeof setTimeout> | null>(null);

  useEffect(() => {
    if (debounce.current) clearTimeout(debounce.current);
    if (!query.trim()) { setSuggestions([]); setOpen(false); return; }
    debounce.current = setTimeout(async () => {
      const matches = await searchGenes(query);
      setSuggestions(matches);
      setOpen(matches.length > 0);
      setHighlighted(-1);
    }, 180);
  }, [query]);

  function commit(symbol: string) {
    setQuery(symbol);
    setOpen(false);
    onSelect(symbol);
  }

  function onKeyDown(e: React.KeyboardEvent) {
    if (!open) return;
    if (e.key === "ArrowDown") { e.preventDefault(); setHighlighted(h => Math.min(h + 1, suggestions.length - 1)); }
    if (e.key === "ArrowUp")   { e.preventDefault(); setHighlighted(h => Math.max(h - 1, 0)); }
    if (e.key === "Enter") {
      e.preventDefault();
      const target = highlighted >= 0 ? suggestions[highlighted] : suggestions[0];
      if (target) commit(target.symbol);
    }
    if (e.key === "Escape") setOpen(false);
  }

  return (
    <div className="search-wrap">
      <div className="search-box">
        <span className="search-icon">⌕</span>
        <input
          type="text"
          placeholder="Search gene symbol, e.g. TAGLN"
          value={query}
          onChange={e => setQuery(e.target.value)}
          onKeyDown={onKeyDown}
          onFocus={() => suggestions.length > 0 && setOpen(true)}
          onBlur={() => setTimeout(() => setOpen(false), 150)}
          autoComplete="off"
          spellCheck={false}
        />
      </div>
      {open && (
        <ul className="suggestions">
          {suggestions.map((m, i) => (
            <li
              key={m.symbol}
              className={i === highlighted ? "active" : ""}
              onMouseDown={() => commit(m.symbol)}
            >
              <strong>{m.symbol}</strong>
              <span className="gene-name">{m.name}</span>
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}
