"""Lightweight reporting: a self-contained HTML summary (+ optional plots).

Working scaffold — renders the gene-level ASE table and headline counts to HTML
with no third-party templating. If matplotlib is importable, it also writes a
per-gene allelic-ratio scatter; otherwise it degrades gracefully.
"""
from __future__ import annotations

import html
import os
from typing import Dict, List, Optional

from .tables import read_table


def _summary_counts(genes: List[Dict]) -> Dict[str, int]:
    calls = sum(1 for g in genes if g.get("ase_call") in ("True", True))
    var = sum(1 for g in genes if g.get("direction") == "variable"
              and g.get("ase_call") in ("True", True))
    fix = sum(1 for g in genes if g.get("direction") == "fixed"
              and g.get("ase_call") in ("True", True))
    return {"genes_tested": len(genes), "ase_genes": calls,
            "variable_biased": var, "fixed_biased": fix}


def _maybe_plot(genes: List[Dict], path: str) -> Optional[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    xs, ys, cols = [], [], []
    for g in genes:
        try:
            lr = float(g["log2_ratio"])
            q = float(g["q_value"])
        except (ValueError, KeyError, TypeError):
            continue
        xs.append(lr)
        ys.append(-_safe_log10(q))
        cols.append("crimson" if g.get("ase_call") in ("True", True) else "grey")
    if not xs:
        return None
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(xs, ys, c=cols, s=12, alpha=0.7)
    ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("log2 allelic ratio (variable / fixed)")
    ax.set_ylabel("-log10 q-value")
    ax.set_title("F1 allele-specific expression")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def _safe_log10(x: float) -> float:
    import math
    return math.log10(x) if x > 0 else -300.0


def write_report(gene_ase_tsv: str, out_html: str, *,
                 title: str = "ASErtain report") -> str:
    genes = read_table(gene_ase_tsv)
    counts = _summary_counts(genes)
    plot_path = os.path.join(os.path.dirname(out_html) or ".",
                             "ase_volcano.png")
    plotted = _maybe_plot(genes, plot_path)

    rows = "\n".join(
        "<tr>" + "".join(f"<td>{html.escape(str(g.get(c, '')))}</td>"
                         for c in ("gene_id", "gene_name", "n_snps", "n_plants",
                                   "log2_ratio", "p_primary", "q_value",
                                   "direction", "consistent_backgrounds",
                                   "fixed_allele_seen", "ase_call")) + "</tr>"
        for g in genes[:500]
    )
    img = (f'<img src="{os.path.basename(plotted)}" style="max-width:500px">'
           if plotted else "<p><em>matplotlib not available — plot skipped.</em></p>")

    doc = f"""<!doctype html><html><head><meta charset="utf-8">
<title>{html.escape(title)}</title>
<style>body{{font-family:system-ui,sans-serif;margin:2rem;}}
table{{border-collapse:collapse;font-size:13px}}
td,th{{border:1px solid #ccc;padding:2px 6px}}
th{{background:#f0f0f0}} .k{{color:#666}}</style></head><body>
<h1>{html.escape(title)}</h1>
<p class="k">Genes tested: {counts['genes_tested']} &nbsp;|&nbsp;
ASE genes: <b>{counts['ase_genes']}</b> &nbsp;|&nbsp;
variable-biased: {counts['variable_biased']} &nbsp;|&nbsp;
fixed-biased: {counts['fixed_biased']}</p>
{img}
<h2>Per-gene results (first 500)</h2>
<table><tr><th>gene_id</th><th>name</th><th>nSNP</th><th>nPlant</th>
<th>log2</th><th>p</th><th>q</th><th>dir</th><th>consistent</th>
<th>fixedSeen</th><th>ASE</th></tr>
{rows}</table></body></html>"""
    with open(out_html, "w") as fh:
        fh.write(doc)
    return out_html
