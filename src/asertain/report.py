"""Lightweight reporting: a self-contained HTML summary (+ optional plots).

Working scaffold — renders the gene-level ASE table and headline counts to HTML
with no third-party templating. If matplotlib is importable, it also writes a
per-gene allelic-ratio scatter; otherwise it degrades gracefully.
"""
from __future__ import annotations

import html
import os
from typing import Dict, List, Optional

from .labels import Labels
from .tables import read_table


def _summary_counts(genes: List[Dict]) -> Dict[str, int]:
    calls = sum(1 for g in genes if g.get("ase_call") in ("True", True))
    var = sum(1 for g in genes if g.get("direction") == "variable"
              and g.get("ase_call") in ("True", True))
    fix = sum(1 for g in genes if g.get("direction") == "fixed"
              and g.get("ase_call") in ("True", True))
    return {"genes_tested": len(genes), "ase_genes": calls,
            "variable_biased": var, "fixed_biased": fix}


def _maybe_plot(genes: List[Dict], path: str,
                labels: Labels = Labels()) -> Optional[str]:
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
    ax.set_xlabel(f"log2 allelic ratio ({labels.variable} / {labels.fixed})")
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
                 title: str = "ASErtain report",
                 labels: Optional[Labels] = None) -> str:
    labels = labels or Labels()
    genes = read_table(gene_ase_tsv)
    counts = _summary_counts(genes)
    plot_path = os.path.join(os.path.dirname(out_html) or ".",
                             "ase_volcano.png")
    plotted = _maybe_plot(genes, plot_path, labels)

    def _dir_label(g: Dict) -> str:
        return labels.value_to_display(str(g.get("direction", "")))

    def _concordant(g: Dict) -> str:
        # The directional-concordance gate that actually decides the call differs
        # by aggregation: snp_concordant (maxsnp) vs phase_concordant (plant).
        key = ("snp_concordant" if str(g.get("agg_method", "")) == "maxsnp"
               else "phase_concordant")
        return g.get(key, "")

    rows = "\n".join(
        "<tr>" + "".join(f"<td>{html.escape(str(val))}</td>" for val in (
            g.get("gene_id", ""), g.get("gene_name", ""), g.get("n_snps", ""),
            g.get("n_plants", ""), g.get("log2_ratio", ""), g.get("p_primary", ""),
            g.get("q_value", ""), _dir_label(g), g.get("consistent_backgrounds", ""),
            _concordant(g), g.get("low_ambiguity", ""),
            g.get("ambiguous_fraction", ""), g.get("n_snps_same_dir", ""),
            g.get("fixed_allele_seen", ""), g.get("ase_call", ""))) + "</tr>"
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
<p class="k">Lineages: variable = <b>{html.escape(labels.variable)}</b>,
fixed = <b>{html.escape(labels.fixed)}</b></p>
<p class="k">Genes tested: {counts['genes_tested']} &nbsp;|&nbsp;
ASE genes: <b>{counts['ase_genes']}</b> &nbsp;|&nbsp;
{html.escape(labels.variable)}-biased: {counts['variable_biased']} &nbsp;|&nbsp;
{html.escape(labels.fixed)}-biased: {counts['fixed_biased']}</p>
{img}
<h2>Per-gene results (first 500)</h2>
<p class="k">An <b>ASE</b> call requires ALL of:
q &lt; alpha &nbsp;·&nbsp; |log2| &ge; min-effect &nbsp;·&nbsp;
<b>consistent</b> (same direction across backgrounds) &nbsp;·&nbsp;
nPlant &ge; min-plants &nbsp;·&nbsp;
<b>concordant</b> (per-SNP directions agree: snp_concordant for maxsnp,
phase_concordant for the plant/haplotype path) &nbsp;·&nbsp;
<b>lowAmbig</b> (ambiguous-read fraction &le; ceiling).
<b>{html.escape(labels.fixed)}Seen</b> is provenance only — it is NOT part of the
gate. A gene that is FALSE despite a small q fails one of the bold columns.</p>
<table><tr><th>gene_id</th><th>name</th><th>nSNP</th><th>nPlant</th>
<th>log2({html.escape(labels.variable)}/{html.escape(labels.fixed)})</th>
<th>p</th><th>q</th><th>dir</th><th>consistent</th><th>concordant</th>
<th>lowAmbig</th><th>ambigFrac</th><th>nSameDir</th>
<th>{html.escape(labels.fixed)}Seen</th><th>ASE</th></tr>
{rows}</table></body></html>"""
    with open(out_html, "w") as fh:
        fh.write(doc)
    return out_html
