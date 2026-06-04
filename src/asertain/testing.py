"""Gene-level ASE testing with nested replication (flower -> plant).

Replication hierarchy:  background -> F1 plant -> flower (RNA sample).
Flowers from one plant share a genome, so they are NOT independent biological
replicates. We therefore collapse each plant's flowers into one (variable, fixed)
count per gene, and treat the **plant** as the unit. The primary test is a
beta-binomial across plants, whose overdispersion parameter absorbs plant-to-
plant biological variation; a per-plant logit t-test is reported alongside.

RNA-only caveat surfaced as a column: when the reference is one parent, the
fixed-lineage allele can be lost to mapping bias, which masquerades as "complete
ASE". `fixed_allele_seen` / `n_plants_fixed_seen` let you tell a real complete
imbalance from an allele that was never observable — pair with --bias-mode
nmask/wasp before trusting extreme ratios.
"""
from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, List, Tuple

from . import stats as st
from .stats import cont_ratio

GENE_COLS = [
    "gene_id", "gene_name", "n_snps", "n_plants", "n_flowers",
    "variable_reads", "fixed_reads", "mean_variable_ratio", "log2_ratio",
    "null_p", "p_betabinom", "bb_rho", "t_stat", "p_ttest",
    "p_pooled_binom_anticons", "p_primary", "q_value",
    "direction", "per_background_log2", "consistent_backgrounds",
    "fixed_allele_seen", "n_plants_fixed_seen", "n_both_hom_snps", "ase_call",
]


def _gene_null(records: List[Dict]) -> float:
    num = sum(r["null_p"] * r["total_depth"] for r in records)
    den = sum(r["total_depth"] for r in records)
    return num / den if den else 0.5


def _plant_counts(records: List[Dict]) -> Dict[str, Tuple[int, int, str]]:
    """plant -> (variable_sum, total_sum, background), summing across flowers+SNPs."""
    agg: Dict[str, List] = defaultdict(lambda: [0, 0, None])
    for r in records:
        a = agg[r["plant"]]
        a[0] += r["variable_count"]
        a[1] += r["variable_count"] + r["fixed_count"]
        a[2] = r["background"]
    return {p: (v, n, bg) for p, (v, n, bg) in agg.items()}


def _background_log2(plant_counts: Dict[str, Tuple[int, int, str]]) -> Dict[str, float]:
    by_bg: Dict[str, List[float]] = defaultdict(list)
    for _, (v, n, bg) in plant_counts.items():
        if n > 0:
            by_bg[bg].append(cont_ratio(v, n))
    return {bg: math.log2((m := sum(r) / len(r)) / (1 - m))
            for bg, r in by_bg.items()}


def test_genes(count_records: List[Dict], *,
               alpha: float = 0.05,
               min_effect_log2: float = 0.0,
               min_plants: int = 2) -> List[Dict]:
    by_gene: Dict[str, List[Dict]] = defaultdict(list)
    for r in count_records:
        if r.get("gene_id") in (None, "", "intergenic"):
            continue
        by_gene[r["gene_id"]].append(r)

    results: List[Dict] = []
    for gene_id, recs in by_gene.items():
        plant_counts = _plant_counts(recs)
        null_p = _gene_null(recs)
        pairs = [(v, n) for (v, n, _) in plant_counts.values()]

        bb = st.betabinom_test(pairs, null_p=null_p)
        rep = st.replicate_logit_test(pairs, null_p=null_p)

        v_tot = sum(v for v, _ in pairs)
        n_tot = sum(n for _, n in pairs)
        f_tot = n_tot - v_tot
        p_pooled = st.binomial_p(v_tot, n_tot, null_p) if n_tot else 1.0
        mean_ratio = (v_tot / n_tot) if n_tot else 0.5
        log2_ratio = (math.log2(mean_ratio / (1 - mean_ratio))
                      if 0 < mean_ratio < 1 else float("nan"))

        bg_log2 = _background_log2(plant_counts)
        n_plants_fixed = sum(1 for (v, n, _) in plant_counts.values() if (n - v) > 0)

        if bb is not None:
            p_primary = bb.p_value
        elif rep is not None:
            p_primary = rep.p_value
        else:
            p_primary = 1.0

        results.append({
            "gene_id": gene_id,
            "gene_name": recs[0].get("gene_name", gene_id),
            "n_snps": len({r["snp_id"] for r in recs}),
            "n_plants": len(pairs),
            "n_flowers": len({r["flower"] for r in recs}),
            "variable_reads": v_tot,
            "fixed_reads": f_tot,
            "mean_variable_ratio": round(mean_ratio, 4),
            "log2_ratio": round(log2_ratio, 4) if not math.isnan(log2_ratio) else "NA",
            "null_p": round(null_p, 4),
            "p_betabinom": bb.p_value if bb else "NA",
            "bb_rho": round(bb.rho, 4) if bb else "NA",
            "t_stat": round(rep.t_stat, 4) if rep else "NA",
            "p_ttest": rep.p_value if rep else "NA",
            "p_pooled_binom_anticons": p_pooled,
            "p_primary": p_primary,
            "q_value": 1.0,
            "direction": st.direction(log2_ratio if not math.isnan(log2_ratio) else 0.0),
            "per_background_log2": ";".join(f"{k}={v:.3f}" for k, v in sorted(bg_log2.items())),
            "consistent_backgrounds": st.consistent_across(bg_log2),
            "fixed_allele_seen": f_tot > 0,
            "n_plants_fixed_seen": n_plants_fixed,
            "n_both_hom_snps": len({r["snp_id"] for r in recs if r.get("tier") == "both_hom"}),
            "_min_plants": min_plants,
            "ase_call": False,
        })

    qvals = st.bh_adjust([r["p_primary"] for r in results])
    for r, q in zip(results, qvals):
        r["q_value"] = q
        effect_ok = (r["log2_ratio"] != "NA"
                     and abs(r["log2_ratio"]) >= min_effect_log2)
        r["ase_call"] = bool(q < alpha and effect_ok
                             and r["consistent_backgrounds"]
                             and r["n_plants"] >= r["_min_plants"])
        r.pop("_min_plants", None)
    results.sort(key=lambda r: r["p_primary"])
    return results
