"""Gene-level ASE testing from per-SNP allele counts.

Collapses count records to one (variable, fixed) pair *per replicate per gene*,
then applies the replicate-aware tests from `stats`. The replicate is the unit
of inference; SNPs within a gene are summed (not treated as independent), and
the two kunthii backgrounds are summarised separately for a consistency check.
"""
from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, List, Tuple

from . import stats as st
from .stats import cont_ratio

GENE_COLS = [
    "gene_id", "gene_name", "n_snps", "n_reps",
    "variable_reads", "fixed_reads", "mean_variable_ratio", "log2_ratio",
    "null_p", "t_stat", "p_ttest",
    "bb_mu", "bb_rho", "p_betabinom",
    "p_pooled_binom_anticons", "p_primary", "q_value",
    "direction", "per_background_log2", "consistent_backgrounds", "ase_call",
]


def _gene_null(records: List[Dict]) -> float:
    """Depth-weighted mean per-SNP null expectation across a gene's records."""
    num = sum(r["null_p"] * r["total_depth"] for r in records)
    den = sum(r["total_depth"] for r in records)
    return num / den if den else 0.5


def _replicate_counts(records: List[Dict]) -> Dict[str, Tuple[int, int, str]]:
    """sample -> (variable_sum, total_sum, background) across a gene's SNPs."""
    agg: Dict[str, List] = defaultdict(lambda: [0, 0, None])
    for r in records:
        a = agg[r["f1_sample"]]
        a[0] += r["variable_count"]
        a[1] += r["variable_count"] + r["fixed_count"]
        a[2] = r["background"]
    return {s: (v, n, bg) for s, (v, n, bg) in agg.items()}


def _background_log2(rep_counts: Dict[str, Tuple[int, int, str]]) -> Dict[str, float]:
    by_bg: Dict[str, List[float]] = defaultdict(list)
    for _, (v, n, bg) in rep_counts.items():
        if n > 0:
            by_bg[bg].append(cont_ratio(v, n))
    out: Dict[str, float] = {}
    for bg, ratios in by_bg.items():
        m = sum(ratios) / len(ratios)
        out[bg] = math.log2(m / (1 - m))
    return out


def test_genes(count_records: List[Dict], *,
               alpha: float = 0.05,
               min_effect_log2: float = 0.0,
               min_reps: int = 2) -> List[Dict]:
    """Return one gene-level ASE record per gene (excluding intergenic)."""
    by_gene: Dict[str, List[Dict]] = defaultdict(list)
    for r in count_records:
        if r.get("gene_id") in (None, "", "intergenic"):
            continue
        by_gene[r["gene_id"]].append(r)

    results: List[Dict] = []
    for gene_id, recs in by_gene.items():
        rep_counts = _replicate_counts(recs)
        null_p = _gene_null(recs)
        pairs = [(v, n) for (v, n, _) in rep_counts.values()]

        rep_test = st.replicate_logit_test(pairs, null_p=null_p)
        bb = st.betabinom_test(pairs, null_p=null_p)

        # Pooled binomial — descriptive only, flagged anti-conservative.
        v_tot = sum(v for v, _ in pairs)
        n_tot = sum(n for _, n in pairs)
        p_pooled = st.binomial_p(v_tot, n_tot, null_p) if n_tot else 1.0
        mean_ratio = (v_tot / n_tot) if n_tot else 0.5
        log2_ratio = (math.log2(mean_ratio / (1 - mean_ratio))
                      if 0 < mean_ratio < 1 else float("nan"))

        bg_log2 = _background_log2(rep_counts)

        # Primary p-value: replicate t-test if enough reps, else beta-binomial.
        if rep_test is not None:
            p_primary = rep_test.p_value
        elif bb is not None:
            p_primary = bb.p_value
        else:
            p_primary = 1.0

        results.append({
            "gene_id": gene_id,
            "gene_name": recs[0].get("gene_name", gene_id),
            "n_snps": len({r["snp_id"] for r in recs}),
            "n_reps": len(pairs),
            "variable_reads": v_tot,
            "fixed_reads": n_tot - v_tot,
            "mean_variable_ratio": round(mean_ratio, 4),
            "log2_ratio": round(log2_ratio, 4) if not math.isnan(log2_ratio) else "NA",
            "null_p": round(null_p, 4),
            "t_stat": round(rep_test.t_stat, 4) if rep_test else "NA",
            "p_ttest": rep_test.p_value if rep_test else "NA",
            "bb_mu": round(bb.mu, 4) if bb else "NA",
            "bb_rho": round(bb.rho, 4) if bb else "NA",
            "p_betabinom": bb.p_value if bb else "NA",
            "p_pooled_binom_anticons": p_pooled,
            "p_primary": p_primary,
            "q_value": 1.0,                     # filled after BH below
            "direction": st.direction(log2_ratio if not math.isnan(log2_ratio) else 0.0),
            "per_background_log2": ";".join(f"{k}={v:.3f}"
                                            for k, v in sorted(bg_log2.items())),
            "consistent_backgrounds": st.consistent_across(bg_log2),
            "_bg_log2": bg_log2,
            "_min_reps": min_reps,
            "ase_call": False,                  # filled below
        })

    # BH across genes on the primary p-value.
    pvals = [r["p_primary"] for r in results]
    qvals = st.bh_adjust(pvals)
    for r, q in zip(results, qvals):
        r["q_value"] = q
        effect_ok = (r["log2_ratio"] != "NA"
                     and abs(r["log2_ratio"]) >= min_effect_log2)
        r["ase_call"] = bool(q < alpha and effect_ok
                             and r["consistent_backgrounds"]
                             and r["n_reps"] >= r["_min_reps"])
        r.pop("_bg_log2", None)
        r.pop("_min_reps", None)
    results.sort(key=lambda r: r["p_primary"])
    return results
