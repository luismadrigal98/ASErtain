"""Gene-level ASE testing with nested replication and an n-honest combine rule.

Hierarchy:  background -> F1 plant (biological replicate) -> flower (technical).

Procedure per gene:
  1. Sum each plant's flowers into per-SNP (variable, fixed) counts.
  2. Test EACH plant on its own (beta-binomial across its SNPs when >= MIN_SNPS,
     else a pooled binomial). This is a per-individual ASE test.
  3. Combine plants by intersection–union (max-p): the gene's p is the largest
     per-plant p, so a call requires EVERY contributing plant to be individually
     significant — valid and conservative at small n, and it bakes in the
     cross-background consistency requirement (audit C2/C3/M4).

Flags surfaced for honesty under a single-parent reference (audit M3, M5):
  * fixed_allele_seen / n_plants_fixed_seen — distinguishes real complete ASE
    from an allele lost to reference-mapping bias.
  * possible_ref_bias — variable==reference and the fixed allele never observed.
  * phase_concordant — whether the gene's per-SNP directions agree within plants
    (a `phased`-tier mis-call would otherwise contaminate the summed ratio).
"""
from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from . import stats as st
from .stats import cont_ratio

MIN_SNPS_FOR_BETABINOM = 3

GENE_COLS = [
    "gene_id", "gene_name", "n_snps", "n_plants", "n_backgrounds", "n_flowers",
    "variable_reads", "fixed_reads", "mean_variable_ratio", "log2_ratio",
    "null_p", "p_primary", "q_value", "method",
    "per_plant_log2", "per_plant_p", "direction",
    "consistent_backgrounds", "phase_concordant",
    "fixed_allele_seen", "n_plants_fixed_seen", "possible_ref_bias",
    "n_both_hom_snps", "ase_call",
]


def _gene_null(records: List[Dict]) -> float:
    """Allelic-depth-weighted mean per-SNP null expectation."""
    num = sum(r["null_p"] * (r["variable_count"] + r["fixed_count"]) for r in records)
    den = sum(r["variable_count"] + r["fixed_count"] for r in records)
    return num / den if den else 0.5


def _plant_snp_counts(records: List[Dict]) -> Dict[str, List[Tuple[int, int]]]:
    """plant -> list of (variable, total) per SNP, summing flowers within plant."""
    per: Dict[str, Dict[str, List[int]]] = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in records:
        slot = per[r["plant"]][r["snp_id"]]
        slot[0] += r["variable_count"]
        slot[1] += r["variable_count"] + r["fixed_count"]
    return {p: [(v, n) for v, n in snps.values()] for p, snps in per.items()}


def _plant_test(pairs: List[Tuple[int, int]], null_p: float) -> Optional[Dict]:
    """Per-plant test over its SNP (k, n) pairs. Returns p, log2, method, k, n."""
    pairs = [(k, n) for k, n in pairs if n > 0]
    if not pairs:
        return None
    k_tot = sum(k for k, _ in pairs)
    n_tot = sum(n for _, n in pairs)
    if len(pairs) >= MIN_SNPS_FOR_BETABINOM:
        bb = st.betabinom_test(pairs, null_p=null_p)
        if bb is not None and bb.converged:
            return {"p": bb.p_value, "log2": bb.log2_ratio,
                    "method": "betabinom", "k": k_tot, "n": n_tot}
    ratio = cont_ratio(k_tot, n_tot)
    return {"p": st.binomial_p(k_tot, n_tot, null_p),
            "log2": math.log2(ratio / (1 - ratio)),
            "method": "binomial", "k": k_tot, "n": n_tot}


def test_genes(count_records: List[Dict], *,
               alpha: float = 0.05,
               min_effect_log2: float = 0.0,
               min_plants: int = 2,
               ref_is_variable: bool = False) -> List[Dict]:
    by_gene: Dict[str, List[Dict]] = defaultdict(list)
    for r in count_records:
        if r.get("gene_id") in (None, "", "intergenic"):
            continue
        by_gene[r["gene_id"]].append(r)

    bg_of_plant = {r["plant"]: r["background"] for r in count_records}
    results: List[Dict] = []

    for gene_id, recs in by_gene.items():
        null_p = _gene_null(recs)
        plant_pairs = _plant_snp_counts(recs)

        plant_res: Dict[str, Dict] = {}
        for plant, pairs in plant_pairs.items():
            res = _plant_test(pairs, null_p)
            if res is not None:
                plant_res[plant] = res
        if not plant_res:
            continue

        # Per-plant direction (relative to the null) and effect.
        plant_dir = {p: st.direction(cont_ratio(r["k"], r["n"]), null_p)
                     for p, r in plant_res.items()}
        bg_log2 = {}
        for p, r in plant_res.items():
            bg_log2.setdefault(bg_of_plant[p], []).append(r["log2"])
        bg_mean_log2 = {bg: sum(v) / len(v) for bg, v in bg_log2.items()}

        n_backgrounds = len(bg_mean_log2)
        # Intersection–union: gene p = worst (max) per-plant p.
        p_iut = max(r["p"] for r in plant_res.values())
        methods = sorted({r["method"] for r in plant_res.values()})

        # Pooled effect for display.
        v_tot = sum(r["variable_count"] for r in recs)
        f_tot = sum(r["fixed_count"] for r in recs)
        n_tot = v_tot + f_tot
        mean_ratio = (v_tot / n_tot) if n_tot else null_p
        log2_ratio = (math.log2(mean_ratio / (1 - mean_ratio))
                      if 0 < mean_ratio < 1 else float("nan"))

        dirs = {d for d in plant_dir.values() if d != "balanced"}
        consistent = (n_backgrounds >= 2 and len(dirs) == 1)

        # Within-plant per-SNP direction concordance (audit M5).
        phase_conc = _phase_concordant(recs, null_p)

        n_plants_fixed = sum(1 for r in plant_res.values()
                             if (r["n"] - r["k"]) > 0)
        fixed_seen = f_tot > 0
        possible_ref_bias = bool(ref_is_variable and not fixed_seen)

        results.append({
            "gene_id": gene_id,
            "gene_name": recs[0].get("gene_name", gene_id),
            "n_snps": len({r["snp_id"] for r in recs}),
            "n_plants": len(plant_res),
            "n_backgrounds": n_backgrounds,
            "n_flowers": len({r["flower"] for r in recs}),
            "variable_reads": v_tot,
            "fixed_reads": f_tot,
            "mean_variable_ratio": round(mean_ratio, 4),
            "log2_ratio": round(log2_ratio, 4) if not math.isnan(log2_ratio) else "NA",
            "null_p": round(null_p, 4),
            "p_primary": p_iut,
            "q_value": 1.0,
            "method": "+".join(methods),
            "per_plant_log2": ";".join(f"{p}={r['log2']:.3f}"
                                       for p, r in sorted(plant_res.items())),
            "per_plant_p": ";".join(f"{p}={r['p']:.2e}"
                                    for p, r in sorted(plant_res.items())),
            "direction": st.direction(mean_ratio, null_p),
            "consistent_backgrounds": consistent,
            "phase_concordant": phase_conc,
            "fixed_allele_seen": fixed_seen,
            "n_plants_fixed_seen": n_plants_fixed,
            "possible_ref_bias": possible_ref_bias,
            "n_both_hom_snps": len({r["snp_id"] for r in recs
                                    if r.get("tier") == "both_hom"}),
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
                             and r["n_plants"] >= r["_min_plants"]
                             and r["phase_concordant"])
        r.pop("_min_plants", None)
    results.sort(key=lambda r: r["p_primary"])
    return results


def _phase_concordant(records: List[Dict], null_p: float) -> bool:
    """Within each plant, do the gene's per-SNP imbalance directions agree?

    A `phased`-tier SNP whose parent was mis-called homozygous (parental ASE)
    can be summed in the wrong orientation; if a plant's SNPs disagree on
    direction the gene ratio is suspect. Returns True if every plant's SNPs are
    directionally concordant (ignoring near-balanced SNPs), else False.
    """
    per_plant_snp: Dict[str, Dict[str, List[int]]] = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in records:
        slot = per_plant_snp[r["plant"]][r["snp_id"]]
        slot[0] += r["variable_count"]
        slot[1] += r["variable_count"] + r["fixed_count"]
    for snps in per_plant_snp.values():
        dirs = set()
        for v, n in snps.values():
            if n < 1:
                continue
            d = st.direction(cont_ratio(v, n), null_p, tol=0.1)
            if d != "balanced":
                dirs.add(d)
        if len(dirs) > 1:
            return False
    return True
