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

Flowers contribute unequally (a deeply sequenced flower would otherwise dominate
its plant's pooled allelic ratio). Before summing, each flower is rescaled by a
per-plant **size factor** so every flower of a plant contributes comparably —
see `flower_size_factors`. This is on by default (`flower_norm="equalize"`) and
preserves each flower's own allelic ratio (both alleles scale together).

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
FLOWER_NORM_MODES = ("equalize", "none")

GENE_COLS = [
    "gene_id", "gene_name", "n_snps", "n_plants", "n_backgrounds", "n_flowers",
    "variable_reads", "fixed_reads", "other_reads", "mean_variable_ratio",
    "log2_ratio", "null_p", "p_primary", "q_value", "method",
    "per_plant_log2", "per_plant_p", "direction",
    "consistent_backgrounds", "phase_concordant",
    "ambiguous_fraction", "low_ambiguity",
    "variable_allele_seen", "n_plants_variable_seen",
    "fixed_allele_seen", "n_plants_fixed_seen", "possible_ref_bias",
    "n_both_hom_snps", "ase_call",
    "agg_method", "top_snp", "top_snp_p", "n_snps_same_dir",
    "snp_concordant", "within_gene_correction",
]

# Default ceiling on the fraction of allele-overlapping reads that match neither
# clean haplotype/allele. For --counter haplotype these are *ambiguous* fragments
# (carrying BOTH a variable and a fixed allele) — a direct phasing-quality signal;
# for --counter pileup they are third-allele / error reads. A gene above the
# ceiling is flagged (low_ambiguity=False) and not called.
MAX_OTHER_FRACTION = 0.10

# Aggregation strategies for turning a gene's SNPs into one call.
#   "plant"  : the default/haplotype path -- pool a plant's SNPs, test each
#              plant, combine plants (LD handled by the haplotype counter or
#              absorbed as dispersion by the beta-binomial).
#   "maxsnp" : the advisor's alternative -- a plain binomial PER SNP, take the
#              strongest-signal SNP per gene, and require the gene's SNPs to
#              agree in direction (all point to the same parent). Avoids LD by
#              never combining correlated SNPs; the selection of the best of m
#              SNPs is corrected with WITHIN_GENE_CORRECTION below.
GENE_AGGREGATION_MODES = ("plant", "maxsnp")
WITHIN_GENE_CORRECTION_MODES = ("sidak", "bonferroni", "none")


def _gene_null(records: List[Dict]) -> float:
    """Allelic-depth-weighted mean per-SNP null expectation."""
    num = sum(r["null_p"] * (r["variable_count"] + r["fixed_count"]) for r in records)
    den = sum(r["variable_count"] + r["fixed_count"] for r in records)
    return num / den if den else 0.5


# ---------------------------------------------------------------------------
# Flower-contribution normalisation (audit/task: technical replicates differ in
# depth, so a single deep flower must not dominate its plant's allelic ratio)
# ---------------------------------------------------------------------------

def flower_size_factors(records: List[Dict], *,
                        mode: str = "equalize") -> Dict[Tuple[str, str], float]:
    """Per-(plant, flower) size factors that equalise flower contribution.

    A flower's *library* depth is summed across **all** SNPs of its plant (so
    the factor is a stable library-level quantity, not a noisy per-gene one) and
    divided by the geometric mean depth of that plant's flowers. Dividing a
    flower's counts by this factor rescales every flower of a plant to the same
    effective depth, so deep and shallow flowers contribute comparably. Because
    both alleles are divided by the *same* factor, each flower's allelic ratio is
    untouched — only its weight in the pooled (k, n) changes.

    `mode="none"` returns an empty map (callers then treat every factor as 1.0,
    i.e. the previous raw-summing behaviour).

    Honesty note: equalising *weight* means each flower's ratio is trusted
    equally regardless of its depth — a deep flower is down-weighted, a shallow
    one up-weighted. That is the intended behaviour (a flower is an observation,
    not a depth quota) but it reweights statistical confidence rather than
    keeping it strictly depth-proportional. The plant remains the unit of
    inference, and the across-plant max-p rule plus the cross-background
    consistency requirement are the real guards against any single noisy flower.
    """
    if mode == "none":
        return {}
    if mode not in FLOWER_NORM_MODES:
        raise ValueError(f"unknown flower_norm mode '{mode}'; "
                         f"choose from {FLOWER_NORM_MODES}")
    depth: Dict[Tuple[str, str], float] = defaultdict(float)
    flowers_of_plant: Dict[str, set] = defaultdict(set)
    for r in records:
        d = r["variable_count"] + r["fixed_count"]
        depth[(r["plant"], r["flower"])] += d
        flowers_of_plant[r["plant"]].add(r["flower"])

    factors: Dict[Tuple[str, str], float] = {}
    for plant, flowers in flowers_of_plant.items():
        positive = [depth[(plant, f)] for f in flowers if depth[(plant, f)] > 0]
        if len(positive) < 2:
            # One (informative) flower: nothing to equalise.
            for f in flowers:
                factors[(plant, f)] = 1.0
            continue
        gm = math.exp(sum(math.log(d) for d in positive) / len(positive))
        for f in flowers:
            d = depth[(plant, f)]
            factors[(plant, f)] = (d / gm) if d > 0 else 1.0
    return factors


def _plant_snp_counts(records: List[Dict],
                      size_factors: Optional[Dict[Tuple[str, str], float]] = None,
                      *, with_null: bool = False):
    """plant -> per-SNP allele counts, summing size-factor-rescaled flowers.

    With `size_factors=None` (or empty) this is the plain raw sum. By default
    each SNP is a (variable, total) pair. With `with_null=True` it is a
    (variable, total, null_p) triple — the per-SNP null expectation, needed to
    test SNPs against their OWN nulls when those differ across a gene (audit C1,
    --bias-mode null-shift). The null is constant for a (plant, SNP) across its
    flowers, so the first record's value is used.

    The flower size factor only ever *reduces or preserves* a plant's pooled
    total N: by AM>=GM, sum_f(depth_f / w_f) = F * geomean(depths) <= sum_f
    depth_f. So equalising flower weight cannot inflate the per-plant N above the
    reads actually observed (it reweights each flower's ratio, not the test's
    total confidence)."""
    sf = size_factors or {}
    per: Dict[str, Dict[str, List[float]]] = defaultdict(
        lambda: defaultdict(lambda: [0.0, 0.0]))
    null_of: Dict[Tuple[str, str], float] = {}
    for r in records:
        w = sf.get((r["plant"], r["flower"]), 1.0) or 1.0
        slot = per[r["plant"]][r["snp_id"]]
        slot[0] += r["variable_count"] / w
        slot[1] += (r["variable_count"] + r["fixed_count"]) / w
        null_of.setdefault((r["plant"], r["snp_id"]), r.get("null_p", 0.5))
    # Round the rescaled effective counts back to integers for the
    # binomial / beta-binomial likelihoods. Rounding is monotone so v<=n holds.
    if not with_null:
        return {p: [(int(round(v)), int(round(n))) for v, n in snps.values()]
                for p, snps in per.items()}
    return {p: [(int(round(v)), int(round(n)), null_of[(p, sid)])
                for sid, (v, n) in snps.items()]
            for p, snps in per.items()}


def _plant_test(triples: List[Tuple[int, int, float]]) -> Optional[Dict]:
    """Per-plant test over its SNP (k, n, null_p) triples.

    Returns p, log2, dir, method, k, n, n_snps, rho (rho is None for the
    binomial). When the SNPs share a null (the usual case — null is 0.5 for every
    mode except a per-SNP null-shift control), the counts are pooled and tested
    against that common null with the validated beta-binomial / exact binomial.
    When the nulls DIFFER across SNPs, pooling would be invalid (audit C1), so a
    per-SNP-null logit-shift LRT is used instead (`stats.shift_test`).
    """
    triples = [(k, n, p0) for k, n, p0 in triples if n > 0]
    if not triples:
        return None
    k_tot = sum(k for k, _, _ in triples)
    n_tot = sum(n for _, n, _ in triples)
    n_snps = len(triples)
    nulls = [p0 for _, _, p0 in triples]
    homogeneous = (max(nulls) - min(nulls)) < 1e-9
    common_null = nulls[0]

    if not homogeneous:
        # Heterogeneous per-SNP nulls: never pool. Fit one logit shift.
        sh = st.shift_test(triples)
        if sh is not None:
            direction = ("variable" if sh.delta > 0
                         else "fixed" if sh.delta < 0 else "balanced")
            return {"p": sh.p_value, "log2": sh.log2_ratio, "dir": direction,
                    "method": sh.method, "k": k_tot, "n": n_tot,
                    "n_snps": n_snps, "rho": sh.rho}

    pairs = [(k, n) for k, n, _ in triples]
    if n_snps >= MIN_SNPS_FOR_BETABINOM:
        bb = st.betabinom_test(pairs, null_p=common_null)
        if bb is not None and bb.converged:
            return {"p": bb.p_value, "log2": bb.log2_ratio,
                    "dir": st.direction(bb.mu, common_null),
                    "method": "betabinom", "k": k_tot, "n": n_tot,
                    "n_snps": n_snps, "rho": bb.rho}
    ratio = cont_ratio(k_tot, n_tot)
    return {"p": st.binomial_p(k_tot, n_tot, common_null),
            "log2": math.log2(ratio / (1 - ratio)),
            "dir": st.direction(ratio, common_null),
            "method": "binomial", "k": k_tot, "n": n_tot,
            "n_snps": n_snps, "rho": None}


def _score_records(recs, *, size_factors, combine):
    """Statistical aggregate for a set of count records treated as ONE unit.

    The unit is a whole gene for ``gene_aggregation="plant"`` (the default /
    haplotype path) or a single SNP for ``gene_aggregation="maxsnp"``. In both
    cases the per-plant test + plant-combine machinery is identical — only the
    set of records handed in differs. Returns None if no plant is testable.
    """
    plant_pairs = _plant_snp_counts(recs, size_factors, with_null=True)
    plant_res: Dict[str, Dict] = {}
    for plant, triples in plant_pairs.items():
        res = _plant_test(triples)
        if res is not None:
            plant_res[plant] = res
    if not plant_res:
        return None

    null_p = _gene_null(recs)
    bg_of_plant = {r["plant"]: r["background"] for r in recs}
    plant_dir = {p: r["dir"] for p, r in plant_res.items()}
    bg_log2: Dict[str, List[float]] = {}
    for p, r in plant_res.items():
        bg_log2.setdefault(bg_of_plant[p], []).append(r["log2"])
    n_backgrounds = len(bg_log2)

    p_maxp = max(r["p"] for r in plant_res.values())
    if combine == "stouffer":
        p_primary_val = st.stouffer_combine(
            [(r["p"], plant_dir[p]) for p, r in plant_res.items()])
    else:
        p_primary_val = p_maxp
    methods = sorted({r["method"] for r in plant_res.values()})

    norm_k = sum(k for trips in plant_pairs.values() for k, _, _ in trips)
    norm_n = sum(n for trips in plant_pairs.values() for _, n, _ in trips)
    mean_ratio = (norm_k / norm_n) if norm_n else null_p
    log2_ratio = (math.log2(mean_ratio / (1 - mean_ratio))
                  if 0 < mean_ratio < 1 else float("nan"))

    dirs = {d for d in plant_dir.values() if d != "balanced"}
    consistent = (n_backgrounds >= 2 and len(dirs) == 1)
    n_plants_fixed = sum(1 for r in plant_res.values() if (r["n"] - r["k"]) > 0)
    n_plants_variable = sum(1 for r in plant_res.values() if r["k"] > 0)

    return {
        "p_primary": p_primary_val, "p_maxp": p_maxp,
        "mean_ratio": mean_ratio, "log2_ratio": log2_ratio,
        "direction": st.direction(mean_ratio, null_p),
        "method": "+".join(methods),
        "per_plant_log2": ";".join(f"{p}={r['log2']:.3f}"
                                   for p, r in sorted(plant_res.items())),
        "per_plant_p": ";".join(f"{p}={r['p']:.2e}"
                                for p, r in sorted(plant_res.items())),
        "n_plants": len(plant_res), "n_backgrounds": n_backgrounds,
        "consistent_backgrounds": consistent,
        "n_plants_variable_seen": n_plants_variable,
        "n_plants_fixed_seen": n_plants_fixed,
    }


def _gene_qc(recs, *, ref_lineage, max_other_fraction):
    """Gene-level QC / provenance, always computed over ALL the gene's records
    (independent of the aggregation strategy used for the call)."""
    null_p = _gene_null(recs)
    v_tot = sum(r["variable_count"] for r in recs)
    f_tot = sum(r["fixed_count"] for r in recs)
    o_tot = sum(r.get("other_count", 0) for r in recs)
    total_obs = v_tot + f_tot + o_tot
    other_fraction = (o_tot / total_obs) if total_obs else 0.0
    fixed_seen = f_tot > 0
    variable_seen = v_tot > 0
    possible_ref_bias = bool(
        (ref_lineage == "variable" and not fixed_seen)
        or (ref_lineage == "fixed" and not variable_seen))
    hap_counts = [r["n_hap_snps"] for r in recs
                  if r.get("n_hap_snps") not in (None, "")]
    n_snps = max(hap_counts) if hap_counts else len({r["snp_id"] for r in recs})
    return {
        "null_p": null_p, "variable_reads": v_tot, "fixed_reads": f_tot,
        "other_reads": o_tot, "ambiguous_fraction": other_fraction,
        "low_ambiguity": other_fraction <= max_other_fraction,
        "variable_allele_seen": variable_seen, "fixed_allele_seen": fixed_seen,
        "possible_ref_bias": possible_ref_bias,
        "phase_concordant": _phase_concordant(recs),
        "n_snps": n_snps,
        "n_both_hom_snps": len({r["snp_id"] for r in recs
                                if r.get("tier") == "both_hom"}),
    }


def _aggregate_maxsnp(recs, *, size_factors, combine, correction):
    """The advisor's alternative: a plain binomial PER SNP, then take the
    strongest-signal SNP for the gene and validate that the gene's SNPs all
    point to the same parent.

    Rationale. SNPs in one gene are in LD, so combining them (pileup pooling)
    pseudo-replicates the same molecule and inflates significance. Rather than
    model that correlation, we sidestep it: each SNP is tested on its own
    (one independent binomial per SNP, replicate-aware via the same per-plant
    max-p combine), and the gene takes its **best** SNP. Because we picked the
    best of ``m`` SNPs, the selected p is corrected for ``m`` (Sidak by
    default; under LD the SNPs are positively correlated so this is
    conservative). The validation the advisor asked for is the directional
    **concordance** requirement: every informative SNP of the gene must point
    the same way (same parent); a gene whose SNPs disagree is not called.
    """
    by_snp: Dict[str, List[Dict]] = defaultdict(list)
    for r in recs:
        by_snp[r["snp_id"]].append(r)

    snp_units: Dict[str, Dict] = {}
    for sid, srecs in by_snp.items():
        u = _score_records(srecs, size_factors=size_factors, combine=combine)
        if u is not None:
            snp_units[sid] = u
    if not snp_units:
        return None, None

    dirset = {u["direction"] for u in snp_units.values()
              if u["direction"] in ("variable", "fixed")}
    snp_concordant = len(dirset) <= 1
    agreed_dir = next(iter(dirset)) if len(dirset) == 1 else None

    # Strongest signal among the SNPs in the agreed direction (or all SNPs if
    # every SNP looks balanced / there is no agreed direction).
    candidates = {sid: u for sid, u in snp_units.items()
                  if agreed_dir is None or u["direction"] == agreed_dir}         or snp_units
    top_sid, top_u = min(candidates.items(), key=lambda kv: kv[1]["p_primary"])

    m = len(snp_units)
    p_raw = float(top_u["p_primary"])
    if correction == "sidak":
        p_gene = 1.0 - (1.0 - min(max(p_raw, 0.0), 1.0)) ** m
    elif correction == "bonferroni":
        p_gene = min(1.0, p_raw * m)
    else:
        p_gene = p_raw

    n_same_dir = sum(1 for u in snp_units.values()
                     if agreed_dir and u["direction"] == agreed_dir)
    unit = dict(top_u)
    unit["p_primary"] = p_gene
    extra = {
        "agg_method": "maxsnp", "top_snp": top_sid,
        "top_snp_p": round(p_raw, 6),
        "n_snps_same_dir": n_same_dir, "snp_concordant": snp_concordant,
        "within_gene_correction": correction,
    }
    return unit, extra


def test_genes(count_records: List[Dict], *,
               alpha: float = 0.05,
               min_effect_log2: float = 0.0,
               min_plants: int = 2,
               ref_is_variable: bool = False,
               ref_lineage: Optional[str] = None,
               flower_norm: str = "equalize",
               combine: str = "maxp",
               gene_aggregation: str = "plant",
               within_gene_correction: str = "sidak",
               max_other_fraction: float = MAX_OTHER_FRACTION) -> List[Dict]:
    # `ref_lineage` ('variable'/'fixed'/None) is the lineage the reference equals,
    # and drives the symmetric reference-bias flag. `ref_is_variable=True` is the
    # backward-compatible alias for ref_lineage='variable'.
    #
    # `gene_aggregation` selects how a gene's SNPs become one call:
    #   "plant"  (default) pool a plant's SNPs, test each plant, combine plants —
    #            the path the read-backed haplotype counter feeds (one (K,N)/gene).
    #   "maxsnp" plain binomial per SNP, take the strongest SNP, require the
    #            gene's SNPs to agree in direction. `within_gene_correction`
    #            ('sidak'|'bonferroni'|'none') corrects the picked-best-of-m p.
    if ref_lineage is None and ref_is_variable:
        ref_lineage = "variable"
    if combine not in ("maxp", "stouffer"):
        raise ValueError(f"unknown combine rule '{combine}'; choose maxp|stouffer")
    if gene_aggregation not in GENE_AGGREGATION_MODES:
        raise ValueError(f"unknown gene_aggregation '{gene_aggregation}'; "
                         f"choose from {GENE_AGGREGATION_MODES}")
    if within_gene_correction not in WITHIN_GENE_CORRECTION_MODES:
        raise ValueError(f"unknown within_gene_correction "
                         f"'{within_gene_correction}'; choose from "
                         f"{WITHIN_GENE_CORRECTION_MODES}")
    by_gene: Dict[str, List[Dict]] = defaultdict(list)
    for r in count_records:
        if r.get("gene_id") in (None, "", "intergenic"):
            continue
        by_gene[r["gene_id"]].append(r)

    size_factors = flower_size_factors(count_records, mode=flower_norm)
    results: List[Dict] = []

    for gene_id, recs in by_gene.items():
        qc = _gene_qc(recs, ref_lineage=ref_lineage,
                      max_other_fraction=max_other_fraction)
        if gene_aggregation == "maxsnp":
            unit, extra = _aggregate_maxsnp(
                recs, size_factors=size_factors, combine=combine,
                correction=within_gene_correction)
        else:
            unit = _score_records(recs, size_factors=size_factors,
                                  combine=combine)
            extra = {"agg_method": "plant", "top_snp": "NA", "top_snp_p": "NA",
                     "n_snps_same_dir": "NA",
                     "snp_concordant": qc["phase_concordant"],
                     "within_gene_correction": "none"}
        if unit is None:
            continue

        mean_ratio = unit["mean_ratio"]
        log2_ratio = unit["log2_ratio"]
        # The validation gate: the maxsnp path is gated on cross-SNP direction
        # concordance (the advisor's requirement); the plant path keeps the
        # within-plant phase concordance check.
        concordant_for_call = (extra["snp_concordant"]
                               if gene_aggregation == "maxsnp"
                               else qc["phase_concordant"])

        results.append({
            "gene_id": gene_id,
            "gene_name": recs[0].get("gene_name", gene_id),
            "n_snps": qc["n_snps"],
            "n_plants": unit["n_plants"],
            "n_backgrounds": unit["n_backgrounds"],
            "n_flowers": len({r["flower"] for r in recs}),
            "variable_reads": qc["variable_reads"],
            "fixed_reads": qc["fixed_reads"],
            "other_reads": qc["other_reads"],
            "mean_variable_ratio": round(mean_ratio, 4),
            "log2_ratio": round(log2_ratio, 4) if not math.isnan(log2_ratio) else "NA",
            "null_p": round(qc["null_p"], 4),
            "p_primary": unit["p_primary"],
            "q_value": 1.0,
            "method": unit["method"],
            "per_plant_log2": unit["per_plant_log2"],
            "per_plant_p": unit["per_plant_p"],
            "direction": unit["direction"],
            "consistent_backgrounds": unit["consistent_backgrounds"],
            "phase_concordant": qc["phase_concordant"],
            "ambiguous_fraction": round(qc["ambiguous_fraction"], 4),
            "low_ambiguity": qc["low_ambiguity"],
            "variable_allele_seen": qc["variable_allele_seen"],
            "n_plants_variable_seen": unit["n_plants_variable_seen"],
            "fixed_allele_seen": qc["fixed_allele_seen"],
            "n_plants_fixed_seen": unit["n_plants_fixed_seen"],
            "possible_ref_bias": qc["possible_ref_bias"],
            "n_both_hom_snps": qc["n_both_hom_snps"],
            "agg_method": extra["agg_method"],
            "top_snp": extra["top_snp"],
            "top_snp_p": extra["top_snp_p"],
            "n_snps_same_dir": extra["n_snps_same_dir"],
            "snp_concordant": extra["snp_concordant"],
            "within_gene_correction": extra["within_gene_correction"],
            "_min_plants": min_plants,
            "_concordant_for_call": concordant_for_call,
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
                             and r["_concordant_for_call"]
                             and r["low_ambiguity"])
        r.pop("_min_plants", None)
        r.pop("_concordant_for_call", None)
    results.sort(key=lambda r: r["p_primary"])
    return results


# ---------------------------------------------------------------------------
# Verbose audit tables (written when --verbose) — let you trace every number
# from raw per-flower counts up to the gene-level call.
# ---------------------------------------------------------------------------

SNP_DETAIL_COLS = [
    "gene_id", "gene_name", "chrom", "pos", "snp_id", "tier",
    "plant", "background", "variable_allele", "fixed_allele", "variable_is_ref",
    "n_flowers", "variable_count", "fixed_count", "other_count",
    "allelic_depth", "variable_ratio", "null_p",
]

# Per gene × SNP (plants and flowers collapsed) — the fundamental ASE evidence
# unit. Lets a reader verify the counts behind every gene without expanding to
# the full gene×SNP×plant table.
SNP_GENE_COLS = [
    "gene_id", "gene_name", "chrom", "pos", "snp_id", "tier",
    "n_plants", "n_flowers", "variable_count", "fixed_count", "other_count",
    "allelic_depth", "variable_ratio", "per_plant_counts", "null_p",
]

PLANT_DETAIL_COLS = [
    "gene_id", "gene_name", "plant", "background", "n_snps",
    "variable_reads", "fixed_reads", "allelic_depth",
    "variable_ratio", "log2_ratio", "rho", "method", "p_plant", "null_p",
]


def snp_plant_detail(count_records: List[Dict]) -> List[Dict]:
    """One row per (gene, SNP, plant): flowers summed, with the per-SNP ratio.

    This is the fully expanded evidence table — every gene×SNP×plant combination
    and its allele counts — for auditing exactly what drove each gene call.
    """
    by_gene: Dict[str, List[Dict]] = defaultdict(list)
    for r in count_records:
        if r.get("gene_id") in (None, "", "intergenic"):
            continue
        by_gene[r["gene_id"]].append(r)

    rows: List[Dict] = []
    for gene_id, recs in by_gene.items():
        null_p = _gene_null(recs)
        agg: Dict[Tuple[str, str], Dict] = {}
        for r in recs:
            key = (r["plant"], r["snp_id"])
            a = agg.setdefault(key, {"v": 0, "f": 0, "o": 0, "flowers": set(), "meta": r})
            a["v"] += r["variable_count"]
            a["f"] += r["fixed_count"]
            a["o"] += r["other_count"]
            a["flowers"].add(r["flower"])
        for (plant, snp_id), a in agg.items():
            m = a["meta"]
            depth = a["v"] + a["f"]
            rows.append({
                "gene_id": gene_id, "gene_name": m.get("gene_name", gene_id),
                "chrom": m["chrom"], "pos": m["pos"], "snp_id": snp_id,
                "tier": m.get("tier", ""), "plant": plant,
                "background": m["background"],
                "variable_allele": m["variable_allele"],
                "fixed_allele": m["fixed_allele"],
                "variable_is_ref": m["variable_is_ref"],
                "n_flowers": len(a["flowers"]),
                "variable_count": a["v"], "fixed_count": a["f"], "other_count": a["o"],
                "allelic_depth": depth,
                "variable_ratio": round(a["v"] / depth, 4) if depth else "NA",
                "null_p": round(null_p, 4),
            })
    rows.sort(key=lambda r: (r["gene_id"], r["chrom"], r["pos"], r["plant"]))
    return rows


def snp_gene_summary(count_records: List[Dict]) -> List[Dict]:
    """One row per (gene, SNP): flowers and plants collapsed.

    The fundamental ASE evidence unit — the raw (observed) variable/fixed counts
    behind each gene, with a compact per-plant breakdown so the plant-level split
    is still visible without the full gene×SNP×plant expansion.
    """
    by_gene: Dict[str, List[Dict]] = defaultdict(list)
    for r in count_records:
        if r.get("gene_id") in (None, "", "intergenic"):
            continue
        by_gene[r["gene_id"]].append(r)

    rows: List[Dict] = []
    for gene_id, recs in by_gene.items():
        null_p = _gene_null(recs)
        by_snp: Dict[str, List[Dict]] = defaultdict(list)
        for r in recs:
            by_snp[r["snp_id"]].append(r)
        for snp_id, srecs in by_snp.items():
            m = srecs[0]
            v = sum(r["variable_count"] for r in srecs)
            f = sum(r["fixed_count"] for r in srecs)
            o = sum(r["other_count"] for r in srecs)
            depth = v + f
            per_plant: Dict[str, List[int]] = defaultdict(lambda: [0, 0])
            for r in srecs:
                per_plant[r["plant"]][0] += r["variable_count"]
                per_plant[r["plant"]][1] += r["variable_count"] + r["fixed_count"]
            pp = ";".join(f"{p}={vv}/{nn}"
                          for p, (vv, nn) in sorted(per_plant.items()))
            rows.append({
                "gene_id": gene_id, "gene_name": m.get("gene_name", gene_id),
                "chrom": m["chrom"], "pos": m["pos"], "snp_id": snp_id,
                "tier": m.get("tier", ""),
                "n_plants": len(per_plant),
                "n_flowers": len({r["flower"] for r in srecs}),
                "variable_count": v, "fixed_count": f, "other_count": o,
                "allelic_depth": depth,
                "variable_ratio": round(v / depth, 4) if depth else "NA",
                "per_plant_counts": pp,
                "null_p": round(null_p, 4),
            })
    rows.sort(key=lambda r: (r["gene_id"], r["chrom"], r["pos"]))
    return rows


def plant_gene_detail(count_records: List[Dict], *,
                      flower_norm: str = "equalize") -> List[Dict]:
    """One row per (gene, plant): the exact input and output of each per-plant
    test (K_P, N_P, n_snps, rho, method, p) that the max-p combine consumes.

    Uses the same flower size-factor normalisation as `test_genes`, so these
    K/N are precisely the numbers the per-plant test was run on."""
    by_gene: Dict[str, List[Dict]] = defaultdict(list)
    for r in count_records:
        if r.get("gene_id") in (None, "", "intergenic"):
            continue
        by_gene[r["gene_id"]].append(r)

    bg_of_plant = {r["plant"]: r["background"] for r in count_records}
    name_of_gene = {r["gene_id"]: r.get("gene_name", r["gene_id"]) for r in count_records}
    size_factors = flower_size_factors(count_records, mode=flower_norm)

    rows: List[Dict] = []
    for gene_id, recs in by_gene.items():
        null_p = _gene_null(recs)
        for plant, triples in _plant_snp_counts(recs, size_factors,
                                                with_null=True).items():
            res = _plant_test(triples)
            if res is None:
                continue
            k, n = res["k"], res["n"]
            rows.append({
                "gene_id": gene_id, "gene_name": name_of_gene.get(gene_id, gene_id),
                "plant": plant, "background": bg_of_plant.get(plant, ""),
                "n_snps": res["n_snps"],
                "variable_reads": k, "fixed_reads": n - k, "allelic_depth": n,
                "variable_ratio": round(k / n, 4) if n else "NA",
                "log2_ratio": round(res["log2"], 4),
                "rho": round(res["rho"], 4) if res["rho"] is not None else "NA",
                "method": res["method"], "p_plant": res["p"],
                "null_p": round(null_p, 4),
            })
    rows.sort(key=lambda r: (r["gene_id"], r["plant"]))
    return rows


def _phase_concordant(records: List[Dict]) -> bool:
    """Within each plant, do the gene's per-SNP imbalance directions agree?

    A `phased`-tier SNP whose parent was mis-called homozygous (parental ASE)
    can be summed in the wrong orientation; if a plant's SNPs disagree on
    direction the gene ratio is suspect. Returns True if every plant's SNPs are
    directionally concordant (ignoring near-balanced SNPs), else False.

    Each SNP is judged against its OWN null (which differs across SNPs under a
    per-SNP null-shift control), not a single gene-level null — otherwise SNPs
    that are all consistently shifted relative to their individual nulls would
    look discordant against an averaged null (audit C1 interaction).
    """
    per_plant_snp: Dict[str, Dict[str, List[float]]] = defaultdict(
        lambda: defaultdict(lambda: [0, 0, 0.5]))
    for r in records:
        slot = per_plant_snp[r["plant"]][r["snp_id"]]
        slot[0] += r["variable_count"]
        slot[1] += r["variable_count"] + r["fixed_count"]
        slot[2] = r.get("null_p", 0.5)           # per-SNP null
    for snps in per_plant_snp.values():
        dirs = set()
        for v, n, snp_null in snps.values():
            if n < 1:
                continue
            d = st.direction(cont_ratio(int(v), int(n)), snp_null, tol=0.1)
            if d != "balanced":
                dirs.add(d)
        if len(dirs) > 1:
            return False
    return True
