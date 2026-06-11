"""Benchmark gene-aggregation strategies on simulated F1 ASE data.

Runs three ways of turning a gene's informative SNPs into one ASE call on the
SAME simulated fragments:

  P  pileup  + plant-aggregation  (per-SNP beta-binomial; within-gene SNPs in
             LD are pseudo-replicated -> inflated depth)         [LD-naive]
  H  haplotype + plant-aggregation (read-backed: one independent (K,N)/gene)
             -- the method the lab wants to keep                 [LD-correct]
  M  pileup  + maxsnp-aggregation (plain binomial per SNP, strongest SNP per
             gene, require SNPs to agree in direction)           [advisor's]

and reports, against the known truth, each method's sensitivity / empirical FDR
and the pairwise overlap of the gene sets they call.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))

from asertain import testing
from ase_simulate import ASESimulator, Design


def called_set(genes_result):
    return {g["gene_id"] for g in genes_result if g["ase_call"]}


def direction_of(genes_result):
    return {g["gene_id"]: g["direction"] for g in genes_result}


def score(name, called, truth_ase, truth_dir, pred_dir):
    tp = called & truth_ase
    fp = called - truth_ase
    sens = len(tp) / len(truth_ase) if truth_ase else float("nan")
    fdr = len(fp) / len(called) if called else 0.0
    # direction accuracy among true positives
    dir_ok = sum(1 for g in tp if pred_dir.get(g) == truth_dir[g])
    dir_acc = dir_ok / len(tp) if tp else float("nan")
    print(f"  {name:36s}  called={len(called):4d}  "
          f"sensitivity={sens:5.1%}  emp.FDR={fdr:5.1%}  "
          f"dir.acc={dir_acc:5.1%}")
    return dict(called=called, sens=sens, fdr=fdr)


def overlap(a, b):
    inter = len(a & b)
    union = len(a | b)
    jacc = inter / union if union else 1.0
    omin = inter / min(len(a), len(b)) if min(len(a), len(b)) else 1.0
    return jacc, omin, inter, union


def run(n_genes=400, frac_ase=0.4, seed=1, alpha=0.05,
        min_effect_log2=0.5, correction="sidak", verbose=True):
    sim = ASESimulator(Design(), seed=seed)
    genes = sim.make_genes(n_genes, frac_ase)
    pileup, hap = sim.simulate(genes)

    truth_ase = {g.gene_id for g in genes if g.is_ase}
    truth_dir = {g.gene_id: g.direction for g in genes}

    P = testing.test_genes(pileup, alpha=alpha, min_effect_log2=min_effect_log2,
                           gene_aggregation="plant")
    H = testing.test_genes(hap, alpha=alpha, min_effect_log2=min_effect_log2,
                           gene_aggregation="plant")
    M = testing.test_genes(pileup, alpha=alpha, min_effect_log2=min_effect_log2,
                           gene_aggregation="maxsnp",
                           within_gene_correction=correction)

    cP, cH, cM = called_set(P), called_set(H), called_set(M)
    if verbose:
        print(f"\n=== Simulation: {n_genes} genes, {len(truth_ase)} truly ASE, "
              f"seed={seed}, alpha={alpha}, min|log2|={min_effect_log2}, "
              f"maxsnp-corr={correction} ===")
        print("Per-method performance vs ground truth:")
        score("P  pileup + plant (LD-naive)", cP, truth_ase, truth_dir, direction_of(P))
        score("H  haplotype + plant (keep)", cH, truth_ase, truth_dir, direction_of(H))
        score("M  pileup + maxsnp (advisor)", cM, truth_ase, truth_dir, direction_of(M))

        print("\nPairwise overlap of called gene sets:")
        for name, a, b in [("H vs M (haplotype vs max-SNP)", cH, cM),
                           ("H vs P (haplotype vs pileup)", cH, cP),
                           ("P vs M (pileup vs max-SNP)", cP, cM)]:
            j, om, inter, union = overlap(a, b)
            print(f"  {name:34s}  Jaccard={j:5.1%}  "
                  f"overlap(min)={om:5.1%}  shared={inter}  union={union}")
    return dict(genes=genes, pileup=pileup, hap=hap,
                P=P, H=H, M=M, cP=cP, cH=cH, cM=cM, truth_ase=truth_ase)


if __name__ == "__main__":
    # Aggregate across several seeds for a stable read on the overlap.
    import statistics as stx
    jHM, jHP, jPM = [], [], []
    sH, sM, sP = [], [], []
    fH, fM, fP = [], [], []
    for sd in range(1, 9):
        r = run(seed=sd, verbose=(sd == 1))
        jHM.append(overlap(r["cH"], r["cM"])[0])
        jHP.append(overlap(r["cH"], r["cP"])[0])
        jPM.append(overlap(r["cP"], r["cM"])[0])
        def metr(c):
            tp = len(c & r["truth_ase"]); fp = len(c - r["truth_ase"])
            return (tp / len(r["truth_ase"]), fp / len(c) if c else 0.0)
        s, f = metr(r["cH"]); sH.append(s); fH.append(f)
        s, f = metr(r["cM"]); sM.append(s); fM.append(f)
        s, f = metr(r["cP"]); sP.append(s); fP.append(f)

    print("\n=== Across 8 seeds (mean) ===")
    print(f"  Sensitivity   H={stx.mean(sH):.1%}  M={stx.mean(sM):.1%}  P={stx.mean(sP):.1%}")
    print(f"  Empirical FDR H={stx.mean(fH):.1%}  M={stx.mean(fM):.1%}  P={stx.mean(fP):.1%}")
    print(f"  Jaccard overlap  H~M={stx.mean(jHM):.1%}  H~P={stx.mean(jHP):.1%}  P~M={stx.mean(jPM):.1%}")
