"""Multi-seed benchmark stratified by true effect strength.

H (haplotype+plant) and M (pileup+maxsnp) are fast and run over many seeds for a
stable overlap estimate; the LD-naive baseline P (pileup+plant beta-binomial) is
slower, so it is summarised on a single seed to illustrate the inflation problem.
"""
import os, sys, statistics as stx
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.dirname(__file__))
from asertain import testing
from ase_simulate import ASESimulator, Design

OUT = sys.argv[1] if len(sys.argv) > 1 else "/tmp/bench_report.txt"
SEEDS = list(range(1, 13))
NG = 140
lines = []
def pr(s): lines.append(s); print(s, flush=True)
def called(res): return {g["gene_id"] for g in res if g["ase_call"]}

agg = {"H": {"sens": [], "fdr": []}, "M": {"sens": [], "fdr": []}}
tier_sens = {"H": {}, "M": {}}
jacc_all, omin_all, jacc_sm, omin_sm = [], [], [], []

for sd in SEEDS:
    sim = ASESimulator(Design(), seed=sd)
    genes = sim.make_genes(NG, 0.4)
    pileup, hap = sim.simulate(genes)
    truth = {g.gene_id for g in genes if g.is_ase}
    strong_mod = {g.gene_id for g in genes if g.strength in ("strong", "moderate")}
    tiers = {}
    for g in genes:
        tiers.setdefault(g.strength, set()).add(g.gene_id)
    res = {
        "H": testing.test_genes(hap, min_effect_log2=0.5, gene_aggregation="plant"),
        "M": testing.test_genes(pileup, min_effect_log2=0.5,
                                gene_aggregation="maxsnp", within_gene_correction="none"),
    }
    c = {k: called(v) for k, v in res.items()}
    for k in res:
        cc = c[k]
        agg[k]["sens"].append(len(cc & truth) / len(truth))
        agg[k]["fdr"].append(len(cc - truth) / len(cc) if cc else 0.0)
        for t, gids in tiers.items():
            if t == "null": continue
            tier_sens[k].setdefault(t, []).append(len(cc & gids) / len(gids))
    cH, cM = c["H"], c["M"]
    inter, uni = len(cH & cM), len(cH | cM)
    jacc_all.append(inter / uni if uni else 1.0)
    omin_all.append(inter / min(len(cH), len(cM)) if min(len(cH), len(cM)) else 1.0)
    cHs, cMs = cH & strong_mod, cM & strong_mod
    i2, u2 = len(cHs & cMs), len(cHs | cMs)
    jacc_sm.append(i2 / u2 if u2 else 1.0)
    omin_sm.append(i2 / min(len(cHs), len(cMs)) if min(len(cHs), len(cMs)) else 1.0)

# LD-naive baseline P on one seed.
sim = ASESimulator(Design(), seed=1); genes = sim.make_genes(NG, 0.4)
pileup, hap = sim.simulate(genes)
truth = {g.gene_id for g in genes if g.is_ase}
P = called(testing.test_genes(pileup, min_effect_log2=0.5, gene_aggregation="plant"))
P_sens = len(P & truth) / len(truth); P_fdr = len(P - truth) / len(P) if P else 0.0

pr(f"Multi-seed ASE benchmark  ({len(SEEDS)} seeds x {NG} genes, 40% truly ASE)")
pr("Same simulated fragments feed every method; maxsnp correction=none, min|log2|=0.5, alpha=0.05\n")
pr("Method                          mean sensitivity   mean empirical-FDR")
pr(f"  H  haplotype + plant   (keep)      {stx.mean(agg['H']['sens']):6.1%}            {stx.mean(agg['H']['fdr']):6.1%}")
pr(f"  M  pileup + maxsnp  (advisor)      {stx.mean(agg['M']['sens']):6.1%}            {stx.mean(agg['M']['fdr']):6.1%}")
pr(f"  P  pileup + plant (LD-naive,1 seed){P_sens:6.1%}            {P_fdr:6.1%}")
pr("\nSensitivity by true effect strength (mean over seeds):")
pr(f"  {'tier':10s}  {'H haplotype':>12s}  {'M maxsnp':>10s}")
for t in ("strong", "moderate", "weak"):
    pr(f"  {t:10s}  {stx.mean(tier_sens['H'][t]):12.1%}  {stx.mean(tier_sens['M'][t]):10.1%}")
pr("\nHaplotype (keep) vs max-SNP (advisor) overlap of called gene sets:")
pr(f"  all ASE calls         Jaccard={stx.mean(jacc_all):6.1%}    overlap(min)={stx.mean(omin_all):6.1%}")
pr(f"  strong+moderate only  Jaccard={stx.mean(jacc_sm):6.1%}    overlap(min)={stx.mean(omin_sm):6.1%}")
with open(OUT, "w") as fh: fh.write("\n".join(lines) + "\n")
print("WROTE", OUT)
