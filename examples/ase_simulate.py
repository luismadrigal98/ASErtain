"""Fragment-level ASE simulator for benchmarking gene-aggregation strategies.

The whole point of the simulator is to reproduce, from a *single ground truth*,
the two count tables the pipeline can build:

  * the per-SNP **pileup** table  (`--counter pileup`):  a read that spans k
    informative SNPs is counted at EACH of them, so within-gene SNPs are in LD
    and the per-gene depth is inflated ~k-fold;
  * the read-backed **haplotype** table (`--counter haplotype`): the same read
    is assigned to one parental haplotype and counted ONCE per gene.

Because both tables come from the same simulated fragments, any difference in
the genes they call is attributable to the aggregation strategy, not to
different data. Each gene has a known truth label (ASE / not, direction,
strength), so we can score sensitivity and false-discovery as well as the
pairwise overlap between approaches.

Records use the exact ASErtain allele-count schema, so `testing.test_genes`
consumes them with no adaptation.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Dict, List, Tuple


# ---------------------------------------------------------------------------
# Design: plants nested in backgrounds, flowers nested in plants
# ---------------------------------------------------------------------------

@dataclass
class Design:
    # background -> list of plant names
    backgrounds: Dict[str, List[str]] = field(default_factory=lambda: {
        "bg1": ["P1", "P2"], "bg2": ["P3", "P4"]})
    flowers_per_plant: int = 3
    # Per-flower library size factor (technical-replicate depth differences):
    # most flowers similar, the occasional deep one.
    flower_lib_lognorm_sigma: float = 0.35

    def plant_background(self) -> Dict[str, str]:
        return {pl: bg for bg, pls in self.backgrounds.items() for pl in pls}

    def plants(self) -> List[str]:
        return [pl for pls in self.backgrounds.values() for pl in pls]


@dataclass
class GeneTruth:
    gene_id: str
    is_ase: bool
    direction: str          # 'variable' | 'fixed' | 'none'
    true_frac: float        # variable-allele fraction (0.5 == balanced)
    n_snps: int
    strength: str           # 'strong' | 'moderate' | 'weak' | 'null'


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

class ASESimulator:
    def __init__(self, design: Design = None, *, seed: int = 0,
                 base_error: float = 0.005,
                 mean_frags_per_gene_flower: float = 45.0,
                 frag_span_min: int = 1, frag_span_max: int = 4):
        self.design = design or Design()
        self.rng = random.Random(seed)
        self.base_error = base_error
        self.mean_frags = mean_frags_per_gene_flower
        self.frag_span_min = frag_span_min
        self.frag_span_max = frag_span_max

    # -- gene truth -------------------------------------------------------
    def make_genes(self, n_genes: int, frac_ase: float) -> List[GeneTruth]:
        genes: List[GeneTruth] = []
        for i in range(n_genes):
            gid = f"g{i:04d}"
            n_snps = self.rng.randint(2, 8)
            if self.rng.random() < frac_ase:
                # Strength tiers -> variable-allele fraction away from 0.5.
                strength = self.rng.choices(
                    ["strong", "moderate", "weak"], weights=[0.45, 0.4, 0.15])[0]
                mag = {"strong": self.rng.uniform(0.78, 0.92),
                       "moderate": self.rng.uniform(0.66, 0.75),
                       "weak": self.rng.uniform(0.58, 0.63)}[strength]
                direction = self.rng.choice(["variable", "fixed"])
                frac = mag if direction == "variable" else 1.0 - mag
                genes.append(GeneTruth(gid, True, direction, frac, n_snps, strength))
            else:
                genes.append(GeneTruth(gid, False, "none", 0.5, n_snps, "null"))
        return genes

    # -- per-flower library factors --------------------------------------
    def _flower_libfactors(self) -> Dict[Tuple[str, str], float]:
        f: Dict[Tuple[str, str], float] = {}
        for pl in self.design.plants():
            for j in range(self.design.flowers_per_plant):
                fl = f"{pl}_fl{j+1}"
                f[(pl, fl)] = self.rng.lognormvariate(
                    0.0, self.design.flower_lib_lognorm_sigma)
        return f

    # -- main: simulate both count tables --------------------------------
    def simulate(self, genes: List[GeneTruth]
                 ) -> Tuple[List[Dict], List[Dict]]:
        """Return (pileup_records, haplotype_records) for all genes/flowers."""
        bg_of = self.design.plant_background()
        libf = self._flower_libfactors()
        pileup: List[Dict] = []
        hap: List[Dict] = []

        for g in genes:
            # Lay SNPs out along the gene; per-plant true fraction (jittered).
            positions = [1000 + 60 * k for k in range(g.n_snps)]
            for pl in self.design.plants():
                bg = bg_of[pl]
                if g.is_ase:
                    pfrac = min(0.97, max(0.03,
                               g.true_frac + self.rng.gauss(0, 0.02)))
                else:
                    pfrac = 0.5
                for j in range(self.design.flowers_per_plant):
                    fl = f"{pl}_fl{j+1}"
                    n_frag = max(1, int(self.rng.gauss(
                        self.mean_frags * libf[(pl, fl)],
                        0.15 * self.mean_frags)))
                    self._simulate_flower(
                        g, positions, pl, bg, fl, pfrac, n_frag, pileup, hap)
        return pileup, hap

    def _simulate_flower(self, g, positions, plant, bg, flower, pfrac,
                         n_frag, pileup, hap):
        n = g.n_snps
        # per-SNP tallies
        snp_var = [0] * n
        snp_fix = [0] * n
        snp_oth = [0] * n
        hap_var = hap_fix = hap_amb = 0

        for _ in range(n_frag):
            # Fragment's true haplotype.
            true_var = self.rng.random() < pfrac
            # Which consecutive SNPs this fragment covers.
            span = self.rng.randint(self.frag_span_min,
                                    min(self.frag_span_max, n))
            start = self.rng.randint(0, n - span)
            covered = range(start, start + span)
            votes_var = votes_fix = 0
            for s in covered:
                # Allele actually observed (base error flips it).
                obs_var = true_var
                if self.rng.random() < self.base_error:
                    obs_var = not obs_var
                if obs_var:
                    snp_var[s] += 1
                    votes_var += 1
                else:
                    snp_fix[s] += 1
                    votes_fix += 1
            # Haplotype assignment for the whole fragment (counted once).
            if votes_var > 0 and votes_fix == 0:
                hap_var += 1
            elif votes_fix > 0 and votes_var == 0:
                hap_fix += 1
            else:
                hap_amb += 1

        # Emit per-SNP pileup records (var allele 'A', fixed 'T').
        for s in range(n):
            depth = snp_var[s] + snp_fix[s] + snp_oth[s]
            if depth == 0:
                continue
            pileup.append({
                "flower": flower, "plant": plant, "background": bg,
                "chrom": "chr1", "pos": positions[s],
                "snp_id": f"{g.gene_id}:{positions[s]}",
                "variable_allele": "A", "fixed_allele": "T",
                "variable_is_ref": "True", "tier": "both_hom",
                "variable_count": snp_var[s], "fixed_count": snp_fix[s],
                "other_count": snp_oth[s], "total_depth": depth,
                "null_p": 0.5, "gene_id": g.gene_id, "gene_name": g.gene_id,
            })
        # Emit one haplotype record for the gene/flower.
        if hap_var + hap_fix > 0:
            hap.append({
                "flower": flower, "plant": plant, "background": bg,
                "chrom": "chr1", "pos": positions[0],
                "snp_id": f"hap:{g.gene_id}",
                "variable_allele": "hap", "fixed_allele": "hap",
                "variable_is_ref": "True", "tier": "haplotype",
                "variable_count": hap_var, "fixed_count": hap_fix,
                "other_count": hap_amb, "total_depth": hap_var + hap_fix + hap_amb,
                "null_p": 0.5, "gene_id": g.gene_id, "gene_name": g.gene_id,
                "n_hap_snps": g.n_snps,
            })
