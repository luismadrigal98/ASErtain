"""Replicate-aware ASE statistics (pure Python / scipy).

Design rationale (the part the IGV-by-eye approach can't give you):

* The unit of biological replication is the **F1 individual**, not the SNP and
  not the read. Linked SNPs in one gene are not independent, and pooling reads
  into a single binomial ignores individual-to-individual overdispersion
  (anti-conservative).
* Primary gene-level test: collapse each gene's diagnostic SNPs to one
  (variable, fixed) count *per replicate*, take the logit allelic ratio per
  replicate, and run a one-sample t-test of those logits against the null
  (0.5, or a bias-adjusted expectation). With a handful of F1 replicates this
  has real power and the right error model.
* Secondary test: a beta-binomial fit across replicates that uses read depth
  while still modelling overdispersion (LRT of mean vs null).
* We also report the naive pooled binomial purely as a descriptive number,
  clearly flagged as anti-conservative.
* Consistency across variable-species backgrounds is reported explicitly: a gene
  biased the same direction in the F1s of every parental background is the robust
  call.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import optimize, stats


# ---------------------------------------------------------------------------
# Basic helpers
# ---------------------------------------------------------------------------

def binomial_p(k: int, n: int, p: float = 0.5) -> float:
    """Two-sided exact binomial test p-value."""
    if n == 0:
        return 1.0
    return float(stats.binomtest(k, n, p, alternative="two-sided").pvalue)


def bh_adjust(pvals: Sequence[float]) -> List[float]:
    """Benjamini–Hochberg FDR-adjusted p-values, preserving input order."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), 1):
        q = pvals[idx] * n / (n - rank + 1)
        prev = min(prev, q)
        adj[idx] = prev
    return adj


def logit(p: float, eps: float = 1e-6) -> float:
    p = min(max(p, eps), 1 - eps)
    return math.log(p / (1 - p))


def cont_ratio(k: int, n: int) -> float:
    """Haldane-corrected proportion, safe at 0 and n."""
    return (k + 0.5) / (n + 1)


# ---------------------------------------------------------------------------
# Per-replicate logit t-test (primary)
# ---------------------------------------------------------------------------

@dataclass
class ReplicateTest:
    n_reps: int
    mean_logit: float
    null_logit: float
    t_stat: float
    p_value: float
    mean_ratio: float            # back-transformed mean replicate ratio
    log2_ratio: float            # log2(variable/fixed) at the mean


def replicate_logit_test(rep_counts: Sequence[Tuple[int, int]],
                         null_p: float = 0.5) -> Optional[ReplicateTest]:
    """One-sample t-test of per-replicate logit ratios vs the null.

    rep_counts: list of (variable_reads, total_reads) — one tuple per replicate.
    Requires >= 2 replicates with non-zero depth.
    """
    pairs = [(k, n) for k, n in rep_counts if n > 0]
    if len(pairs) < 2:
        return None
    logits = np.array([logit(cont_ratio(k, n)) for k, n in pairs])
    null = logit(null_p)
    res = stats.ttest_1samp(logits, popmean=null)
    mean_logit = float(logits.mean())
    mean_ratio = 1 / (1 + math.exp(-mean_logit))
    return ReplicateTest(
        n_reps=len(pairs),
        mean_logit=mean_logit,
        null_logit=null,
        t_stat=float(res.statistic),
        p_value=float(res.pvalue),
        mean_ratio=mean_ratio,
        log2_ratio=math.log2(mean_ratio / (1 - mean_ratio)),
    )


# ---------------------------------------------------------------------------
# Beta-binomial (secondary, depth-aware, overdispersed)
# ---------------------------------------------------------------------------

def _bb_params(mu: float, rho: float) -> Tuple[float, float]:
    """Convert (mean, intra-class correlation) to scipy betabinom (a, b)."""
    rho = min(max(rho, 1e-6), 1 - 1e-6)
    s = (1 - rho) / rho            # alpha + beta
    return mu * s, (1 - mu) * s


def _bb_negloglik(theta, ks, ns) -> float:
    mu = 1 / (1 + math.exp(-theta[0]))
    rho = 1 / (1 + math.exp(-theta[1]))
    a, b = _bb_params(mu, rho)
    ll = stats.betabinom.logpmf(ks, ns, a, b)
    if not np.all(np.isfinite(ll)):
        return 1e12
    return -float(np.sum(ll))


@dataclass
class BetaBinomTest:
    mu: float                 # fitted mean variable-allele fraction
    rho: float                # fitted overdispersion (intra-class correlation)
    log2_ratio: float
    p_value: float            # LRT vs mu == null_p
    converged: bool


def betabinom_test(rep_counts: Sequence[Tuple[int, int]],
                   null_p: float = 0.5) -> Optional[BetaBinomTest]:
    """Beta-binomial LRT of the mean allelic fraction against null_p.

    Uses read depth across replicates and an estimated overdispersion. Returns
    None if fewer than 2 informative replicates.
    """
    pairs = [(k, n) for k, n in rep_counts if n > 0]
    if len(pairs) < 2:
        return None
    ks = np.array([k for k, _ in pairs], dtype=float)
    ns = np.array([n for _, n in pairs], dtype=float)

    pooled = float(ks.sum() / ns.sum())
    x0 = [logit(pooled), 0.0]      # start: pooled mean, rho=0.5

    full = optimize.minimize(_bb_negloglik, x0, args=(ks, ns),
                             method="Nelder-Mead",
                             options={"xatol": 1e-4, "fatol": 1e-4,
                                      "maxiter": 2000})
    mu = 1 / (1 + math.exp(-full.x[0]))
    rho = 1 / (1 + math.exp(-full.x[1]))

    # Null model: mu fixed at null_p, optimise rho only.
    def _null_nll(theta_rho):
        a, b = _bb_params(null_p, 1 / (1 + math.exp(-theta_rho[0])))
        ll = stats.betabinom.logpmf(ks, ns, a, b)
        if not np.all(np.isfinite(ll)):
            return 1e12
        return -float(np.sum(ll))

    null = optimize.minimize(_null_nll, [full.x[1]], method="Nelder-Mead",
                            options={"xatol": 1e-4, "fatol": 1e-4,
                                     "maxiter": 2000})
    lr = 2 * (null.fun - full.fun)
    p = float(stats.chi2.sf(max(lr, 0.0), df=1))
    return BetaBinomTest(
        mu=mu, rho=rho,
        log2_ratio=math.log2(mu / (1 - mu)) if 0 < mu < 1 else float("nan"),
        p_value=p,
        converged=bool(full.success and null.success),
    )


# ---------------------------------------------------------------------------
# Consistency across backgrounds
# ---------------------------------------------------------------------------

def direction(log2_ratio: float, tol: float = 0.0) -> str:
    if log2_ratio > tol:
        return "variable"
    if log2_ratio < -tol:
        return "fixed"
    return "balanced"


def consistent_across(background_log2: Dict[str, float]) -> bool:
    """True if every background with a defined effect points the same way."""
    dirs = {direction(v) for v in background_log2.values()
            if v is not None and not math.isnan(v)}
    dirs.discard("balanced")
    return len(dirs) == 1
