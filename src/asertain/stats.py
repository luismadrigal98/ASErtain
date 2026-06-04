"""Statistics for ASE, sized honestly for few biological replicates.

Design rationale (see DESIGN.md and AUDIT.md):

* The biological unit is the F1 individual (plant); flowers are technical
  sub-samples summed into their plant.
* With only a handful of plants you cannot calibrate a beta-binomial *across
  plants* (its overdispersion is unidentifiable at n=2 — audit finding C2).
  Instead we test **each plant on its own** — a beta-binomial across that
  plant's informative SNPs when there are enough SNPs, else a pooled binomial —
  and combine plants with an **intersection–union (max-p)** rule: a gene is ASE
  only if *every* contributing plant is individually significant in the same
  direction. Max-p is a valid (conservative) test of "all plants imbalanced",
  needs no across-plant variance estimate, and makes the cross-background
  consistency a requirement rather than a footnote.
* The per-plant beta-binomial uses a multi-start optimiser and reports a
  convergence flag; non-converged fits fall back to the binomial.
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
    p = min(max(p, 1e-6), 1 - 1e-6)
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


def direction(ratio: float, null_p: float = 0.5, tol: float = 0.0) -> str:
    """Direction of imbalance of an observed ratio relative to its null."""
    if ratio > null_p + tol:
        return "variable"
    if ratio < null_p - tol:
        return "fixed"
    return "balanced"


def consistent_across(background_log2: Dict[str, float]) -> bool:
    """True if every background with a defined effect points the same way AND at
    least two backgrounds contribute (a single background cannot be 'consistent';
    audit finding M4)."""
    dirs = {("variable" if v > 0 else "fixed")
            for v in background_log2.values()
            if v is not None and not (isinstance(v, float) and math.isnan(v))
            and v != 0}
    return len(background_log2) >= 2 and len(dirs) == 1


# ---------------------------------------------------------------------------
# Per-replicate logit t-test (secondary / descriptive)
# ---------------------------------------------------------------------------

@dataclass
class ReplicateTest:
    n_reps: int
    mean_logit: float
    t_stat: float
    p_value: float
    mean_ratio: float
    log2_ratio: float


def replicate_logit_test(rep_counts: Sequence[Tuple[int, int]],
                         null_p: float = 0.5) -> Optional[ReplicateTest]:
    """One-sample t-test of per-replicate logit ratios vs the null.

    Secondary/descriptive only (the primary path is the per-plant max-p rule).
    A variance floor prevents the t-statistic from diverging to +/-inf (and the
    p-value from collapsing to exactly 0) when replicates agree — audit M1.
    """
    pairs = [(k, n) for k, n in rep_counts if n > 0]
    if len(pairs) < 2:
        return None
    logits = np.array([logit(cont_ratio(k, n)) for k, n in pairs])
    null = logit(null_p)
    mean_logit = float(logits.mean())
    sd = float(logits.std(ddof=1))
    sd = max(sd, 0.05)                      # floor: never claim zero variance
    se = sd / math.sqrt(len(pairs))
    t = (mean_logit - null) / se if se > 0 else 0.0
    p = float(2 * stats.t.sf(abs(t), df=len(pairs) - 1))
    mean_ratio = 1 / (1 + math.exp(-mean_logit))
    return ReplicateTest(
        n_reps=len(pairs), mean_logit=mean_logit, t_stat=t, p_value=p,
        mean_ratio=mean_ratio,
        log2_ratio=math.log2(mean_ratio / (1 - mean_ratio)))


# ---------------------------------------------------------------------------
# Beta-binomial over a set of (k, n) observations (SNPs within one plant)
# ---------------------------------------------------------------------------

def _bb_params(mu: float, rho: float) -> Tuple[float, float]:
    rho = min(max(rho, 1e-6), 1 - 1e-6)
    s = (1 - rho) / rho
    return mu * s, (1 - mu) * s


def _bb_nll(mu: float, rho: float, ks, ns) -> float:
    a, b = _bb_params(mu, rho)
    ll = stats.betabinom.logpmf(ks, ns, a, b)
    if not np.all(np.isfinite(ll)):
        return 1e12
    return -float(np.sum(ll))


_RHO_SEEDS = (0.005, 0.05, 0.2, 0.5, 0.85)


def _fit_full(ks, ns, mu0: float):
    best = None
    for r in _RHO_SEEDS:
        res = optimize.minimize(
            lambda t: _bb_nll(1 / (1 + math.exp(-t[0])),
                              1 / (1 + math.exp(-t[1])), ks, ns),
            [logit(mu0), logit(r)], method="Nelder-Mead",
            options={"xatol": 1e-5, "fatol": 1e-5, "maxiter": 4000})
        if best is None or res.fun < best.fun:
            best = res
    return best


def _fit_null(ks, ns, mu_fixed: float):
    best = None
    for r in _RHO_SEEDS:
        res = optimize.minimize(
            lambda t: _bb_nll(mu_fixed, 1 / (1 + math.exp(-t[0])), ks, ns),
            [logit(r)], method="Nelder-Mead",
            options={"xatol": 1e-5, "fatol": 1e-5, "maxiter": 4000})
        if best is None or res.fun < best.fun:
            best = res
    return best


@dataclass
class BetaBinomTest:
    mu: float
    rho: float
    log2_ratio: float
    p_value: float
    converged: bool


def betabinom_test(obs: Sequence[Tuple[int, int]],
                   null_p: float = 0.5) -> Optional[BetaBinomTest]:
    """Beta-binomial LRT that the mean fraction differs from null_p.

    `obs` is a list of (k, n) — typically one per SNP within a single plant.
    Needs >= 2 informative observations. Uses a multi-start optimiser for both
    the full and (mu-fixed) null models, so the overdispersion that protects
    extreme data is actually found (audit C2). Returns None if < 2 obs.
    """
    pairs = [(k, n) for k, n in obs if n > 0]
    if len(pairs) < 2:
        return None
    ks = np.array([k for k, _ in pairs], dtype=float)
    ns = np.array([n for _, n in pairs], dtype=float)
    pooled = float(min(max(ks.sum() / ns.sum(), 1e-6), 1 - 1e-6))

    full = _fit_full(ks, ns, pooled)
    null = _fit_null(ks, ns, null_p)
    mu = 1 / (1 + math.exp(-full.x[0]))
    rho = 1 / (1 + math.exp(-full.x[1]))
    lr = 2 * (null.fun - full.fun)
    p = float(stats.chi2.sf(max(lr, 0.0), df=1))
    return BetaBinomTest(
        mu=mu, rho=rho,
        log2_ratio=math.log2(mu / (1 - mu)) if 0 < mu < 1 else float("nan"),
        p_value=p,
        converged=bool(full.success and null.success and lr >= -1e-6))
