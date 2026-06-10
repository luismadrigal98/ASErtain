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


def stouffer_combine(items: Sequence[Tuple[float, str]]) -> float:
    """Directional Stouffer combination of per-plant (two-sided p, direction).

    `items` is one (p_value, direction) per F1 plant, direction in
    {'variable','fixed','balanced'}. Each plant is converted to a signed z in its
    own direction (variable = +, fixed = -, balanced = 0) and combined with EQUAL
    weight (the plant is the unit of inference), giving a two-sided combined p.

    This is the scalable alternative to the max-p intersection–union rule: max-p
    requires EVERY plant to be individually significant (honest and conservative
    at n=2, but its power falls as plants are added), whereas Stouffer tests for a
    consistent aggregate shift and gains power with more replicates. ASErtain
    still gates an ASE call on the cross-background consistency requirement, so a
    Stouffer call additionally needs the plants to agree in direction. Use max-p
    for the small-n consistency claim; Stouffer when scaling to several plants.
    """
    pairs = [(min(max(p, 1e-300), 1.0), d) for p, d in items]
    if not pairs:
        return 1.0
    sign = {"variable": 1.0, "fixed": -1.0, "balanced": 0.0}
    zs = [sign.get(d, 0.0) * float(stats.norm.isf(p / 2.0)) for p, d in pairs]
    z = sum(zs) / math.sqrt(len(zs))
    return float(2.0 * stats.norm.sf(abs(z)))


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
    """Convert the interpretable (mean, overdispersion) to scipy's Beta (a, b).

    The same Beta-binomial has three equivalent parametrisations; we interpret
    and optimise in one and evaluate the likelihood in another:

      (a, b)   native Beta SHAPE params that scipy.stats.betabinom needs.
               Read as pseudo-counts: a = prior 'successes' (variable allele),
               b = prior 'failures' (fixed allele); mean = a/(a+b).
               Drawback: a and b are entangled — changing a moves the mean AND
               the spread, which is bad for a test about the mean alone.

      (mu, s)  mu = a/(a+b)  is the MEAN (the biology: variable-allele fraction).
               s  = a + b    is the CONCENTRATION (total pseudo-count): large s
               => tight spike at mu (little overdispersion; s->inf is the plain
               binomial); small s => broad (much overdispersion). Inverting:
                   a = mu * s,   b = (1 - mu) * s        <-- this function.

      (mu, rho) rho = 1/(s+1) in (0,1) is the bounded OVERDISPERSION /
               intra-class correlation (correlation of two reads from one SNP);
               Var(p_s) = mu(1-mu)*rho. Hence  s = (1 - rho)/rho.

    Why three: scipy needs (a, b); the hypothesis is about mu only, so we hold
    mu fixed under H0 and let rho float (trivial in (mu, rho), awkward in (a, b))
    — this orthogonalises the parameter of interest from the nuisance. The
    optimiser works in logit(mu), logit(rho) (unconstrained reals), then maps
    logit -> (mu, rho) -> s=(1-rho)/rho -> (a, b) here -> betabinom.logpmf.
    """
    rho = min(max(rho, 1e-6), 1 - 1e-6)
    s = (1 - rho) / rho                # concentration a+b from overdispersion
    return mu * s, (1 - mu) * s        # a = mu*s (successes), b = (1-mu)*s (failures)


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


# ---------------------------------------------------------------------------
# Per-SNP-null shift test (for heterogeneous nulls, e.g. --bias-mode null-shift)
# ---------------------------------------------------------------------------
#
# When a gene's SNPs do NOT share a null (each SNP has its own expected
# variable-allele fraction p0_i, as under a per-SNP balanced-control table),
# you must NOT pool the counts and test the pooled ratio against an averaged
# null — that conflates real imbalance with heterogeneity of the nulls
# themselves (audit C1). Instead we fit a single *logit shift* delta shared
# across the plant's SNPs:
#
#     logit(p_i) = logit(p0_i) + delta            (H1: delta free)
#     logit(p_i) = logit(p0_i)                    (H0: delta = 0; each SNP at
#                                                  its own null)
#
# and test delta = 0 by a likelihood-ratio test (df = 1). delta is the common
# allelic shift on the log-odds scale, so a single coherent ASE direction is
# estimated while each SNP is correctly centred on its own bias. Overdispersion
# (rho) is carried when there are enough SNPs, exactly as in `betabinom_test`;
# otherwise the binomial (rho -> 0) form is used.


@dataclass
class ShiftTest:
    delta: float          # logit-scale shift shared across the SNPs
    rho: Optional[float]  # overdispersion (None for the binomial form)
    log2_ratio: float     # delta expressed as a log2-odds shift (sign = direction)
    p_value: float
    method: str           # 'betabinom_shift' | 'binom_shift'
    converged: bool


def _logit_np(p):
    p = np.clip(np.asarray(p, dtype=float), 1e-9, 1 - 1e-9)
    return np.log(p / (1 - p))


def _shift_binom_nll(delta: float, ks, ns, base_logit) -> float:
    p = 1.0 / (1.0 + np.exp(-(base_logit + delta)))
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -float(np.sum(ks * np.log(p) + (ns - ks) * np.log(1 - p)))


def _shift_bb_nll(delta: float, rho: float, ks, ns, base_logit) -> float:
    rho = min(max(rho, 1e-6), 1 - 1e-6)
    s = (1 - rho) / rho
    mu = 1.0 / (1.0 + np.exp(-(base_logit + delta)))
    a, b = mu * s, (1 - mu) * s
    ll = stats.betabinom.logpmf(ks, ns, a, b)
    if not np.all(np.isfinite(ll)):
        return 1e12
    return -float(np.sum(ll))


def shift_test(obs: Sequence[Tuple[int, int, float]]) -> Optional[ShiftTest]:
    """LRT for a common logit shift away from each observation's own null.

    `obs` is a list of (k, n, null_p) — one per SNP within a plant, where the
    nulls may differ. Needs >= 1 observation with n > 0. Uses the beta-binomial
    (overdispersed) form when there are >= MIN_SNPS_FOR_BETABINOM SNPs and it
    converges, else the binomial form. Returns None if there is nothing to test.
    """
    triples = [(k, n, p0) for k, n, p0 in obs if n > 0]
    if not triples:
        return None
    ks = np.array([k for k, _, _ in triples], dtype=float)
    ns = np.array([n for _, n, _ in triples], dtype=float)
    base = _logit_np([p0 for _, _, p0 in triples])

    use_bb = len(triples) >= 3      # mirrors MIN_SNPS_FOR_BETABINOM in testing
    if use_bb:
        best = None
        for r in _RHO_SEEDS:
            res = optimize.minimize(
                lambda t: _shift_bb_nll(t[0], 1 / (1 + math.exp(-t[1])),
                                        ks, ns, base),
                [0.0, logit(r)], method="Nelder-Mead",
                options={"xatol": 1e-5, "fatol": 1e-5, "maxiter": 4000})
            if best is None or res.fun < best.fun:
                best = res
        # Refit rho under H0 (delta fixed at 0) for a proper nested LRT.
        null_fit = None
        for r in _RHO_SEEDS:
            res0 = optimize.minimize(
                lambda t: _shift_bb_nll(0.0, 1 / (1 + math.exp(-t[0])),
                                        ks, ns, base),
                [logit(r)], method="Nelder-Mead",
                options={"xatol": 1e-5, "fatol": 1e-5, "maxiter": 4000})
            if null_fit is None or res0.fun < null_fit.fun:
                null_fit = res0
        if best.success and null_fit.success:
            delta = float(best.x[0])
            rho = 1 / (1 + math.exp(-best.x[1]))
            lr = 2 * (null_fit.fun - best.fun)
            p = float(stats.chi2.sf(max(lr, 0.0), df=1))
            return ShiftTest(delta=delta, rho=rho,
                             log2_ratio=delta / math.log(2), p_value=p,
                             method="betabinom_shift",
                             converged=bool(lr >= -1e-6))
        # fall through to the binomial form if the bb fit did not converge

    full = optimize.minimize(
        lambda t: _shift_binom_nll(t[0], ks, ns, base),
        [0.0], method="Nelder-Mead",
        options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 4000})
    delta = float(full.x[0])
    lr = 2 * (_shift_binom_nll(0.0, ks, ns, base) - full.fun)
    p = float(stats.chi2.sf(max(lr, 0.0), df=1))
    return ShiftTest(delta=delta, rho=None, log2_ratio=delta / math.log(2),
                     p_value=p, method="binom_shift",
                     converged=bool(full.success and lr >= -1e-6))
