"""Parental differential expression — the *total* (cis + trans) divergence.

The F1 allelic ratio measures **cis** only. To classify regulatory divergence
and to *sanity-check* ASE candidates you also need the parental expression
difference (variable lineage vs fixed lineage). This module computes that from
the parental RNA-seq libraries, pure-Python (numpy/scipy + samtools):

    1. count reads per gene per parental library      (`count_parental_expression`)
    2. library-size normalise (DESeq median-of-ratios) (`size_factors`)
    3. per-gene test variable vs fixed                 (`differential_expression`)

The output is oriented variable/fixed (a positive log2 fold change = higher in
the variable lineage) and column-named so it drops straight into
`asertain contrast` (`log2FoldChange`, `padj`).

Honesty, by design. Two limitations are surfaced, never hidden:

* **Gene-region counts, not exon-union.** `samtools view -c` over the gene
  interval is a proxy for a featureCounts/HTSeq exon-union count; it includes
  intronic overlap. Fine for a direction/magnitude sanity check; for a
  publication DE table, run DESeq2/edgeR and pass it to `contrast` directly.
* **Pseudoreplication.** If a lineage is represented by a single genotype
  (e.g. one fixed-lineage parent sampled as several flowers), the flowers are
  technical — not biological — replicates, so the p-value is anticonservative.
  We warn loudly and report `n_genotypes_*`; trust the *direction* more than the
  p-value in that case.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import stats

from . import external, stats as st
from .annotation import GeneIndex
from .config import CrossConfig

DE_COLS = [
    "gene_id", "gene_name", "baseMean",
    "mean_variable", "mean_fixed", "log2FoldChange",
    "pvalue", "padj",
    "n_variable", "n_fixed",                 # biological units (genotypes)
    "n_variable_flowers", "n_fixed_flowers",  # technical sub-samples
    "higher_in", "method",
]
# `higher_in` holds a role token (variable/fixed), so it is relabelled like a
# direction column when written with user labels.
DE_DIRECTION_COLS = ("higher_in",)


# ---------------------------------------------------------------------------
# 1. Per-gene read counting from the parental BAMs
# ---------------------------------------------------------------------------

def count_parental_expression(cfg: CrossConfig, gene_index: GeneIndex, *,
                              gene_ids: Optional[set] = None,
                              min_mapq: int = 20,
                              samtools: str = "samtools",
                              progress: bool = True
                              ) -> Tuple[Dict[str, Dict[str, int]],
                                         Dict[str, str], Dict[str, str],
                                         Dict[str, int], Dict[str, str]]:
    """Count reads per gene per parental library (flower).

    `gene_ids`, if given, restricts counting to those genes (e.g. the ASE
    candidate genes) — important because per-gene `samtools` counting is one
    subprocess per gene×sample, so a transcriptome-wide run is slow. For
    genome-wide DE, use a dedicated counter (featureCounts/HTSeq) + DESeq2 and
    pass that table to `contrast` instead.

    Returns (counts, sample_lineage, gene_names, library_sizes, sample_genotype):
      counts[gene_id][sample] = read count
      sample_lineage[sample]  = 'variable' | 'fixed'
      gene_names[gene_id]     = display name
      library_sizes[sample]   = total mapped reads (for depth normalisation)
      sample_genotype[sample] = parent name (the biological unit the flower nests in)
    """
    import os
    samples: List[Tuple[str, str, str, str]] = []   # (sample, bam, lineage, genotype)
    for p in cfg.parents:
        for fl in p.flowers:
            samples.append((fl.name, fl.bam, p.lineage, p.name))
    if not samples:
        raise ValueError(
            "No parental RNA libraries found. Add `flowers:` (with bam paths) "
            "under the parents in the config to enable the parental-DE stage, "
            "or supply an external DE table to `contrast`.")

    for s, bam, _, _ in samples:
        if not os.path.exists(bam):
            raise FileNotFoundError(f"parental BAM for '{s}' not found: {bam}")
        external.ensure_bam_index(bam, samtools=samtools)

    genes = [(chrom, g) for chrom, g in gene_index.iter_genes()
             if gene_ids is None or g.gene_id in gene_ids]
    if not genes:
        raise ValueError("No genes to count for parental DE (gene_ids filter "
                         "matched nothing in the annotation).")
    gene_names = {g.gene_id: g.gene_name for _, g in genes}
    counts: Dict[str, Dict[str, int]] = {g.gene_id: {} for _, g in genes}
    sample_lineage = {s: lin for s, _, lin, _ in samples}
    sample_genotype = {s: geno for s, _, _, geno in samples}
    library_sizes: Dict[str, int] = {}

    for si, (sample, bam, _lin, _geno) in enumerate(samples, 1):
        if progress:
            print(f"  [{si}/{len(samples)}] counting {sample} over "
                  f"{len(genes)} genes", flush=True)
        library_sizes[sample] = external.samtools_total_mapped(bam, samtools=samtools)
        for chrom, g in genes:
            region = f"{chrom}:{g.start}-{g.end}"
            counts[g.gene_id][sample] = external.samtools_count(
                bam, region, min_mapq=min_mapq, samtools=samtools)
    return counts, sample_lineage, gene_names, library_sizes, sample_genotype


# ---------------------------------------------------------------------------
# 2. Library-size normalisation
# ---------------------------------------------------------------------------

def size_factors(matrix: np.ndarray) -> np.ndarray:
    """DESeq median-of-ratios size factors for a genes×samples count matrix.

    Falls back to total-count factors when too few genes are non-zero in every
    sample (common for small candidate-gene panels)."""
    n_samples = matrix.shape[1]
    if matrix.size == 0:
        return np.ones(n_samples)
    with np.errstate(divide="ignore"):
        logmat = np.log(matrix)
    loggeo = logmat.mean(axis=1)               # per-gene log geometric mean
    usable = np.isfinite(loggeo)               # genes positive in ALL samples
    if usable.sum() >= 1:
        log_ratios = logmat[usable] - loggeo[usable][:, None]
        sf = np.exp(np.median(log_ratios, axis=0))
        if np.all(np.isfinite(sf)) and np.all(sf > 0):
            return sf
    # Fallback: library size (total counts) relative to the mean library.
    totals = matrix.sum(axis=0).astype(float)
    mean_total = totals.mean() if totals.mean() > 0 else 1.0
    sf = totals / mean_total
    sf[sf <= 0] = 1.0
    return sf


# ---------------------------------------------------------------------------
# 3. Per-gene differential expression (variable vs fixed)
# ---------------------------------------------------------------------------

def differential_expression(counts: Dict[str, Dict[str, int]],
                            sample_lineage: Dict[str, str],
                            gene_names: Dict[str, str], *,
                            library_sizes: Optional[Dict[str, int]] = None,
                            sample_genotype: Optional[Dict[str, str]] = None,
                            pseudocount: float = 1.0,
                            min_per_group: int = 2) -> List[Dict]:
    """Variable-vs-fixed DE on depth-normalised, genotype-collapsed gene counts.

    **Depth normalisation.** If `library_sizes` (total mapped reads per flower,
    from `samtools idxstats` over the WHOLE BAM) is given, each flower's counts
    are scaled to counts-per-million-style equivalents — robust even for a
    handful of candidate genes, as long as the library size reflects
    whole-transcriptome depth (true for real RNA-seq). It is total-count, not
    composition-robust median-of-ratios; fine for a candidate direction check,
    not a transcriptome-wide DESeq2 replacement. With no `library_sizes`, a DESeq
    median-of-ratios size factor is estimated from the matrix (needs many genes).

    **Nested replication.** Flowers are technical sub-samples of a *genotype*
    (the biological replicate), exactly as flowers nest in an F1 plant on the ASE
    side. So depth-normalised flowers are first **collapsed to their genotype**
    (mean per genotype) — this makes the fold change weight each genotype equally
    regardless of how many flowers it has (e.g. amphorellae's 6 flowers do not
    outvote k2's 4 + k3's 3), and removes flower-level pseudoreplication from the
    test. `sample_genotype` provides the flower→genotype map; without it each
    flower is treated as its own genotype (legacy behaviour).

    **Test.** A Welch t-test on log2(genotype-mean + pseudocount) when BOTH
    lineages have >= `min_per_group` genotypes (proper biological replication,
    `method=welch_genotype`). If a lineage has only one genotype (so a valid
    across-genotype test is impossible — the common case when one parent is
    sampled), it falls back to a flower-level Welch flagged
    `method=welch_flower_pseudorep`, and the caller warns. log2 fold change is
    oriented variable/fixed; BH over the tested genes.
    """
    samples = sorted(sample_lineage)
    if sample_genotype is None:
        sample_genotype = {s: s for s in samples}   # each flower its own genotype

    gene_ids = sorted(counts)
    matrix = np.array([[counts[g].get(s, 0) for s in samples] for g in gene_ids],
                      dtype=float)
    if library_sizes is not None:
        libs = np.array([max(library_sizes.get(s, 0), 1) for s in samples], dtype=float)
        sf = libs / libs.mean()                # size factor = lib / mean lib
    else:
        sf = size_factors(matrix)
    norm = matrix / sf[None, :]                # genes × flowers, depth-normalised

    # Genotype -> its flower column indices, and genotype -> lineage.
    geno_cols: Dict[str, List[int]] = {}
    geno_lineage: Dict[str, str] = {}
    for i, s in enumerate(samples):
        g = sample_genotype.get(s, s)
        geno_cols.setdefault(g, []).append(i)
        geno_lineage[g] = sample_lineage[s]
    var_genos = sorted(g for g, lin in geno_lineage.items() if lin == "variable")
    fix_genos = sorted(g for g, lin in geno_lineage.items() if lin == "fixed")
    var_flowers = sum(len(geno_cols[g]) for g in var_genos)
    fix_flowers = sum(len(geno_cols[g]) for g in fix_genos)
    genotype_testable = len(var_genos) >= min_per_group and len(fix_genos) >= min_per_group

    rows: List[Dict] = []
    p_for_bh: List[Tuple[int, float]] = []     # (row index, p)
    for gi, gid in enumerate(gene_ids):
        # Collapse flowers -> genotype means (the biological unit).
        geno_mean = {g: float(norm[gi, geno_cols[g]].mean()) for g in geno_cols}
        gv = np.array([geno_mean[g] for g in var_genos])
        gf = np.array([geno_mean[g] for g in fix_genos])
        mean_v = float(gv.mean()) if len(gv) else float("nan")
        mean_f = float(gf.mean()) if len(gf) else float("nan")
        base = float(np.array(list(geno_mean.values())).mean())
        log2fc = math.log2((mean_v + pseudocount) / (mean_f + pseudocount))
        higher = "variable" if log2fc > 0 else ("fixed" if log2fc < 0 else "balanced")

        pval: Optional[float] = None
        method = "none"
        if genotype_testable:
            method = "welch_genotype"
            a, b = np.log2(gv + pseudocount), np.log2(gf + pseudocount)
        else:
            # Fall back to flower level (pseudoreplicated) so triage still has a p.
            method = "welch_flower_pseudorep"
            a = np.log2(norm[gi, [i for g in var_genos for i in geno_cols[g]]] + pseudocount)
            b = np.log2(norm[gi, [i for g in fix_genos for i in geno_cols[g]]] + pseudocount)
        if len(a) >= min_per_group and len(b) >= min_per_group:
            if a.std() == 0 and b.std() == 0 and a.mean() == b.mean():
                pval = 1.0
            else:
                t, p = stats.ttest_ind(a, b, equal_var=False)
                pval = float(p) if np.isfinite(p) else 1.0
        else:
            method = "no_test"

        rows.append({
            "gene_id": gid, "gene_name": gene_names.get(gid, gid),
            "baseMean": round(base, 3),
            "mean_variable": round(mean_v, 3) if not math.isnan(mean_v) else "NA",
            "mean_fixed": round(mean_f, 3) if not math.isnan(mean_f) else "NA",
            "log2FoldChange": round(log2fc, 4),
            "pvalue": pval if pval is not None else "NA",
            "padj": "NA",
            "n_variable": len(var_genos), "n_fixed": len(fix_genos),
            "n_variable_flowers": var_flowers, "n_fixed_flowers": fix_flowers,
            "higher_in": higher, "method": method,
        })
        if pval is not None:
            p_for_bh.append((len(rows) - 1, pval))

    # BH over the testable genes only.
    if p_for_bh:
        qs = st.bh_adjust([p for _, p in p_for_bh])
        for (ri, _), q in zip(p_for_bh, qs):
            rows[ri]["padj"] = round(float(q), 5)

    rows.sort(key=lambda r: (r["padj"] if isinstance(r["padj"], float) else 1.0,
                             -abs(r["log2FoldChange"])))
    return rows


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_parental_de(cfg: CrossConfig, gene_index: GeneIndex, *,
                    gene_ids: Optional[set] = None,
                    min_mapq: int = 20, samtools: str = "samtools",
                    pseudocount: float = 1.0,
                    progress: bool = True) -> List[Dict]:
    """Count parental libraries and compute variable-vs-fixed DE, with a
    pseudoreplication warning when a lineage has a single genotype."""
    n_var_geno = sum(1 for p in cfg.variable_parents if p.flowers)
    n_fix_geno = sum(1 for p in cfg.fixed_parents if p.flowers)
    if min(n_var_geno, n_fix_geno) < 2:
        print(f"  WARNING: parental DE is pseudoreplicated — variable lineage "
              f"has {n_var_geno} genotype(s) with RNA, fixed has {n_fix_geno}. "
              f"Flowers of one genotype are technical replicates, so p-values "
              f"are anticonservative. Trust the fold-change DIRECTION over the "
              f"p-value, or supply a replicated external DE table to `contrast`.")
    counts, sample_lineage, gene_names, library_sizes, sample_genotype = count_parental_expression(
        cfg, gene_index, gene_ids=gene_ids, min_mapq=min_mapq,
        samtools=samtools, progress=progress)
    return differential_expression(counts, sample_lineage, gene_names,
                                   library_sizes=library_sizes,
                                   sample_genotype=sample_genotype,
                                   pseudocount=pseudocount)
