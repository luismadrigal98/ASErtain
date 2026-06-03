"""Diagnostic-SNP discovery from exact-parent genotypes.

This is the module that answers the reviewer's central worry — *which SNPs do we
trust?* Instead of pooling a whole species and calling fixation by frequency, we
genotype each named parent individual separately and keep only sites that are:

    1. homozygous and concordant across all variable-species parents,
    2. homozygous in the fixed-species parent(s), and
    3. fixed for *different* alleles between the two species.

A site where the variable-species parents disagree (or one is heterozygous) is
NOT discarded silently — it is reported as 'background_specific', usable only for
the F1s descending from the parent for which it is cleanly diagnostic. That keeps
the robust 'shared' set for combined analysis while preserving power per
background.

Genotype calling tolerates sequencing noise: a parent is called homozygous when the
minor-allele fraction (from AD) is below `maf_threshold`, not strictly 0.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple

from .annotation import GeneIndex
from .config import CrossConfig
from .vcf import SampleCall, Variant, iter_variants

HOM_REF, HOM_ALT, HET, MISSING = "hom_ref", "hom_alt", "het", "missing"


@dataclass
class DiagnosticSNP:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    fixed_allele: str                       # nucleotide fixed in fixed-species parent(s)
    variable_allele_shared: Optional[str]   # nucleotide for 'shared' sites, else None
    diagnostic_class: str                   # 'shared' | 'background_specific'
    backgrounds: List[str]                  # variable-parent backgrounds this SNP serves
    bg_variable_allele: Dict[str, str]      # background -> expected variable nucleotide
    parent_states: Dict[str, str] = field(default_factory=dict)
    parent_depths: Dict[str, Optional[int]] = field(default_factory=dict)
    gene_id: str = "intergenic"
    gene_name: str = "intergenic"
    location: str = "intergenic"

    def variable_is_ref_for(self, background: str) -> bool:
        allele = (self.variable_allele_shared
                  or self.bg_variable_allele.get(background))
        return allele == self.ref


@dataclass
class DiagnoseStats:
    total: int = 0
    biallelic_snp: int = 0
    parents_callable: int = 0
    fixed_in_fixed_species: int = 0
    shared: int = 0
    background_specific: int = 0


def call_state(call: SampleCall, *, min_depth: int, maf_threshold: float) -> str:
    """Classify one parent's genotype, tolerating low-level sequencing noise.

    Prefers allelic depth (AD) when present; otherwise falls back to the GT
    field gated on total depth. Returns HOM_REF / HOM_ALT / HET / MISSING.
    """
    if call.ad is not None and len(call.ad) >= 2:
        ref_d, alt_d = call.ad[0], call.ad[1]
        total = ref_d + alt_d
        if total < min_depth:
            return MISSING
        alt_frac = alt_d / total
        if alt_frac <= maf_threshold:
            return HOM_REF
        if (1 - alt_frac) <= maf_threshold:
            return HOM_ALT
        return HET

    # No AD: use GT, but require depth evidence if we have it.
    if call.dp is not None and call.dp < min_depth:
        return MISSING
    if call.gt is None:
        return MISSING
    a, b = call.gt
    if a == b:
        return HOM_REF if a == "0" else HOM_ALT
    return HET


def _hom_index(state: str) -> Optional[str]:
    return {HOM_REF: "0", HOM_ALT: "1"}.get(state)


def _nuc(idx: str, ref: str, alt: str) -> str:
    return ref if idx == "0" else alt


def classify_site(variant: Variant, cfg: CrossConfig, *,
                  min_depth: int, maf_threshold: float,
                  stats: DiagnoseStats) -> Optional[DiagnosticSNP]:
    """Apply the diagnostic criteria to one biallelic SNP. None if not usable."""
    alt = variant.alt[0]

    # --- fixed-species parent(s): must be homozygous, concordant ----------
    fixed_idx: Optional[str] = None
    parent_states: Dict[str, str] = {}
    parent_depths: Dict[str, Optional[int]] = {}
    for p in cfg.fixed_parents:
        c = variant.call(p.vcf_sample)
        st = call_state(c, min_depth=min_depth, maf_threshold=maf_threshold)
        parent_states[p.name] = st
        parent_depths[p.name] = c.dp
        idx = _hom_index(st)
        if idx is None:
            return None                      # het or missing in fixed parent
        if fixed_idx is None:
            fixed_idx = idx
        elif fixed_idx != idx:
            return None                      # fixed parents disagree
    stats.fixed_in_fixed_species += 1

    # --- variable-species parents: per-parent homozygous index ------------
    var_idx: Dict[str, Optional[str]] = {}
    any_callable = False
    for p in cfg.variable_parents:
        c = variant.call(p.vcf_sample)
        st = call_state(c, min_depth=min_depth, maf_threshold=maf_threshold)
        parent_states[p.name] = st
        parent_depths[p.name] = c.dp
        idx = _hom_index(st)
        var_idx[p.name] = idx
        if idx is not None:
            any_callable = True
    if not any_callable:
        return None

    # Which backgrounds are cleanly diagnostic (variable hom, != fixed)?
    bg_variable_allele: Dict[str, str] = {}
    for name, idx in var_idx.items():
        if idx is not None and idx != fixed_idx:
            bg_variable_allele[name] = _nuc(idx, variant.ref, alt)
    if not bg_variable_allele:
        return None

    fixed_nuc = _nuc(fixed_idx, variant.ref, alt)
    callable_names = [p.name for p in cfg.variable_parents
                      if var_idx[p.name] is not None]
    # 'shared' (safe for every F1 background) requires EVERY variable parent to
    # be callable (homozygous) AND concordant on one diagnostic allele. If any
    # variable parent is heterozygous or missing, the site cannot be trusted for
    # that parent's F1s, so it is background-specific instead.
    shared = (len(callable_names) == len(cfg.variable_parents)
              and set(bg_variable_allele) == set(callable_names)
              and len(set(bg_variable_allele.values())) == 1)

    if shared:
        stats.shared += 1
        variable_shared = next(iter(bg_variable_allele.values()))
        diag_class = "shared"
        backgrounds = sorted(bg_variable_allele)
    else:
        stats.background_specific += 1
        variable_shared = None
        diag_class = "background_specific"
        backgrounds = sorted(bg_variable_allele)

    return DiagnosticSNP(
        chrom=variant.chrom, pos=variant.pos, ref=variant.ref, alt=alt,
        qual=variant.qual,
        fixed_allele=fixed_nuc,
        variable_allele_shared=variable_shared,
        diagnostic_class=diag_class,
        backgrounds=backgrounds,
        bg_variable_allele=bg_variable_allele,
        parent_states=parent_states,
        parent_depths=parent_depths,
    )


def find_diagnostic_snps(cfg: CrossConfig, vcf_path: str, *,
                         min_depth: int = 8,
                         maf_threshold: float = 0.10,
                         min_qual: float = 30.0,
                         snps_only: bool = True,
                         chrom_filter: Optional[str] = None,
                         gene_index: Optional[GeneIndex] = None,
                         ) -> Tuple[List[DiagnosticSNP], DiagnoseStats]:
    """Scan a multi-sample VCF and return diagnostic SNPs + run statistics."""
    out: List[DiagnosticSNP] = []
    stats = DiagnoseStats()
    window = cfg.annotation_window

    for v in iter_variants(vcf_path, snps_only=snps_only,
                           chrom_filter=chrom_filter, min_qual=min_qual):
        stats.total += 1
        if not v.is_biallelic_snp:
            continue
        stats.biallelic_snp += 1
        snp = classify_site(v, cfg, min_depth=min_depth,
                            maf_threshold=maf_threshold, stats=stats)
        if snp is None:
            continue
        if gene_index is not None:
            hit = gene_index.annotate(snp.chrom, snp.pos, window=window)
            snp.gene_id, snp.gene_name, snp.location = (
                hit.gene_id, hit.gene_name, hit.location)
        out.append(snp)
    stats.parents_callable = len(out)
    return out, stats
