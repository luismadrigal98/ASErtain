"""Allele-specific read counting in F1 BAMs, with pluggable reference-bias modes.

Because the F1 reads may be mapped to either parent's reference or to a
third-party reference, *every* bias strategy is selectable by flag and the rest
of the pipeline is agnostic to the choice:

    --bias-mode none       null expectation = 0.5; do nothing
    --bias-mode report     null = 0.5, but record variable_is_ref so downstream
                           can stratify and detect systematic reference pull
    --bias-mode null-shift use a per-SNP balanced-control table to set the null
                           expectation (e.g. F1 gDNA reference-allele fraction)
    --bias-mode wasp       expect WASP-filtered BAMs (or invoke external.wasp_filter)
    --bias-mode nmask      expect BAMs aligned to an N-masked reference

Modes that remove bias at the alignment step (wasp/nmask) simply count and keep
null = 0.5; modes that handle it statistically (null-shift) write a per-record
`null_p` consumed by the `test` stage.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

from . import external
from .config import CrossConfig
from .genotypes import DiagnosticSNP


# ---------------------------------------------------------------------------
# mpileup parsing (refactored from the original ase_read_counter.py)
# ---------------------------------------------------------------------------

def parse_pileup_bases(ref_base: str, bases: str) -> List[str]:
    """Expand an mpileup base string into a list of uppercase base calls."""
    out: List[str] = []
    i = 0
    n = len(bases)
    while i < n:
        c = bases[i]
        if c in ".,":
            out.append(ref_base.upper())
            i += 1
        elif c == "^":          # read-start marker: skip the mapping-quality char
            i += 2
        elif c == "$":          # read-end marker
            i += 1
        elif c in "+-":         # indel: +2AC / -1G  -> skip length+bases
            i += 1
            num = ""
            while i < n and bases[i].isdigit():
                num += bases[i]
                i += 1
            i += int(num) if num else 0
        elif c in "ACGTNacgtn*":
            if c != "*":
                out.append(c.upper())
            i += 1
        else:
            i += 1
    return out


def count_alleles(bam: str, chrom: str, pos: int,
                  variable_allele: str, fixed_allele: str, *,
                  min_mapq: int, min_baseq: int,
                  reference: Optional[str] = None,
                  samtools: str = "samtools") -> Dict[str, int]:
    """Count variable/fixed/other reads at one position via samtools mpileup."""
    empty = {"variable_count": 0, "fixed_count": 0,
             "other_count": 0, "total_depth": 0}
    out = external.samtools_mpileup(
        bam, f"{chrom}:{pos}-{pos}",
        min_mapq=min_mapq, min_baseq=min_baseq,
        reference=reference, samtools=samtools)
    if not out:
        return empty
    fields = out.split("\t")
    if len(fields) < 5:
        return empty
    ref_base, depth, pile = fields[2], int(fields[3]), fields[4]
    calls = parse_pileup_bases(ref_base, pile)
    v = calls.count(variable_allele.upper())
    f = calls.count(fixed_allele.upper())
    return {"variable_count": v, "fixed_count": f,
            "other_count": len(calls) - v - f, "total_depth": depth}


# ---------------------------------------------------------------------------
# Bias-control table (for --bias-mode null-shift)
# ---------------------------------------------------------------------------

def load_control_bias(path: str) -> Dict[Tuple[str, int], float]:
    """Load a balanced-control table: chrom<TAB>pos<TAB>ref_fraction.

    `ref_fraction` is the observed reference-allele fraction at that site in a
    50/50 control (e.g. F1 genomic DNA, or pooled reciprocal crosses). Used to
    set a per-SNP null expectation that absorbs site-specific mapping bias.
    """
    table: Dict[Tuple[str, int], float] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c, p, frac = line.split("\t")[:3]
            table[(c, int(p))] = float(frac)
    return table


def null_expectation(snp: DiagnosticSNP, background: str, bias_mode: str,
                     control: Optional[Dict[Tuple[str, int], float]]) -> float:
    """Expected *variable-allele* fraction under the null for this SNP."""
    if bias_mode == "null-shift" and control is not None:
        ref_frac = control.get((snp.chrom, snp.pos))
        if ref_frac is not None:
            # variable_is_ref ? expect ref_frac : expect 1 - ref_frac
            return ref_frac if snp.variable_is_ref_for(background) else 1 - ref_frac
    return 0.5


# ---------------------------------------------------------------------------
# Applicability of a diagnostic SNP to a given F1 background
# ---------------------------------------------------------------------------

def applicable_allele(snp: DiagnosticSNP, background: str) -> Optional[str]:
    """Return the expected variable-allele nucleotide for this background, or
    None if the SNP is not diagnostic for it."""
    if snp.diagnostic_class == "shared":
        return snp.variable_allele_shared
    return snp.bg_variable_allele.get(background)


# ---------------------------------------------------------------------------
# Main counting driver
# ---------------------------------------------------------------------------

def count_f1_samples(cfg: CrossConfig, snps: List[DiagnosticSNP], *,
                     bias_mode: str = "report",
                     control_table: Optional[str] = None,
                     min_mapq: int = 20, min_baseq: int = 20,
                     min_depth: int = 10,
                     samtools: str = "samtools",
                     progress: bool = True) -> List[Dict]:
    """Count alleles for every (F1 replicate × applicable diagnostic SNP)."""
    reference = cfg.reference.fasta
    control = (load_control_bias(control_table)
               if bias_mode == "null-shift" and control_table else None)
    if bias_mode == "null-shift" and control is None:
        raise ValueError("--bias-mode null-shift requires --control-table")

    records: List[Dict] = []
    for ri, rep in enumerate(cfg.f1, 1):
        if not os.path.exists(rep.bam):
            raise FileNotFoundError(f"BAM for F1 '{rep.name}' not found: {rep.bam}")
        external.ensure_bam_index(rep.bam, samtools=samtools)
        usable = [s for s in snps if applicable_allele(s, rep.background)]
        if progress:
            print(f"  [{ri}/{len(cfg.f1)}] {rep.name} "
                  f"(bg={rep.background}): {len(usable)} usable SNPs", flush=True)
        for snp in usable:
            var_nuc = applicable_allele(snp, rep.background)
            counts = count_alleles(
                rep.bam, snp.chrom, snp.pos, var_nuc, snp.fixed_allele,
                min_mapq=min_mapq, min_baseq=min_baseq,
                reference=reference, samtools=samtools)
            if counts["total_depth"] < min_depth:
                continue
            records.append({
                "f1_sample": rep.name,
                "background": rep.background,
                "chrom": snp.chrom,
                "pos": snp.pos,
                "snp_id": f"{snp.chrom}:{snp.pos}",
                "variable_allele": var_nuc,
                "fixed_allele": snp.fixed_allele,
                "variable_is_ref": snp.variable_is_ref_for(rep.background),
                "variable_count": counts["variable_count"],
                "fixed_count": counts["fixed_count"],
                "other_count": counts["other_count"],
                "total_depth": counts["total_depth"],
                "null_p": round(null_expectation(snp, rep.background,
                                                 bias_mode, control), 5),
                "gene_id": snp.gene_id,
                "gene_name": snp.gene_name,
            })
    return records
