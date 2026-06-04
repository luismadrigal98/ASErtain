"""Allele-specific read counting in F1 BAMs, per flower, with bias modes.

Counting is done at the *flower* (RNA sample) level; flowers are later collapsed
to their plant in the test stage. For each flower we only consider SNPs that are
informative for that flower's plant, and we expect that plant's variable/fixed
alleles.

Reference-bias modes (all flag-selectable, because the reads may be mapped to
one parent, the other, or a third reference):

    none       null = 0.5, no bookkeeping
    report     null = 0.5, record variable_is_ref so systematic pull is visible
    null-shift per-SNP null from a balanced-control table
    nmask      expect BAMs realigned to an N-masked reference (see `asertain
               mask-reference`); count normally, null = 0.5
    wasp       expect WASP-filtered BAMs (see external.wasp_filter); null = 0.5
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

from . import external
from .config import CrossConfig
from .genotypes import InformativeSNP


# ---------------------------------------------------------------------------
# mpileup parsing (RNA-aware: handles spliced-read reference skips)
# ---------------------------------------------------------------------------

def parse_pileup_bases(ref_base: str, bases: str) -> List[str]:
    """Expand an mpileup base string into a list of uppercase base calls.

    Handles match (./,), mismatches, insertions/deletions notation, read
    start/end markers, deletion placeholders (*), and — important for RNA-seq —
    reference-skip markers (< >) produced by reads spanning introns.
    """
    out: List[str] = []
    i, n = 0, len(bases)
    while i < n:
        c = bases[i]
        if c in ".,":
            out.append(ref_base.upper())
            i += 1
        elif c == "^":            # read-start marker: skip the mapping-quality char
            i += 2
        elif c in "$*<>":         # read-end / deletion / reference-skip: no base call
            i += 1
        elif c in "+-":           # indel: +2AC / -1G  -> skip length + bases
            i += 1
            num = ""
            while i < n and bases[i].isdigit():
                num += bases[i]
                i += 1
            i += int(num) if num else 0
        elif c in "ACGTNacgtn":
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
    ref_base, raw_depth, pile = fields[2], int(fields[3]), fields[4]
    calls = parse_pileup_bases(ref_base, pile)
    v = calls.count(variable_allele.upper())
    f = calls.count(fixed_allele.upper())
    # total_depth is the number of usable base calls (excludes intron skips <>,
    # deletions * and indel padding that parse_pileup_bases drops). The raw
    # mpileup column-4 depth would over-count for spliced RNA reads (audit C1).
    return {"variable_count": v, "fixed_count": f,
            "other_count": len(calls) - v - f,
            "total_depth": len(calls), "raw_depth": raw_depth}


# ---------------------------------------------------------------------------
# Bias-control table (for --bias-mode null-shift)
# ---------------------------------------------------------------------------

def load_control_bias(path: str) -> Dict[Tuple[str, int], float]:
    """Load chrom<TAB>pos<TAB>ref_fraction from a balanced control."""
    table: Dict[Tuple[str, int], float] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c, p, frac = line.split("\t")[:3]
            table[(c, int(p))] = float(frac)
    return table


def null_expectation(snp: InformativeSNP, plant: str, bias_mode: str,
                     control: Optional[Dict[Tuple[str, int], float]]) -> float:
    """Expected *variable-allele* fraction under the null for this SNP."""
    if bias_mode == "null-shift" and control is not None:
        ref_frac = control.get((snp.chrom, snp.pos))
        if ref_frac is not None:
            return ref_frac if snp.variable_is_ref(plant) else 1 - ref_frac
    return 0.5


# ---------------------------------------------------------------------------
# Counting driver
# ---------------------------------------------------------------------------

def count_flowers(cfg: CrossConfig, snps: List[InformativeSNP], *,
                  bias_mode: str = "report",
                  control_table: Optional[str] = None,
                  min_mapq: int = 20, min_baseq: int = 20,
                  min_depth: int = 10,
                  samtools: str = "samtools",
                  progress: bool = True) -> List[Dict]:
    """Count alleles for every (flower × SNP informative for its plant)."""
    reference = cfg.reference.fasta
    if not reference:
        raise ValueError(
            "Counting requires a reference FASTA (config reference.fasta). "
            "Without it, mpileup '.'/',' match symbols cannot be resolved to a "
            "base and the reference allele would be silently miscounted.")
    control = (load_control_bias(control_table)
               if bias_mode == "null-shift" and control_table else None)
    if bias_mode == "null-shift" and control is None:
        raise ValueError("--bias-mode null-shift requires --control-table")

    # Index SNPs by plant for quick per-flower lookup.
    snps_by_plant: Dict[str, List[InformativeSNP]] = {}
    for pl in cfg.f1_plants:
        snps_by_plant[pl.name] = [s for s in snps if pl.name in s.per_plant]

    records: List[Dict] = []
    for pi, pl in enumerate(cfg.f1_plants, 1):
        usable = snps_by_plant[pl.name]
        if progress:
            print(f"  [{pi}/{len(cfg.f1_plants)}] plant {pl.name} "
                  f"(bg={pl.bg}, {len(pl.flowers)} flowers): "
                  f"{len(usable)} informative SNPs", flush=True)
        for fl in pl.flowers:
            if not os.path.exists(fl.bam):
                raise FileNotFoundError(
                    f"BAM for flower '{fl.name}' (plant {pl.name}) not found: {fl.bam}")
            external.ensure_bam_index(fl.bam, samtools=samtools)
            for snp in usable:
                pa = snp.per_plant[pl.name]
                counts = count_alleles(
                    fl.bam, snp.chrom, snp.pos, pa.variable, pa.fixed,
                    min_mapq=min_mapq, min_baseq=min_baseq,
                    reference=reference, samtools=samtools)
                # Filter on allele-bearing reads, not raw/usable depth: a SNP is
                # only informative through reads that carry one of the two alleles.
                if counts["variable_count"] + counts["fixed_count"] < min_depth:
                    continue
                records.append({
                    "flower": fl.name,
                    "plant": pl.name,
                    "background": pl.bg,
                    "chrom": snp.chrom,
                    "pos": snp.pos,
                    "snp_id": f"{snp.chrom}:{snp.pos}",
                    "variable_allele": pa.variable,
                    "fixed_allele": pa.fixed,
                    "variable_is_ref": snp.variable_is_ref(pl.name),
                    "tier": pa.tier,
                    "variable_count": counts["variable_count"],
                    "fixed_count": counts["fixed_count"],
                    "other_count": counts["other_count"],
                    "total_depth": counts["total_depth"],
                    "null_p": round(null_expectation(snp, pl.name,
                                                     bias_mode, control), 5),
                    "gene_id": snp.gene_id,
                    "gene_name": snp.gene_name,
                })
    return records
