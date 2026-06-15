"""Read-backed haplotype counting (phASER-style) — independent reads per gene.

The per-SNP pileup counter treats every informative SNP as a separate
observation, but SNPs in one gene are **not** independent: a single read/fragment
that spans several SNPs is counted once at *each* of them, so the same molecule
is double-counted and the effective depth (and the per-plant p-value) is
inflated. The beta-binomial absorbs SNP-to-SNP dispersion but not this read-level
correlation.

The fix is to count at the level of the **read**, not the SNP. For each fragment
overlapping a gene's informative SNPs we read the allele it carries at every such
SNP (walking the CIGAR), assign the *whole fragment* to the variable or fixed
haplotype, and count it **once** per gene. The result is a single (K, N) of
genuinely independent reads per gene × flower — so the downstream per-plant test
is a clean binomial over independent reads, with no SNP pseudo-replication.

A fragment whose informative SNPs disagree (it carries *both* a variable and a
fixed allele — a sequencing error, a mis-phased SNP, or a recombinant molecule)
is conservatively called **ambiguous** and excluded from K/N (kept in
`other_count` for QC). Paired mates share a QNAME, so their votes are pooled and
the fragment is counted once, which also avoids double-counting the mate overlap.

This needs no reference FASTA (the read's own base is used, not an mpileup match
symbol) and is compatible with the `nmask` / `wasp` de-biasing BAMs. It does
require the informative SNPs to be gene-annotated (it groups by gene).
"""
from __future__ import annotations

import re
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

from . import external
from .config import CrossConfig
from .genotypes import InformativeSNP

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


# ---------------------------------------------------------------------------
# CIGAR walk: the read's base at each target reference position
# ---------------------------------------------------------------------------

def read_bases_at(pos: int, cigar: str, seq: str, qual: str,
                  targets: Sequence[int]) -> Dict[int, Tuple[str, int]]:
    """Bases a read carries at the requested 1-based reference positions.

    `pos` is the read's 1-based leftmost mapped position, `cigar`/`seq`/`qual`
    its SAM fields. Returns {ref_pos: (BASE, base_quality)} only for targets the
    read actually covers with an aligned base (M/=/X); positions falling in a
    deletion or intron skip (D/N) are correctly absent.
    """
    out: Dict[int, Tuple[str, int]] = {}
    if not targets or cigar == "*" or seq == "*":
        return out
    tset = [t for t in targets if t >= pos]
    if not tset:
        return out
    ref = pos          # 1-based reference cursor
    qi = 0             # 0-based query cursor
    has_qual = qual and qual != "*"
    for m in _CIGAR.finditer(cigar):
        length = int(m.group(1))
        op = m.group(2)
        if op in ("M", "=", "X"):
            end = ref + length
            for t in tset:
                if ref <= t < end:
                    qp = qi + (t - ref)
                    if 0 <= qp < len(seq):
                        bq = (ord(qual[qp]) - 33) if has_qual and qp < len(qual) else 60
                        out[t] = (seq[qp].upper(), bq)
            ref += length
            qi += length
        elif op in ("D", "N"):
            ref += length
        elif op in ("I", "S"):
            qi += length
        # H, P consume neither reference nor query
    return out


# ---------------------------------------------------------------------------
# Per-gene haplotype read counting in one BAM
# ---------------------------------------------------------------------------

def count_gene_haplotypes(bam: str, chrom: str, span_start: int, span_end: int,
                          snps: Sequence[Tuple[int, str, str]], *,
                          min_mapq: int = 20, min_baseq: int = 20,
                          samtools: str = "samtools") -> Tuple[int, int, int, int]:
    """Assign each fragment over [span_start, span_end] to variable/fixed/ambiguous.

    `snps` is a list of (pos, variable_allele, fixed_allele) for this gene/plant.
    Returns (variable_reads, fixed_reads, ambiguous_reads, n_fragments_seen).
    """
    region = f"{chrom}:{span_start}-{span_end}"
    targets = [p for p, _, _ in snps]
    allele = {p: (v.upper(), f.upper()) for p, v, f in snps}
    sam = external.samtools_view(bam, region, min_mapq=min_mapq, samtools=samtools)

    frag: Dict[str, List[int]] = defaultdict(lambda: [0, 0, 0])  # qname -> [var, fix, other]
    for line in sam.splitlines():
        if not line or line[0] == "@":
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        qname, pos, cigar, seq, qual = f[0], int(f[3]), f[5], f[9], f[10]
        bases = read_bases_at(pos, cigar, seq, qual, targets)
        if not bases:
            continue
        slot = frag[qname]
        for t, (base, bq) in bases.items():
            if bq < min_baseq:
                continue
            va, fx = allele[t]
            if base == va:
                slot[0] += 1
            elif base == fx:
                slot[1] += 1
            else:
                slot[2] += 1

    var = fix = amb = 0
    for V, F, _O in frag.values():
        if V > 0 and F == 0:
            var += 1
        elif F > 0 and V == 0:
            fix += 1
        elif V > 0 and F > 0:
            amb += 1                 # carries both haplotypes -> ambiguous
        # V==0 and F==0: covered SNP sites but matched neither allele -> drop
    return var, fix, amb, len(frag)


# ---------------------------------------------------------------------------
# Driver: per (flower × gene) haplotype counts, emitted as one pseudo-SNP/gene
# ---------------------------------------------------------------------------

def count_flowers_haplotype(cfg: CrossConfig, snps: List[InformativeSNP], *,
                            min_mapq: int = 20, min_baseq: int = 20,
                            min_depth: int = 10,
                            filter_secondary: bool = False,
                            samtools: str = "samtools",
                            progress: bool = True) -> List[Dict]:
    """Read-backed counts for every (flower × gene informative for its plant).

    Emits records in the same schema as the pileup counter, but one **per gene**
    (snp_id = ``hap:<gene_id>``, tier = ``haplotype``), so the existing test
    stage consumes them unchanged — one independent (K, N) per gene × plant.
    """
    import os

    # Group each plant's gene-annotated informative SNPs by gene.
    plant_genes: Dict[str, Dict[Tuple[str, str, str], List[Tuple]]] = {}
    for pl in cfg.f1_plants:
        by_gene: Dict[Tuple[str, str, str], List[Tuple]] = defaultdict(list)
        for s in snps:
            if pl.name not in s.per_plant:
                continue
            if s.gene_id in (None, "", "intergenic"):
                continue
            pa = s.per_plant[pl.name]
            by_gene[(s.gene_id, s.gene_name, s.chrom)].append((s.pos, pa, s))
        plant_genes[pl.name] = by_gene

    if not any(plant_genes.values()):
        raise ValueError(
            "Read-backed haplotype counting needs gene-annotated informative "
            "SNPs (none found). Provide a `gtf` in the config so SNPs are "
            "assigned to genes, or use --counter pileup.")

    records: List[Dict] = []
    for pi, pl in enumerate(cfg.f1_plants, 1):
        by_gene = plant_genes[pl.name]
        if progress:
            print(f"  [{pi}/{len(cfg.f1_plants)}] plant {pl.name} "
                  f"(bg={pl.bg}, {len(pl.flowers)} flowers): "
                  f"{len(by_gene)} genes with informative SNPs", flush=True)
        for fl in pl.flowers:
            if not os.path.exists(fl.bam):
                raise FileNotFoundError(
                    f"BAM for flower '{fl.name}' (plant {pl.name}) not found: {fl.bam}")
            bam = external.prepare_bam(
                fl.bam, filter_secondary=filter_secondary, samtools=samtools)
            external.ensure_bam_index(bam, samtools=samtools)
            for (gene_id, gene_name, chrom), gsnps in by_gene.items():
                positions = [p for p, _, _ in gsnps]
                triples = [(p, pa.variable, pa.fixed) for p, pa, _ in gsnps]
                var, fix, amb, _n = count_gene_haplotypes(
                    bam, chrom, min(positions), max(positions), triples,
                    min_mapq=min_mapq, min_baseq=min_baseq, samtools=samtools)
                if var + fix < min_depth:
                    continue
                # Gene-level variable_is_ref: majority over its informative SNPs.
                ref_votes = sum(1 for _, pa, s in gsnps if pa.variable == s.ref)
                var_is_ref = (ref_votes * 2 >= len(gsnps))
                records.append({
                    "flower": fl.name, "plant": pl.name, "background": pl.bg,
                    "chrom": chrom, "pos": min(positions),
                    "snp_id": f"hap:{gene_id}",
                    "variable_allele": "hap", "fixed_allele": "hap",
                    "variable_is_ref": var_is_ref, "tier": "haplotype",
                    "variable_count": var, "fixed_count": fix,
                    "other_count": amb, "total_depth": var + fix + amb,
                    "null_p": 0.5, "gene_id": gene_id, "gene_name": gene_name,
                    "n_hap_snps": len(gsnps),
                })
    return records
