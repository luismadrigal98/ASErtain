"""Gene annotation from GTF/GFF3, used to assign diagnostic SNPs to genes.

Lifted from the original ase scripts and condensed. Supports both GTF
(key "value") and GFF3 (key=value) attribute encodings, and an optional
window to capture UTR-proximal sites.
"""
from __future__ import annotations

import gzip
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class Gene:
    gene_id: str
    gene_name: str
    biotype: str
    start: int
    end: int
    strand: str


@dataclass
class Hit:
    gene_id: str = "intergenic"
    gene_name: str = "intergenic"
    biotype: str = "intergenic"
    location: str = "intergenic"   # genic | upstream | downstream | intergenic


class GeneIndex:
    """Per-chromosome list of genes with simple linear overlap lookup."""

    def __init__(self, by_chrom: Dict[str, List[Gene]]):
        self._by_chrom = by_chrom

    @property
    def n_genes(self) -> int:
        return sum(len(v) for v in self._by_chrom.values())

    def genes(self) -> List[Gene]:
        """Every gene in the index (used by the parental-expression stage to
        define per-gene counting intervals)."""
        return [g for genes in self._by_chrom.values() for g in genes]

    def iter_genes(self):
        """Yield (chrom, Gene) for every gene — chrom is needed to build the
        samtools counting region."""
        for chrom, genes in self._by_chrom.items():
            for g in genes:
                yield chrom, g

    @classmethod
    def from_file(cls, path: str) -> "GeneIndex":
        is_gff3 = path.endswith((".gff3", ".gff", ".gff3.gz", ".gff.gz"))
        opener = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
        by_chrom: Dict[str, List[Gene]] = defaultdict(list)
        with opener as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "gene":
                    continue
                attrs = _parse_attrs(f[8], is_gff3)
                gid = (attrs.get("ID") or attrs.get("gene_id")
                       or attrs.get("Name"))
                if not gid:
                    continue
                by_chrom[f[0]].append(Gene(
                    gene_id=gid,
                    gene_name=attrs.get("Name") or attrs.get("gene_name") or gid,
                    biotype=(attrs.get("biotype") or attrs.get("gene_biotype")
                             or attrs.get("gene_type") or "unknown"),
                    start=int(f[3]), end=int(f[4]), strand=f[6],
                ))
        return cls(dict(by_chrom))

    def annotate(self, chrom: str, pos: int, window: int = 0) -> Hit:
        hit = Hit()
        for g in self._by_chrom.get(chrom, ()):
            gs, ge = max(1, g.start - window), g.end + window
            if gs <= pos <= ge:
                hit.gene_id, hit.gene_name, hit.biotype = (
                    g.gene_id, g.gene_name, g.biotype)
                if g.start <= pos <= g.end:
                    hit.location = "genic"
                elif g.strand == "+":
                    hit.location = "upstream" if pos < g.start else "downstream"
                else:
                    hit.location = "downstream" if pos < g.start else "upstream"
                break
        return hit

def _parse_attrs(field: str, is_gff3: bool) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if is_gff3:
        for kv in field.split(";"):
            if "=" in kv:
                k, v = kv.strip().split("=", 1)
                out[k] = v
    else:
        for attr in field.split(";"):
            attr = attr.strip()
            if not attr:
                continue
            if '"' in attr:
                key = attr.split()[0]
                out[key] = attr.split('"')[1]
            else:
                bits = attr.split()
                if len(bits) >= 2:
                    out[bits[0]] = bits[1]
    return out
