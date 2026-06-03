"""Minimal, dependency-free VCF reading tailored to ASE genotype extraction.

Refactored from the shared helpers in the original ase_diagnostic_snps.py /
ase_heterozygous_pipeline.py scripts. We parse only what we need: per-sample
genotype indices and allelic/total depth at biallelic SNP sites.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Tuple


def open_text(path: str):
    """Open a plain or gzipped text file transparently."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_header(handle) -> Tuple[List[str], Dict[str, int]]:
    """Consume `##`/`#CHROM` lines, returning (sample_names, name->col_index)."""
    for line in handle:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            return samples, {s: i for i, s in enumerate(samples)}
    raise ValueError("No #CHROM header line found in VCF")


def sample_names(path: str) -> List[str]:
    """Return just the sample names from a VCF header."""
    with open_text(path) as fh:
        samples, _ = read_header(fh)
    return samples


@dataclass
class SampleCall:
    """Parsed per-sample fields at one site."""
    gt: Optional[Tuple[str, str]]      # allele indices, e.g. ('0', '1'); None if missing
    ad: Optional[Tuple[int, ...]]      # allelic depths (ref, alt, ...) if present
    dp: Optional[int]                  # total depth if available

    @property
    def is_missing(self) -> bool:
        return self.gt is None


@dataclass
class Variant:
    """A single VCF record with lazily-parsed sample calls."""
    chrom: str
    pos: int
    ref: str
    alt: List[str]
    qual: float
    fmt: List[str]
    raw_samples: List[str]
    sample_indices: Dict[str, int]

    @property
    def is_biallelic_snp(self) -> bool:
        return (len(self.ref) == 1 and len(self.alt) == 1
                and len(self.alt[0]) == 1)

    def call(self, sample: str) -> SampleCall:
        """Parse the FORMAT fields for one sample by name."""
        idx = self.sample_indices.get(sample)
        if idx is None or idx >= len(self.raw_samples):
            return SampleCall(None, None, None)
        return _parse_sample_field(self.fmt, self.raw_samples[idx])


def _parse_gt(token: str) -> Optional[Tuple[str, str]]:
    sep = "/" if "/" in token else ("|" if "|" in token else None)
    if sep is None:
        return None
    parts = token.split(sep)
    if len(parts) != 2 or "." in parts:
        return None
    return parts[0], parts[1]


def _parse_sample_field(fmt: List[str], field: str) -> SampleCall:
    vals = field.split(":")
    tag = {t: i for i, t in enumerate(fmt)}

    gt = None
    if "GT" in tag and tag["GT"] < len(vals):
        gt = _parse_gt(vals[tag["GT"]])

    ad = None
    if "AD" in tag and tag["AD"] < len(vals):
        try:
            ad = tuple(int(x) for x in vals[tag["AD"]].split(",")
                       if x not in (".", ""))
        except ValueError:
            ad = None

    dp = None
    if "DP" in tag and tag["DP"] < len(vals):
        try:
            dp = int(vals[tag["DP"]])
        except ValueError:
            dp = None
    if dp is None and ad is not None:
        dp = sum(ad)

    return SampleCall(gt, ad, dp)


def iter_variants(path: str, *, snps_only: bool = False,
                  chrom_filter: Optional[str] = None,
                  min_qual: float = 0.0) -> Iterator[Variant]:
    """Stream VCF records as Variant objects, applying cheap site filters."""
    with open_text(path) as fh:
        samples, indices = read_header(fh)
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom = f[0]
            if chrom_filter and chrom_filter not in chrom:
                continue
            try:
                qual = float(f[5])
            except ValueError:
                qual = 0.0
            if qual < min_qual:
                continue
            alt = f[4].split(",")
            ref = f[3]
            if snps_only and (len(ref) > 1 or any(len(a) > 1 for a in alt)
                              or len(alt) > 1):
                continue
            yield Variant(
                chrom=chrom, pos=int(f[1]), ref=ref, alt=alt,
                qual=qual, fmt=f[8].split(":"),
                raw_samples=f[9:], sample_indices=indices,
            )
