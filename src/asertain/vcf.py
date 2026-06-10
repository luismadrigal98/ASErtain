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
    gt: Optional[Tuple[str, ...]]      # allele indices, e.g. ('0', '1') or
                                       # ('0', '0', '1') for a tetraploid; None if
                                       # missing. Any ploidy is preserved; callers
                                       # that only need presence/absence of each
                                       # allele use the distinct-index set.
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


def _parse_gt(token: str) -> Optional[Tuple[str, ...]]:
    """Parse a GT token to a tuple of allele indices, for ANY ploidy.

    Diploid ('0/1'), haploid ('0'), or polyploid ('0/0/1') are all accepted; the
    genotype-calling logic downstream uses the set of distinct alleles, so this
    keeps polyploid samples usable instead of silently dropping them. Any missing
    allele ('.') makes the whole call missing (None)."""
    sep = "/" if "/" in token else ("|" if "|" in token else None)
    if sep is None:                              # haploid single-allele GT
        return (token,) if token.isdigit() else None
    parts = token.split(sep)
    if not parts or any(p in (".", "") for p in parts):
        return None
    return tuple(parts)


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
            # A missing QUAL ('.') means the caller did not score the site, NOT
            # that it scored zero. Many common callers (GATK GenotypeGVCFs,
            # bcftools in some modes, hard-filtered VCFs) emit '.' for QUAL while
            # still carrying confident genotypes. Coercing it to 0.0 and applying
            # the QUAL filter would silently drop every such site. Store missing
            # QUAL as NaN and exempt it from the filter (only a real numeric QUAL
            # below the threshold is dropped).
            raw_qual = f[5]
            if raw_qual in (".", ""):
                qual = float("nan")
            else:
                try:
                    qual = float(raw_qual)
                except ValueError:
                    qual = float("nan")
            if qual == qual and qual < min_qual:   # qual==qual is False for NaN
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
