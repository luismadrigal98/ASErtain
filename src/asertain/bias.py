"""Reference-bias mitigation helpers.

When the reference equals one parent (e.g. one kunthii variant), reads carrying
that parent's allele map preferentially and the other allele can be lost
entirely — a serious confound for ASE. Two pre-alignment remedies are supported
here; both require you to RE-ALIGN the F1 reads afterwards:

    * N-masking   replace the reference base at every informative SNP with N so
                  neither allele is favoured at the SNP position.
    * WASP        write per-chromosome SNP files for the WASP remap-and-filter
                  workflow (also corrects read-level multi-mismatch bias).

The statistical alternative (no re-alignment) is `--bias-mode null-shift` with a
balanced-control table, handled in `counting`.
"""
from __future__ import annotations

import gzip
import os
from collections import defaultdict
from typing import Dict, Iterator, List, Set, Tuple

from .genotypes import InformativeSNP


# ---------------------------------------------------------------------------
# FASTA streaming (one record at a time -> bounded memory)
# ---------------------------------------------------------------------------

def _iter_fasta(path: str) -> Iterator[Tuple[str, List[str]]]:
    opener = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    name, seq = None, []
    with opener as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, seq
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, seq


def _write_wrapped(fh, name: str, seq: str, width: int = 60) -> None:
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def nmask_reference(reference_fasta: str, snps: List[InformativeSNP],
                    out_fasta: str) -> Dict[str, int]:
    """Write an N-masked copy of `reference_fasta` at all informative positions.

    Returns {chrom: n_sites_masked}. Re-align the F1 reads to `out_fasta` and
    run `asertain count --bias-mode nmask` against the new BAMs.
    """
    positions: Dict[str, Set[int]] = defaultdict(set)
    for s in snps:
        positions[s.chrom].add(s.pos)            # 1-based

    masked: Dict[str, int] = {}
    with open(out_fasta, "w") as out:
        for name, lines in _iter_fasta(reference_fasta):
            chars = list("".join(lines))
            n = 0
            for pos in positions.get(name, ()):
                idx = pos - 1
                if 0 <= idx < len(chars):
                    chars[idx] = "N"
                    n += 1
            masked[name] = n
            _write_wrapped(out, name, "".join(chars))
    return masked


def write_wasp_snp_files(snps: List[InformativeSNP], out_dir: str) -> List[str]:
    """Write per-chromosome SNP files for WASP find_intersecting_snps.py.

    Format per line: <pos> <ref_allele> <alt_allele>  (gzipped, one file per
    chromosome named <chrom>.snps.txt.gz).
    """
    os.makedirs(out_dir, exist_ok=True)
    by_chrom: Dict[str, List[InformativeSNP]] = defaultdict(list)
    for s in snps:
        by_chrom[s.chrom].append(s)
    written: List[str] = []
    for chrom, group in by_chrom.items():
        path = os.path.join(out_dir, f"{chrom}.snps.txt.gz")
        with gzip.open(path, "wt") as fh:
            for s in sorted(group, key=lambda x: x.pos):
                fh.write(f"{s.pos}\t{s.ref}\t{s.alt}\n")
        written.append(path)
    return written
