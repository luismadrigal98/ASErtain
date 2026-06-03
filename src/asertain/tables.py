"""Read/write contracts for the TSV files passed between pipeline stages.

Centralising the column layouts here keeps `diagnose`, `count`, `test`, and
`contrast` decoupled: each stage only needs to agree on these readers/writers.
All tables are tab-separated with a single header line and a leading comment
block (#) describing provenance.
"""
from __future__ import annotations

import csv
from typing import Dict, Iterator, List, Optional

from .genotypes import DiagnosticSNP


# ---------------------------------------------------------------------------
# Diagnostic SNPs
# ---------------------------------------------------------------------------

_DIAG_COLS = [
    "chrom", "pos", "ref", "alt", "qual",
    "fixed_allele", "variable_allele_shared", "diagnostic_class",
    "backgrounds", "bg_variable_allele",
    "parent_states", "gene_id", "gene_name", "location",
]


def _encode_bg_alleles(d: Dict[str, str]) -> str:
    return ";".join(f"{k}:{v}" for k, v in sorted(d.items())) or "."


def _decode_bg_alleles(s: str) -> Dict[str, str]:
    if not s or s == ".":
        return {}
    return dict(kv.split(":", 1) for kv in s.split(";") if ":" in kv)


def write_diagnostic_snps(snps: List[DiagnosticSNP], path: str,
                          *, source_vcf: str = "") -> None:
    with open(path, "w", newline="") as fh:
        fh.write(f"# ASErtain diagnostic SNPs\n# source_vcf: {source_vcf}\n")
        w = csv.writer(fh, delimiter="\t")
        w.writerow(_DIAG_COLS)
        for s in snps:
            states = ";".join(f"{k}={v}" for k, v in sorted(s.parent_states.items()))
            w.writerow([
                s.chrom, s.pos, s.ref, s.alt, f"{s.qual:.2f}",
                s.fixed_allele, s.variable_allele_shared or ".",
                s.diagnostic_class, ",".join(s.backgrounds),
                _encode_bg_alleles(s.bg_variable_allele),
                states, s.gene_id, s.gene_name, s.location,
            ])


def read_diagnostic_snps(path: str) -> List[Dict]:
    """Return diagnostic SNPs as plain dicts (with parsed compound fields)."""
    rows: List[Dict] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            header = line.rstrip("\n").split("\t")
            break
        r = csv.DictReader((l for l in fh if not l.startswith("#")),
                           fieldnames=header, delimiter="\t")
        for row in r:
            row["pos"] = int(row["pos"])
            row["backgrounds"] = (row["backgrounds"].split(",")
                                  if row["backgrounds"] else [])
            row["bg_variable_allele"] = _decode_bg_alleles(row["bg_variable_allele"])
            rows.append(row)
    return rows


def write_bed(snps: List[DiagnosticSNP], path: str) -> None:
    with open(path, "w") as fh:
        fh.write("# chrom\tstart\tend\tname\tscore\tstrand\n")
        for s in snps:
            name = f"{s.chrom}:{s.pos}_{s.fixed_allele}|{s.diagnostic_class}"
            fh.write(f"{s.chrom}\t{s.pos - 1}\t{s.pos}\t{name}\t{s.qual:.0f}\t.\n")


# ---------------------------------------------------------------------------
# Allele counts (count stage output)
# ---------------------------------------------------------------------------

COUNT_COLS = [
    "f1_sample", "background", "chrom", "pos", "snp_id",
    "variable_allele", "fixed_allele", "variable_is_ref",
    "variable_count", "fixed_count", "other_count", "total_depth",
    "null_p", "gene_id", "gene_name",
]


def write_allele_counts(records: List[Dict], path: str,
                        *, bias_mode: str = "") -> None:
    with open(path, "w", newline="") as fh:
        fh.write(f"# ASErtain allele counts\n# bias_mode: {bias_mode}\n")
        w = csv.DictWriter(fh, fieldnames=COUNT_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for rec in records:
            w.writerow(rec)


def read_allele_counts(path: str) -> List[Dict]:
    rows: List[Dict] = []
    with open(path) as fh:
        data_lines = (l for l in fh if not l.startswith("#"))
        r = csv.DictReader(data_lines, delimiter="\t")
        for row in r:
            for c in ("pos", "variable_count", "fixed_count",
                      "other_count", "total_depth"):
                row[c] = int(row[c])
            row["null_p"] = float(row["null_p"])
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Generic gene-level table (test / contrast)
# ---------------------------------------------------------------------------

def write_table(records: List[Dict], cols: List[str], path: str,
                *, comment: str = "") -> None:
    with open(path, "w", newline="") as fh:
        if comment:
            fh.write(f"# {comment}\n")
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for rec in records:
            w.writerow(rec)


def read_table(path: str) -> List[Dict]:
    with open(path) as fh:
        data = (l for l in fh if not l.startswith("#"))
        return list(csv.DictReader(data, delimiter="\t"))
