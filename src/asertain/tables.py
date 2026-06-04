"""Read/write contracts for the TSV files passed between pipeline stages.

Centralising the column layouts here keeps `diagnose`, `count`, `test`, and
`contrast` decoupled. All tables are tab-separated with a single header line and
a leading comment block (#) describing provenance.
"""
from __future__ import annotations

import csv
from typing import Dict, List

from .genotypes import InformativeSNP, PlantAllele


# ---------------------------------------------------------------------------
# Informative SNPs (diagnose stage)
# ---------------------------------------------------------------------------

_SNP_COLS = [
    "chrom", "pos", "ref", "alt", "qual",
    "classification", "n_plants", "backgrounds", "per_plant",
    "gene_id", "gene_name", "location",
]


def _encode_per_plant(d: Dict[str, PlantAllele]) -> str:
    # plant:variableNuc/fixedNuc/tier ; ...
    return ";".join(f"{name}:{pa.variable}/{pa.fixed}/{pa.tier}"
                    for name, pa in sorted(d.items())) or "."


def _decode_per_plant(s: str) -> Dict[str, PlantAllele]:
    out: Dict[str, PlantAllele] = {}
    if not s or s == ".":
        return out
    for tok in s.split(";"):
        name, rest = tok.split(":", 1)
        var, fix, tier = rest.split("/")
        out[name] = PlantAllele(variable=var, fixed=fix, tier=tier)
    return out


def write_informative_snps(snps: List[InformativeSNP], path: str,
                           *, source_vcf: str = "") -> None:
    with open(path, "w", newline="") as fh:
        fh.write(f"# ASErtain informative SNPs\n# source_vcf: {source_vcf}\n")
        w = csv.writer(fh, delimiter="\t")
        w.writerow(_SNP_COLS)
        for s in snps:
            w.writerow([
                s.chrom, s.pos, s.ref, s.alt, f"{s.qual:.2f}",
                s.classification, len(s.per_plant), ",".join(s.backgrounds),
                _encode_per_plant(s.per_plant),
                s.gene_id, s.gene_name, s.location,
            ])


def read_informative_snps(path: str) -> List[InformativeSNP]:
    out: List[InformativeSNP] = []
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            out.append(InformativeSNP(
                chrom=row["chrom"], pos=int(row["pos"]),
                ref=row["ref"], alt=row["alt"], qual=float(row["qual"]),
                per_plant=_decode_per_plant(row["per_plant"]),
                classification=row["classification"],
                backgrounds=row["backgrounds"].split(",") if row["backgrounds"] else [],
                gene_id=row.get("gene_id", "intergenic"),
                gene_name=row.get("gene_name", "intergenic"),
                location=row.get("location", "intergenic"),
            ))
    return out


def write_bed(snps: List[InformativeSNP], path: str) -> None:
    with open(path, "w") as fh:
        fh.write("# chrom\tstart\tend\tname\tscore\tstrand\n")
        for s in snps:
            name = f"{s.chrom}:{s.pos}_{s.classification}"
            fh.write(f"{s.chrom}\t{s.pos - 1}\t{s.pos}\t{name}\t{s.qual:.0f}\t.\n")


# ---------------------------------------------------------------------------
# Allele counts (count stage)
# ---------------------------------------------------------------------------

COUNT_COLS = [
    "flower", "plant", "background", "chrom", "pos", "snp_id",
    "variable_allele", "fixed_allele", "variable_is_ref", "tier",
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
        r = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
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
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))
