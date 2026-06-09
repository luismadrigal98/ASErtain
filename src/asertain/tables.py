"""Read/write contracts for the TSV files passed between pipeline stages.

Centralising the column layouts here keeps `diagnose`, `count`, `test`, and
`contrast` decoupled. All tables are tab-separated with a single header line and
a leading comment block (#) describing provenance.

Every writer stamps the user's lineage labels into the header and rewrites the
canonical ``variable``/``fixed`` column names (and ``direction`` values) to those
labels (see :mod:`asertain.labels`); every reader recovers the labels from the
header and maps them back to canonical, so the rest of the code only ever sees
``variable``/``fixed``. Files with no label header are read verbatim.
"""
from __future__ import annotations

import csv
from typing import Dict, List, Optional, Sequence

from .genotypes import InformativeSNP, PlantAllele
from .labels import DIRECTION_COLS, Labels, format_comment, parse_comment


# ---------------------------------------------------------------------------
# Label-aware row helpers
# ---------------------------------------------------------------------------

def _write_labeled(fh, records: Sequence[Dict], cols: List[str], labels: Labels,
                   *, direction_cols: Sequence[str] = DIRECTION_COLS) -> None:
    """Write `records` (canonical keys) as a labelled TSV body."""
    disp_cols = [labels.to_display(c) for c in cols]
    dset = set(direction_cols)
    w = csv.DictWriter(fh, fieldnames=disp_cols, delimiter="\t",
                       extrasaction="ignore")
    w.writeheader()
    for rec in records:
        row = {}
        for c in cols:
            if c not in rec:
                continue
            v = rec[c]
            if c in dset:
                v = labels.value_to_display(str(v))
            row[labels.to_display(c)] = v
        w.writerow(row)


def _canon_rows(rows: List[Dict], labels: Optional[Labels],
                *, direction_cols: Sequence[str] = DIRECTION_COLS) -> List[Dict]:
    """Map labelled column names/values in `rows` back to canonical."""
    if labels is None or labels.is_default:
        return rows
    dset = set(direction_cols)
    out: List[Dict] = []
    for row in rows:
        nr: Dict = {}
        for k, v in row.items():
            ck = labels.to_canonical(k)
            if ck in dset and isinstance(v, str):
                v = labels.value_to_canonical(v)
            nr[ck] = v
        out.append(nr)
    return out


def _read_dicts(path: str) -> List[Dict]:
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


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
                           *, source_vcf: str = "",
                           labels: Labels = Labels()) -> None:
    with open(path, "w", newline="") as fh:
        fh.write(f"# ASErtain informative SNPs\n# source_vcf: {source_vcf}\n")
        fh.write(format_comment(labels))
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
    "n_hap_snps",      # read-backed counter only: #SNPs phased into this gene's reads
]


def write_allele_counts(records: List[Dict], path: str,
                        *, bias_mode: str = "",
                        labels: Labels = Labels()) -> None:
    with open(path, "w", newline="") as fh:
        fh.write(f"# ASErtain allele counts\n# bias_mode: {bias_mode}\n")
        fh.write(format_comment(labels))
        _write_labeled(fh, records, COUNT_COLS, labels)


def read_allele_counts(path: str) -> List[Dict]:
    labels = parse_comment(path)
    rows = _canon_rows(_read_dicts(path), labels)
    for row in rows:
        for c in ("pos", "variable_count", "fixed_count",
                  "other_count", "total_depth"):
            row[c] = int(row[c])
        row["null_p"] = float(row["null_p"])
        if row.get("n_hap_snps") not in (None, ""):   # read-backed counter only
            row["n_hap_snps"] = int(row["n_hap_snps"])
    return rows


# ---------------------------------------------------------------------------
# Generic gene-level table (test / contrast / DE)
# ---------------------------------------------------------------------------

def write_table(records: List[Dict], cols: List[str], path: str,
                *, comment: str = "", labels: Labels = Labels(),
                direction_cols: Sequence[str] = DIRECTION_COLS) -> None:
    with open(path, "w", newline="") as fh:
        if comment:
            fh.write(f"# {comment}\n")
        fh.write(format_comment(labels))
        _write_labeled(fh, records, cols, labels, direction_cols=direction_cols)


def read_table(path: str, *, canonicalize: bool = True) -> List[Dict]:
    """Read a TSV to a list of dicts.

    If the file carries an ASErtain label header and `canonicalize` is set, the
    labelled column names / direction values are mapped back to canonical.
    External tables (e.g. a DESeq2 DE table) have no label header and are read
    verbatim.
    """
    rows = _read_dicts(path)
    if not canonicalize:
        return rows
    return _canon_rows(rows, parse_comment(path))
