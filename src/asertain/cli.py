"""`asertain` command-line access point.

A single entry with one subcommand per pipeline stage. Each subcommand reads/
writes the shared TSV contracts in `tables`, so stages can be run individually
or chained with `asertain run`.

    asertain diagnose  --config cfg.yaml --vcf v.vcf.gz --out runs/fls
    asertain count     --config cfg.yaml --snps runs/fls.diagnostic_snps.tsv --out runs/fls
    asertain test      --counts runs/fls.allele_counts.tsv --out runs/fls
    asertain contrast  --gene-ase runs/fls.gene_ase.tsv --parental-de de.tsv --out runs/fls
    asertain report    --gene-ase runs/fls.gene_ase.tsv --out runs/fls
    asertain run       --config cfg.yaml --vcf v.vcf.gz --out runs/fls [--parental-de de.tsv]
    asertain check     --config cfg.yaml
"""
from __future__ import annotations

import argparse
import os
import sys

from . import __version__


# ---------------------------------------------------------------------------
# Shared option helpers
# ---------------------------------------------------------------------------

def _add_diagnose_filters(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("genotype / SNP filters")
    g.add_argument("--min-parent-depth", type=int, default=8,
                   help="Min depth to call a parent genotype (default: 8)")
    g.add_argument("--maf-threshold", type=float, default=0.10,
                   help="Max minor-allele fraction to still call a parent "
                        "homozygous, tolerating sequencing noise (default: 0.10)")
    g.add_argument("--min-qual", type=float, default=30.0,
                   help="Min variant QUAL (default: 30)")
    g.add_argument("--chrom-filter", default=None,
                   help="Substring a contig name must contain (e.g. 'Chr')")


def _add_bias_opts(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("reference-bias handling")
    g.add_argument("--bias-mode", default="report",
                   choices=["none", "report", "null-shift", "wasp", "nmask"],
                   help="How to handle reference mapping bias (default: report)")
    g.add_argument("--control-table", default=None,
                   help="For --bias-mode null-shift: TSV of chrom,pos,ref_fraction "
                        "from a balanced control (e.g. F1 gDNA)")


def _add_count_opts(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("read counting")
    g.add_argument("--min-mapq", type=int, default=20,
                   help="Min mapping quality (default: 20)")
    g.add_argument("--min-baseq", type=int, default=20,
                   help="Min base quality (default: 20)")
    g.add_argument("--min-count-depth", type=int, default=10,
                   help="Min total depth per SNP×replicate to keep (default: 10)")
    g.add_argument("--samtools", default="samtools",
                   help="samtools executable (default: samtools)")


def _add_test_opts(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("statistics")
    g.add_argument("--alpha", type=float, default=0.05,
                   help="FDR threshold for ASE calls (default: 0.05)")
    g.add_argument("--min-effect-log2", type=float, default=0.0,
                   help="Min |log2 allelic ratio| to call ASE (default: 0)")


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def cmd_diagnose(args) -> int:
    from .annotation import GeneIndex
    from .config import load_config
    from .genotypes import find_diagnostic_snps
    from .tables import write_bed, write_diagnostic_snps

    cfg = load_config(args.config)
    _ensure_out_dir(args.out)
    gene_index = GeneIndex.from_file(cfg.gtf) if cfg.gtf else None
    if gene_index:
        print(f"Loaded {gene_index.n_genes} genes for annotation")

    snps, stats = find_diagnostic_snps(
        cfg, args.vcf,
        min_depth=args.min_parent_depth, maf_threshold=args.maf_threshold,
        min_qual=args.min_qual, chrom_filter=args.chrom_filter,
        gene_index=gene_index)

    write_diagnostic_snps(snps, f"{args.out}.diagnostic_snps.tsv", source_vcf=args.vcf)
    write_bed(snps, f"{args.out}.diagnostic_snps.bed")
    print(f"Variants scanned        : {stats.total}")
    print(f"Biallelic SNPs          : {stats.biallelic_snp}")
    print(f"Fixed in fixed-species  : {stats.fixed_in_fixed_species}")
    print(f"Diagnostic SNPs         : {len(snps)}")
    print(f"  shared (all F1s)      : {stats.shared}")
    print(f"  background-specific   : {stats.background_specific}")
    print(f"Wrote {args.out}.diagnostic_snps.tsv / .bed")
    return 0


def cmd_count(args) -> int:
    from .config import load_config
    from .counting import count_f1_samples
    from .tables import read_diagnostic_snps_as_objects, write_allele_counts

    cfg = load_config(args.config)
    _ensure_out_dir(args.out)
    snps = read_diagnostic_snps_as_objects(args.snps)
    print(f"Loaded {len(snps)} diagnostic SNPs from {args.snps}")
    records = count_f1_samples(
        cfg, snps, bias_mode=args.bias_mode, control_table=args.control_table,
        min_mapq=args.min_mapq, min_baseq=args.min_baseq,
        min_depth=args.min_count_depth, samtools=args.samtools)
    write_allele_counts(records, f"{args.out}.allele_counts.tsv",
                        bias_mode=args.bias_mode)
    print(f"Wrote {len(records)} observations to {args.out}.allele_counts.tsv")
    return 0


def cmd_test(args) -> int:
    from .tables import read_allele_counts, write_table
    from .testing import GENE_COLS, test_genes

    _ensure_out_dir(args.out)
    counts = read_allele_counts(args.counts)
    genes = test_genes(counts, alpha=args.alpha,
                       min_effect_log2=args.min_effect_log2)
    write_table(genes, GENE_COLS, f"{args.out}.gene_ase.tsv",
                comment="ASErtain gene-level ASE")
    n_ase = sum(1 for g in genes if g["ase_call"])
    print(f"Genes tested: {len(genes)}  |  ASE calls (q<{args.alpha}): {n_ase}")
    print(f"Wrote {args.out}.gene_ase.tsv")
    return 0


def cmd_contrast(args) -> int:
    from .contrast import CONTRAST_COLS, run_contrast
    from .tables import read_table, write_table

    _ensure_out_dir(args.out)
    genes = read_table(args.gene_ase)
    de = read_table(args.parental_de)
    contrasts = run_contrast(
        genes, de,
        de_gene_col=args.de_gene_col, de_log2_col=args.de_log2_col,
        de_padj_col=args.de_padj_col, ase_alpha=args.alpha,
        de_alpha=args.de_alpha, trans_log2_threshold=args.trans_log2_threshold)
    write_table(contrasts, CONTRAST_COLS, f"{args.out}.cis_trans.tsv",
                comment="ASErtain cis/trans contrast")
    print(f"Classified {len(contrasts)} genes -> {args.out}.cis_trans.tsv")
    return 0


def cmd_report(args) -> int:
    from .report import write_report
    _ensure_out_dir(args.out)
    path = write_report(args.gene_ase, f"{args.out}.report.html",
                        title=args.title)
    print(f"Wrote {path}")
    return 0


def cmd_run(args) -> int:
    from .pipeline import run_pipeline
    run_pipeline(
        args.config, args.vcf, args.out,
        parental_de=args.parental_de,
        bias_mode=args.bias_mode, control_table=args.control_table,
        min_parent_depth=args.min_parent_depth, maf_threshold=args.maf_threshold,
        min_qual=args.min_qual, chrom_filter=args.chrom_filter,
        min_mapq=args.min_mapq, min_baseq=args.min_baseq,
        min_count_depth=args.min_count_depth, alpha=args.alpha,
        min_effect_log2=args.min_effect_log2, samtools=args.samtools)
    return 0


def cmd_check(args) -> int:
    from .config import load_config
    from . import external

    cfg = load_config(args.config)
    print(f"Config OK: project '{cfg.project}'")
    print(f"  variable ({cfg.variable_label}) parents: "
          f"{[p.name for p in cfg.variable_parents]}")
    print(f"  fixed ({cfg.fixed_label}) parents     : "
          f"{[p.name for p in cfg.fixed_parents]}")
    print(f"  F1 backgrounds: {cfg.backgrounds}")
    print(f"  reference identity: {cfg.reference.identity}")
    print("Tool check:")
    external.check_tool(args.samtools)
    missing = [r.name for r in cfg.f1 if not os.path.exists(r.bam)]
    if missing:
        print(f"  ✗ missing BAMs for: {missing}")
        return 1
    print("  ✓ all F1 BAMs present")
    return 0


# ---------------------------------------------------------------------------
# Parser construction
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="asertain",
        description="Allele-specific expression discovery & validation for F1 hybrids.")
    p.add_argument("--version", action="version",
                   version=f"asertain {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    # diagnose
    d = sub.add_parser("diagnose", help="find diagnostic SNPs from parent genotypes")
    d.add_argument("--config", required=True)
    d.add_argument("--vcf", required=True)
    d.add_argument("--out", required=True, help="output prefix")
    _add_diagnose_filters(d)
    d.set_defaults(func=cmd_diagnose)

    # count
    c = sub.add_parser("count", help="count allele-specific reads in F1 BAMs")
    c.add_argument("--config", required=True)
    c.add_argument("--snps", required=True, help="diagnose-stage .diagnostic_snps.tsv")
    c.add_argument("--out", required=True, help="output prefix")
    _add_bias_opts(c)
    _add_count_opts(c)
    c.set_defaults(func=cmd_count)

    # test
    t = sub.add_parser("test", help="replicate-aware gene-level ASE statistics")
    t.add_argument("--counts", required=True, help="count-stage .allele_counts.tsv")
    t.add_argument("--out", required=True, help="output prefix")
    _add_test_opts(t)
    t.set_defaults(func=cmd_test)

    # contrast
    ct = sub.add_parser("contrast", help="cis/trans decomposition vs parental DE")
    ct.add_argument("--gene-ase", required=True, help="test-stage .gene_ase.tsv")
    ct.add_argument("--parental-de", required=True, help="parental DE table (TSV)")
    ct.add_argument("--out", required=True, help="output prefix")
    ct.add_argument("--de-gene-col", default="gene_id")
    ct.add_argument("--de-log2-col", default="log2FoldChange",
                    help="DE column with log2FC oriented variable/fixed")
    ct.add_argument("--de-padj-col", default="padj")
    ct.add_argument("--alpha", type=float, default=0.05, help="ASE FDR threshold")
    ct.add_argument("--de-alpha", type=float, default=0.05, help="DE FDR threshold")
    ct.add_argument("--trans-log2-threshold", type=float, default=1.0,
                    help="|trans log2| above which trans is flagged (approximate)")
    ct.set_defaults(func=cmd_contrast)

    # report
    r = sub.add_parser("report", help="HTML summary + optional plot")
    r.add_argument("--gene-ase", required=True, help="test-stage .gene_ase.tsv")
    r.add_argument("--out", required=True, help="output prefix")
    r.add_argument("--title", default="ASErtain report")
    r.set_defaults(func=cmd_report)

    # run
    ru = sub.add_parser("run", help="run the full pipeline from one config")
    ru.add_argument("--config", required=True)
    ru.add_argument("--vcf", required=True)
    ru.add_argument("--out", required=True, help="output prefix")
    ru.add_argument("--parental-de", default=None,
                    help="optional parental DE table to enable cis/trans contrast")
    _add_diagnose_filters(ru)
    _add_bias_opts(ru)
    _add_count_opts(ru)
    _add_test_opts(ru)
    ru.set_defaults(func=cmd_run)

    # check
    ch = sub.add_parser("check", help="validate config + tool availability")
    ch.add_argument("--config", required=True)
    ch.add_argument("--samtools", default="samtools")
    ch.set_defaults(func=cmd_check)

    return p


def _ensure_out_dir(prefix: str) -> None:
    d = os.path.dirname(prefix)
    if d:
        os.makedirs(d, exist_ok=True)


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
