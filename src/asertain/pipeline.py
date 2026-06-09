"""End-to-end orchestration: `asertain run --config ... --vcf ... --out ...`.

Chains diagnose → count → test → (contrast) → report using a single config and
output prefix. Each stage still writes its own intermediate TSV so a run can be
inspected or resumed by hand.
"""
from __future__ import annotations

import os
from typing import Optional

from . import counting, testing, contrast as contrast_mod, report as report_mod
from .annotation import GeneIndex
from .config import load_config
from .genotypes import find_informative_snps
from .tables import (read_table, write_allele_counts, write_bed,
                     write_informative_snps, write_table)
from .testing import GENE_COLS
from .contrast import CONTRAST_COLS


def run_pipeline(config: str, vcf: str, out_prefix: str, *,
                 parental_de: Optional[str] = None,
                 bias_mode: str = "report",
                 control_table: Optional[str] = None,
                 min_parent_depth: int = 8,
                 maf_threshold: float = 0.10,
                 min_qual: float = 30.0,
                 chrom_filter: Optional[str] = None,
                 min_mapq: int = 20, min_baseq: int = 20,
                 min_count_depth: int = 10,
                 alpha: float = 0.05,
                 min_effect_log2: float = 0.0,
                 min_plants: int = 2,
                 samtools: str = "samtools",
                 verbose: bool = False) -> dict:
    cfg = load_config(config)
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    gene_index = GeneIndex.from_file(cfg.gtf) if cfg.gtf else None

    print("[1/5] diagnose: finding informative SNPs ...")
    snps, dstats = find_informative_snps(
        cfg, vcf, min_depth=min_parent_depth, maf_threshold=maf_threshold,
        min_qual=min_qual, chrom_filter=chrom_filter, gene_index=gene_index)
    write_informative_snps(snps, f"{out_prefix}.informative_snps.tsv", source_vcf=vcf)
    write_bed(snps, f"{out_prefix}.informative_snps.bed")
    print(f"      {len(snps)} informative SNPs "
          f"({dstats.shared} shared, {dstats.plant_specific} plant-specific)")

    print("[2/5] count: allele-specific counting in F1 flowers ...")
    counts = counting.count_flowers(
        cfg, snps, bias_mode=bias_mode, control_table=control_table,
        min_mapq=min_mapq, min_baseq=min_baseq, min_depth=min_count_depth,
        samtools=samtools)
    write_allele_counts(counts, f"{out_prefix}.allele_counts.tsv", bias_mode=bias_mode)
    print(f"      {len(counts)} SNP×flower observations")

    print("[3/5] test: nested (flower→plant) gene-level ASE ...")
    genes = testing.test_genes(
        counts, alpha=alpha, min_effect_log2=min_effect_log2,
        min_plants=min_plants,
        ref_is_variable=(cfg.reference.identity == "variable"))
    write_table(genes, GENE_COLS, f"{out_prefix}.gene_ase.tsv",
                comment="ASErtain gene-level ASE")
    n_ase = sum(1 for g in genes if g["ase_call"])
    print(f"      {len(genes)} genes tested, {n_ase} ASE calls (q<{alpha})")

    if verbose:
        snp_rows = testing.snp_plant_detail(counts)
        plant_rows = testing.plant_gene_detail(counts)
        write_table(snp_rows, testing.SNP_DETAIL_COLS,
                    f"{out_prefix}.snp_gene_counts.tsv",
                    comment="ASErtain per gene×SNP×plant allele counts (flowers summed)")
        write_table(plant_rows, testing.PLANT_DETAIL_COLS,
                    f"{out_prefix}.plant_gene_stats.tsv",
                    comment="ASErtain per gene×plant test inputs/outputs (feeds max-p)")
        print(f"      [verbose] {out_prefix}.snp_gene_counts.tsv ({len(snp_rows)} rows)")
        print(f"      [verbose] {out_prefix}.plant_gene_stats.tsv ({len(plant_rows)} rows)")

    if parental_de:
        print("[4/5] contrast: cis/trans decomposition ...")
        de = read_table(parental_de)
        contrasts = contrast_mod.run_contrast(genes, de, ase_alpha=alpha)
        write_table(contrasts, CONTRAST_COLS, f"{out_prefix}.cis_trans.tsv",
                    comment="ASErtain cis/trans contrast")
        print(f"      {len(contrasts)} genes classified")
    else:
        print("[4/5] contrast: skipped (no --parental-de)")

    print("[5/5] report: writing HTML summary ...")
    report_mod.write_report(f"{out_prefix}.gene_ase.tsv",
                            f"{out_prefix}.report.html",
                            title=f"ASErtain — {cfg.project}")
    print(f"      {out_prefix}.report.html")
    return {"n_snps": len(snps), "n_counts": len(counts),
            "n_genes": len(genes), "n_ase": n_ase}
