"""End-to-end orchestration: `asertain run --config ... --vcf ... --out ...`.

Chains diagnose → count → test → (contrast) → report using a single config and
output prefix. Each stage still writes its own intermediate TSV so a run can be
inspected or resumed by hand.
"""
from __future__ import annotations

import os
from typing import Optional

from . import (counting, testing, contrast as contrast_mod,
               report as report_mod, expression as expr_mod)
from .annotation import GeneIndex
from .config import load_config
from .genotypes import find_informative_snps
from .labels import Labels
from .tables import (read_table, write_allele_counts, write_bed,
                     write_informative_snps, write_table)
from .testing import GENE_COLS
from .contrast import CONTRAST_COLS
from .expression import DE_COLS, DE_DIRECTION_COLS


def run_pipeline(config: str, vcf: str, out_prefix: str, *,
                 parental_de: Optional[str] = None,
                 compute_parental_de: bool = False,
                 bias_mode: str = "report",
                 control_table: Optional[str] = None,
                 min_parent_depth: int = 8,
                 maf_threshold: float = 0.10,
                 min_qual: float = 30.0,
                 chrom_filter: Optional[str] = None,
                 min_mapq: int = 20, min_baseq: int = 20,
                 min_count_depth: int = 10,
                 alpha: float = 0.05,
                 de_alpha: float = 0.05,
                 min_effect_log2: float = 0.0,
                 min_plants: int = 2,
                 flower_norm: str = "equalize",
                 samtools: str = "samtools",
                 verbose: bool = False) -> dict:
    cfg = load_config(config)
    labels = Labels.from_config(cfg)
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    gene_index = GeneIndex.from_file(cfg.gtf) if cfg.gtf else None

    print("[1/5] diagnose: finding informative SNPs ...")
    snps, dstats = find_informative_snps(
        cfg, vcf, min_depth=min_parent_depth, maf_threshold=maf_threshold,
        min_qual=min_qual, chrom_filter=chrom_filter, gene_index=gene_index)
    write_informative_snps(snps, f"{out_prefix}.informative_snps.tsv",
                           source_vcf=vcf, labels=labels)
    write_bed(snps, f"{out_prefix}.informative_snps.bed")
    print(f"      {len(snps)} informative SNPs "
          f"({dstats.shared} shared, {dstats.plant_specific} plant-specific)")

    print("[2/5] count: allele-specific counting in F1 flowers ...")
    counts = counting.count_flowers(
        cfg, snps, bias_mode=bias_mode, control_table=control_table,
        min_mapq=min_mapq, min_baseq=min_baseq, min_depth=min_count_depth,
        samtools=samtools)
    write_allele_counts(counts, f"{out_prefix}.allele_counts.tsv",
                        bias_mode=bias_mode, labels=labels)
    print(f"      {len(counts)} SNP×flower observations")

    print("[3/5] test: nested (flower→plant) gene-level ASE ...")
    genes = testing.test_genes(
        counts, alpha=alpha, min_effect_log2=min_effect_log2,
        min_plants=min_plants, flower_norm=flower_norm,
        ref_is_variable=(cfg.reference.identity == "variable"))
    write_table(genes, GENE_COLS, f"{out_prefix}.gene_ase.tsv",
                comment="ASErtain gene-level ASE", labels=labels)
    n_ase = sum(1 for g in genes if g["ase_call"])
    print(f"      {len(genes)} genes tested, {n_ase} ASE calls (q<{alpha})")

    if verbose:
        snp_rows = testing.snp_plant_detail(counts)
        gene_snp_rows = testing.snp_gene_summary(counts)
        plant_rows = testing.plant_gene_detail(counts, flower_norm=flower_norm)
        write_table(snp_rows, testing.SNP_DETAIL_COLS,
                    f"{out_prefix}.snp_gene_counts.tsv",
                    comment="ASErtain per gene×SNP×plant allele counts (flowers summed)",
                    labels=labels)
        write_table(gene_snp_rows, testing.SNP_GENE_COLS,
                    f"{out_prefix}.gene_snp_counts.tsv",
                    comment="ASErtain per gene×SNP allele counts (plants+flowers collapsed)",
                    labels=labels)
        write_table(plant_rows, testing.PLANT_DETAIL_COLS,
                    f"{out_prefix}.plant_gene_stats.tsv",
                    comment="ASErtain per gene×plant test inputs/outputs (feeds max-p)",
                    labels=labels)
        print(f"      [verbose] {out_prefix}.snp_gene_counts.tsv ({len(snp_rows)} rows)")
        print(f"      [verbose] {out_prefix}.gene_snp_counts.tsv ({len(gene_snp_rows)} rows)")
        print(f"      [verbose] {out_prefix}.plant_gene_stats.tsv ({len(plant_rows)} rows)")

    # Parental DE: an external table wins; otherwise compute it from the parents'
    # RNA BAMs when asked (and possible). Enables the cis/trans contrast and the
    # ASE-direction sanity check.
    de_path = parental_de
    if not de_path and compute_parental_de:
        if gene_index is None:
            print("      parental-de: skipped (config has no gtf)")
        elif not cfg.has_parental_expression():
            print("      parental-de: skipped (parents have no `flowers:` BAMs)")
        else:
            print("[3b] parental-de: differential expression of parental lines ...")
            # Only the ASE-tested (candidate) genes need DE for the sanity check;
            # this keeps per-gene samtools counting fast.
            cand = {g["gene_id"] for g in genes}
            de_rows = expr_mod.run_parental_de(
                cfg, gene_index, gene_ids=cand, min_mapq=min_mapq, samtools=samtools)
            de_path = f"{out_prefix}.parental_de.tsv"
            write_table(de_rows, DE_COLS, de_path,
                        comment="ASErtain parental differential expression (variable/fixed)",
                        labels=labels, direction_cols=DE_DIRECTION_COLS)
            n_de = sum(1 for r in de_rows if isinstance(r.get("padj"), float)
                       and r["padj"] < de_alpha)
            print(f"      {len(de_rows)} genes, {n_de} DE at padj<{de_alpha}")

    if de_path:
        print("[4/5] contrast: cis/trans decomposition + ASE sanity check ...")
        de = read_table(de_path)
        contrasts = contrast_mod.run_contrast(genes, de, ase_alpha=alpha,
                                              de_alpha=de_alpha)
        write_table(contrasts, CONTRAST_COLS, f"{out_prefix}.cis_trans.tsv",
                    comment="ASErtain cis/trans contrast", labels=labels)
        n_disc = sum(1 for c in contrasts
                     if c.get("sanity_check") == "discordant_compensatory")
        print(f"      {len(contrasts)} genes classified "
              f"({n_disc} ASE candidates discordant with DE — inspect)")
    else:
        print("[4/5] contrast: skipped (no parental DE table or BAMs)")

    print("[5/5] report: writing HTML summary ...")
    report_mod.write_report(f"{out_prefix}.gene_ase.tsv",
                            f"{out_prefix}.report.html",
                            title=f"ASErtain — {cfg.project}", labels=labels)
    print(f"      {out_prefix}.report.html")
    return {"n_snps": len(snps), "n_counts": len(counts),
            "n_genes": len(genes), "n_ase": n_ase}
