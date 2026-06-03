# ASErtain

**Allele-specific expression (ASE) discovery and validation for F1 hybrids.**

ASErtain turns the by-eye IGV approach to ASE into a reproducible, statistically
defensible pipeline. Given a multi-sample VCF (parents + F1s) and F1 BAM files,
it finds diagnostic SNPs from *exact-parent* genotypes, counts allele-specific
reads in the F1s, tests for allelic imbalance with **replicate-aware** statistics,
and — if you supply a parental differential-expression table — decomposes each
gene's regulatory divergence into **cis** and **trans** components.

It is organism-agnostic: you label your two parental lineages and F1 groupings in
a config file, and the code works in terms of generic roles
(*variable species* / *fixed species* / *background*). See **DESIGN.md** for the
scientific rationale.

## Install

```bash
pip install -e .            # console script: `asertain`
pip install -e .[plots]     # also enable the report volcano plot (matplotlib)
```

Requires Python ≥ 3.9, `numpy`, `scipy`, `pyyaml`, and `samtools` on `PATH`.

## Quick start

1. Describe your cross in a YAML config (see `example_config.yaml`):
   which exact parent individuals belong to each species, and which F1 replicate maps
   to which parental background.

2. Run the whole pipeline:

```bash
asertain run \
    --config   my_cross.yaml \
    --vcf      variants.vcf.gz \
    --out      runs/study \
    --parental-de parental_DE.tsv \   # optional → enables cis/trans contrast
    --bias-mode report
```

…or run the stages individually:

```bash
asertain diagnose --config my_cross.yaml --vcf variants.vcf.gz --out runs/study
asertain count    --config my_cross.yaml --snps runs/study.diagnostic_snps.tsv --out runs/study --bias-mode report
asertain test     --counts runs/study.allele_counts.tsv --out runs/study
asertain contrast --gene-ase runs/study.gene_ase.tsv --parental-de parental_DE.tsv --out runs/study
asertain report   --gene-ase runs/study.gene_ase.tsv --out runs/study
asertain check    --config my_cross.yaml          # validate config + tools
```

## Pipeline stages

| Subcommand | Input | Output |
|-----------|-------|--------|
| `diagnose` | multi-sample VCF + config | `*.diagnostic_snps.tsv`, `*.bed` |
| `count`    | diagnostic SNPs + F1 BAMs | `*.allele_counts.tsv` |
| `test`     | allele counts | `*.gene_ase.tsv` |
| `contrast` | gene ASE + parental DE | `*.cis_trans.tsv` |
| `report`   | gene ASE | `*.report.html` (+ plot) |
| `run`      | config + VCF (+ DE) | all of the above |

## What makes the calls trustworthy

* **Exact-parent diagnostic SNPs** — sites must be homozygous and concordant
  across the parents of each species, not merely frequent in a pool. Sites where
  one parental background disagrees are kept as *background-specific* and used
  only where they are valid.
* **Flag-driven reference-bias handling** — `none` / `report` / `null-shift` /
  `wasp` / `nmask`, so you can match whatever reference your reads were aligned to.
* **Replicate-aware statistics** — the F1 individual is the unit of inference
  (per-replicate logit t-test + beta-binomial), with a cross-background
  consistency requirement; the anti-conservative pooled binomial is reported only
  as a descriptor.

## Status

Fully implemented: `diagnose`, `count`, `test`, plus the CLI, config, and
file-format layer. Working scaffolds to extend: `contrast` (category logic done;
a formal *trans* significance test is marked for extension), `report`, and the
WASP wrapper in `external.py`.

## Layout

```
src/asertain/
  cli.py          single access point, one subcommand per stage
  config.py       cross-design config (roles, parents, F1 backgrounds)
  vcf.py          minimal VCF reader
  genotypes.py    per-parent genotyping + diagnostic-SNP logic
  counting.py     mpileup allele counting + reference-bias modes
  stats.py        binomial / beta-binomial / BH / per-replicate logit tests
  testing.py      gene-level aggregation and ASE calls
  contrast.py     cis/trans decomposition
  annotation.py   GTF/GFF3 gene index
  tables.py       inter-stage TSV read/write contracts
  external.py     subprocess wrappers (samtools, GATK, WASP) + tool checks
  report.py       HTML/plot summary
  pipeline.py     run-everything orchestration
```
