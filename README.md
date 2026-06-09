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
| `diagnose` | multi-sample VCF + config | `*.informative_snps.tsv`, `*.bed` |
| `count`    | informative SNPs + F1 flower BAMs | `*.allele_counts.tsv` |
| `test`     | allele counts | `*.gene_ase.tsv` |
| `contrast` | gene ASE + parental DE | `*.cis_trans.tsv` |
| `report`   | gene ASE | `*.report.html` (+ plot) |
| `run`      | config + VCF (+ DE) | all of the above |
| `mask-reference` | informative SNPs + reference | N-masked FASTA (+ WASP SNP files) |

Designed for **outbred parents**, **nested replication** (RNA samples within
individuals), and **RNA-seq-only** data.

### Auditing intermediate results

Add `--verbose` (to `test` or `run`) to write the full evidence trail so every
gene call can be traced and explained:

| File | Granularity |
|------|-------------|
| `*.allele_counts.tsv`     | per **flower × SNP** (raw counts, `variable_is_ref`, `tier`) |
| `*.snp_gene_counts.tsv`   | per **gene × SNP × plant** (flowers summed, per-SNP ratio) |
| `*.plant_gene_stats.tsv`  | per **gene × plant** (K, N, n_snps, ρ, method, p) → fed to max-p |
| `*.gene_ase.tsv`          | per **gene** (the call) |

Reading them bottom-up shows exactly how raw reads become a call: flower counts →
per-SNP×plant → per-plant test → `max-p` across plants → gene.

## What makes the calls trustworthy

* **Phased informative SNPs** — for each F1 individual, the maternal (variable)
  and paternal (fixed) allele are resolved from that F1's genotype plus its two
  named parents, so the method works even when parents are heterozygous/outbred.
  Phase is taken from the parents (genetic fact), not from F1 expression, so
  strong-ASE sites are retained rather than miscalled away.
* **Flag-driven reference-bias handling** — `none` / `report` / `null-shift` /
  `wasp` / `nmask`, plus `mask-reference` to build the N-masked reference / WASP
  inputs. Important when the reference is one parent (the other allele can be lost).
* **Nested, replicate-aware statistics** — flowers collapse into their plant; the
  plant (individual) is the unit of inference (beta-binomial across plants +
  per-plant logit t-test), with a cross-background consistency requirement. The
  anti-conservative pooled binomial is reported only as a descriptor, and a
  `fixed_allele_seen` column separates real complete-ASE from mapping dropout.

## Status

Fully implemented and tested (synthetic + live BAMs): `diagnose`, `count`,
`test`, `mask-reference`, plus the CLI, config (new + legacy schema), and
file-format layer. Working scaffolds to extend: `contrast` (category logic done;
a formal *trans* significance test is marked for extension), `report`, and the
WASP remap chain in `external.py` (SNP-file generation done; the aligner step
needs your alignment command).

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
