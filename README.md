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
asertain diagnose    --config my_cross.yaml --vcf variants.vcf.gz --out runs/study
asertain count       --config my_cross.yaml --snps runs/study.informative_snps.tsv --out runs/study --bias-mode report
asertain test        --counts runs/study.allele_counts.tsv --out runs/study   # --flower-norm equalize (default)
asertain parental-de --config my_cross.yaml --out runs/study                  # variable-vs-fixed DE from parent BAMs
asertain contrast    --gene-ase runs/study.gene_ase.tsv --parental-de runs/study.parental_de.tsv --out runs/study
asertain report      --gene-ase runs/study.gene_ase.tsv --out runs/study
asertain check       --config my_cross.yaml          # validate config + tools
```

`asertain run` can build the parental DE itself with `--compute-parental-de`
(needs `flowers:` under the parents in the config) and then use it for the
cis/trans contrast and the ASE-direction sanity check — or pass your own
DESeq2/edgeR table with `--parental-de`.

## Pipeline stages

| Subcommand | Input | Output |
|-----------|-------|--------|
| `diagnose` | multi-sample VCF + config | `*.informative_snps.tsv`, `*.bed` |
| `count`    | informative SNPs + F1 flower BAMs | `*.allele_counts.tsv` (`--counter pileup`\|`haplotype`) |
| `test`     | allele counts | `*.gene_ase.tsv` |
| `parental-de` | config (parent RNA BAMs) + GTF | `*.parental_de.tsv` (variable vs fixed) |
| `contrast` | gene ASE + parental DE | `*.cis_trans.tsv` (+ ASE-vs-DE sanity check) |
| `report`   | gene ASE | `*.report.html` (+ plot) |
| `run`      | config + VCF (+ DE) | all of the above |
| `mask-reference` | informative SNPs + reference | N-masked FASTA (+ WASP SNP files) |

**Labels in every output.** You set `variable_label` / `fixed_label` in the
config (e.g. `kunthii` / `amphorellae`); the canonical `variable`/`fixed` column
names and `direction` values are rewritten to your labels in every table and in
the report (so `variable_count` → `kunthii_count`, `direction = kunthii`). The
labels are stamped into each file's header and mapped back to canonical on read,
so the stages still chain.

Designed for **outbred parents**, **nested replication** (RNA samples within
individuals), and **RNA-seq-only** data.

### Auditing intermediate results

Add `--verbose` (to `test` or `run`) to write the full evidence trail so every
gene call can be traced and explained:

| File | Granularity |
|------|-------------|
| `*.allele_counts.tsv`     | per **flower × SNP** (raw counts, `<var>_is_ref`, `tier`) |
| `*.gene_snp_counts.tsv`   | per **gene × SNP** (plants + flowers collapsed, with a per-plant split) |
| `*.snp_gene_counts.tsv`   | per **gene × SNP × plant** (flowers summed, per-SNP ratio) |
| `*.plant_gene_stats.tsv`  | per **gene × plant** (K, N, n_snps, ρ, method, p) → fed to max-p |
| `*.gene_ase.tsv`          | per **gene** (the call) |

Reading them bottom-up shows exactly how raw reads become a call: flower counts →
per gene×SNP → per gene×SNP×plant → per-plant test → `max-p` across plants → gene.
The per-plant K/N reflect the **flower normalisation** (see below): each flower
is rescaled so a deeply sequenced one cannot dominate its plant's ratio.

## What makes the calls trustworthy

* **Phased informative SNPs** — for each F1 individual, the maternal (variable)
  and paternal (fixed) allele are resolved from that F1's genotype plus its two
  named parents, so the method works even when parents are heterozygous/outbred.
  Phase is taken from the parents (genetic fact), not from F1 expression, so
  strong-ASE sites are retained rather than miscalled away.
* **Read-backed counting (`--counter haplotype`)** — SNPs in one gene aren't
  independent (a read spanning several is counted at each), which inflates the
  per-plant depth. The haplotype counter assigns each *fragment* to a parental
  haplotype across all the SNPs it covers and counts it **once** per gene, giving
  one independent (K, N) per gene×plant and a clean binomial — no SNP
  double-counting. Fragments carrying both alleles are flagged ambiguous; needs
  no reference FASTA and works with `nmask`/`wasp` BAMs.
* **Flag-driven reference-bias handling** — `none` / `report` / `null-shift` /
  `wasp` / `nmask`, plus `mask-reference` to build the N-masked reference / WASP
  inputs. Important when the reference is one parent (the other allele can be lost).
* **Nested, replicate-aware statistics** — flowers collapse into their plant; the
  plant (individual) is the unit of inference (beta-binomial across plants +
  per-plant logit t-test), with a cross-background consistency requirement. The
  anti-conservative pooled binomial is reported only as a descriptor, and a
  `fixed_allele_seen` column separates real complete-ASE from mapping dropout.
* **Flower normalisation** — flowers (technical replicates) differ in depth, so
  before pooling each is rescaled by a per-plant size factor (`--flower-norm
  equalize`, default) and a deeply sequenced flower can no longer dominate its
  plant's allelic ratio. Each flower's own ratio is preserved (both alleles
  scale together); `--flower-norm none` recovers raw summing.
* **Parental-DE sanity check** — `asertain parental-de` computes variable-vs-fixed
  expression from the parents' RNA BAMs (library-size-normalised, pure-Python),
  and `contrast` checks each ASE candidate against it: a real cis shift should
  point toward the more-expressed parent (`sanity_check = concordant`); an ASE
  call opposing the DE (`discordant_compensatory`) is flagged as opposing
  cis/trans, not silently trusted.

## Two LD-robust ways to aggregate a gene's SNPs

SNPs in one gene are in linkage disequilibrium, so they must not be pooled
naively (a read spanning several is counted at each, inflating depth and
shrinking p-values). ASErtain offers two principled routes, selectable with
`asertain test --gene-aggregation`:

* **`plant`** (default) — pool each plant's SNPs, test per plant, combine
  plants. Fed by `--counter haplotype` this gives one *independent* (K, N) per
  gene (each fragment counted once), the most powerful option.
* **`maxsnp`** — a **plain binomial per SNP**, take the **strongest-signal SNP**
  per gene, and require the gene's SNPs to **agree in direction** (all point to
  the same parent). It never pools correlated SNPs, so LD cannot inflate it; the
  best-of-*m*-SNPs selection is corrected with `--within-gene-correction`
  (`sidak`|`bonferroni`|`none`). Extra output columns: `agg_method`, `top_snp`,
  `top_snp_p`, `n_snps_same_dir`, `snp_concordant`.

On simulated data (see `examples/benchmark_approaches.py` /
`benchmark_strata.py`) the two recover the **same ASE genes to a high degree**
(Jaccard ~85–90%; ~98% of `maxsnp` calls are also `plant`/haplotype calls) at
~0% empirical FDR — so `maxsnp` is a strong LD-robust cross-check, while the
haplotype path stays the primary because it is more sensitive on
moderate-effect genes.

**New to the pipeline? Read `ASErtain_pipeline_walkthrough.ipynb`** — an
executable, plain-language tour of every stage with a tiny worked example
(genes, SNPs, and the actual input/intermediate/output files), runnable with no
BAM/VCF/`samtools` (it simulates the counts).

## Status

Fully implemented and tested (synthetic + live BAMs): `diagnose`, `count`
(per-SNP pileup and read-backed haplotype counters), `test` (with flower
normalisation), `parental-de`, `mask-reference`, plus the CLI, config (new +
legacy schema), the label-aware file-format layer, and the ASE-vs-DE sanity
check in `contrast`. Working scaffolds to extend: `contrast`'s
formal *trans* significance test (the category logic and DE-concordance check
are done), `report`, and the WASP remap chain in `external.py` (SNP-file
generation done; the aligner step needs your alignment command).

The parental-DE stage is a lightweight, pure-`samtools` gene-region count +
library-size-normalised Welch test, intended for **candidate-gene panels** and
as a **direction sanity check**. For a genome-wide, publication-grade DE table,
run featureCounts/HTSeq + DESeq2/edgeR and feed it to `contrast --parental-de`.

## Layout

```
src/asertain/
  cli.py          single access point, one subcommand per stage
  config.py       cross-design config (roles, parents, F1 backgrounds)
  vcf.py          minimal VCF reader
  genotypes.py    per-parent genotyping + diagnostic-SNP logic
  counting.py     mpileup per-SNP allele counting + reference-bias modes
  haplotype.py    read-backed haplotype counting (--counter haplotype)
  stats.py        binomial / beta-binomial / BH / per-replicate logit tests
  testing.py      gene-level aggregation, flower normalisation, ASE calls
  expression.py   parental differential expression (variable vs fixed)
  contrast.py     cis/trans decomposition + ASE-vs-DE sanity check
  labels.py       user labels (variable/fixed → display names) for all outputs
  annotation.py   GTF/GFF3 gene index
  tables.py       inter-stage TSV read/write contracts
  external.py     subprocess wrappers (samtools, GATK, WASP) + tool checks
  report.py       HTML/plot summary
  pipeline.py     run-everything orchestration
```
