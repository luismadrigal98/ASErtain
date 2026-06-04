# ASErtain — design notes

ASErtain detects **cis-regulatory divergence** from **allele-specific expression
(ASE)** in F1 hybrids. It is built around a general two-lineage / F1 design with
**outbred parents**, **nested replication**, and **RNA-seq-only** data, and is
not tied to any particular organism.

## Terminology (generic)

The toolkit deliberately uses neutral labels so any cross fits:

| Generic term         | Meaning                                                              |
|----------------------|---------------------------------------------------------------------|
| **variable lineage** | the parental lineage that may have several cross parents and carries the expression difference of interest |
| **fixed lineage**    | the other parental lineage                                          |
| **F1 plant**         | one F1 individual / genome — the **biological replicate**            |
| **flower**           | one RNA sample, **nested within a plant** (technical sub-sample)     |
| **background**       | an F1 grouping, by default its mother (variable-lineage parent)      |
| **informative SNP**  | a site where, *for a given F1*, we can tell which allele is maternal (variable) vs paternal (fixed) — resolved by phasing |

You assign the real-world labels (display names) in the config; the code only
cares about the roles.

## Why F1 ASE answers "is the parental DE caused by cis divergence?"

In an F1 hybrid both alleles share one nucleus and therefore one *trans*
environment. Allelic imbalance in the F1 is thus a clean readout of **cis**
divergence, while the parental expression difference reflects **cis + trans**:

```
cis    = log2 F1 allelic ratio (variable / fixed)
total  = log2 parental expression ratio (variable / fixed)   [your DE result]
trans  = total - cis
```

The `contrast` stage classifies each gene into conserved / cis-only / trans-only
/ cis+trans / cis×trans (opposing) / compensatory.

## Pipeline stages

```
diagnose        VCF  -> informative SNPs    (phased per F1 plant from its parents)
count           BAMs -> per-SNP allele counts (mpileup; flag-driven ref-bias handling)
test            counts -> gene-level ASE    (nested flower->plant statistics)
contrast        +parental DE -> cis/trans classes
report          -> HTML summary (+ optional volcano plot)
run             -> all of the above from one config
mask-reference  informative SNPs -> N-masked reference (+ WASP SNP files) for de-biasing
```

## Four design decisions that matter

### 1. Informative SNPs are phased per F1 individual (outbred-safe)

With outbred parents there is rarely a site "fixed between species", so the
species-pool / fixed-difference approach fails. Instead, for each F1 plant we
resolve *phase* — which allele is maternal (variable) vs paternal (fixed) — from
that F1's own genotype together with its two named parents:

* **`both_hom`** — both parents homozygous for different alleles → inheritance
  fixes both contributed alleles. Needs no F1 genotype.
* **`phased`** — one parent homozygous + the F1's expressed alleles → the
  homozygous parent fixes its contribution, so the F1's other allele is assigned
  to the opposite parent. This rescues sites where the *other* parent is
  heterozygous (the common outbred case).
* otherwise (both parents heterozygous / ambiguous) → uninformative.

Because each F1 has its own parents, the variable/fixed allele identity is
tracked **per plant**. A site informative and concordant for every plant is
`shared`; otherwise it is `plant_specific` and used only for the plants it is
valid for.

### 2. RNA-only genotyping: phase from parents, never demand a heterozygous F1

Genotypes come from *expression*, so a heterozygous site under strong allelic
imbalance is easily miscalled homozygous — and for the F1 those are exactly the
strongest-ASE sites. ASErtain therefore takes phase from the **parents** (a
genetic fact) and uses the F1 genotype only to *observe* which allele it carries.
A `both_hom` site is kept even if the F1 looks homozygous, and a `phased` site is
kept as long as the F1's non-fixed allele is visible. This deliberately retains
the strong-ASE signal instead of discarding it.

Caveat surfaced in the output: parental genotypes are themselves RNA-based, so a
parent with its own allelic imbalance can be miscalled homozygous. Use decent
parent depth (`--min-parent-depth`) and a strict `--maf-threshold`, and treat
extreme gene calls together with the `fixed_allele_seen` column (below).

### 3. Reference-mapping bias is selectable by flag — and matters more here

When the reference is one of the parents (a single parental haplotype), reads
carrying that parent's allele map preferentially and the other allele can be
**lost entirely**,
manufacturing an apparent "complete ASE" toward the reference. Strategies are
flag-selectable because the right one depends on what the reads were mapped to:

| `--bias-mode` | behaviour |
|---------------|-----------|
| `none`        | null expectation = 0.5, no bookkeeping |
| `report`      | null = 0.5, but record `variable_is_ref` so systematic pull is detectable |
| `null-shift`  | per-SNP null from a balanced-control table (e.g. F1 gDNA reference-allele fraction) |
| `wasp`        | expect WASP-filtered BAMs (also corrects read-level multi-mismatch bias) |
| `nmask`       | expect BAMs realigned to an N-masked reference |

`asertain mask-reference` generates the N-masked reference (and per-chromosome
WASP SNP files) from the informative-SNP set; you re-align, then count with
`--bias-mode nmask` (or `wasp`). The `test` stage's `fixed_allele_seen` /
`n_plants_fixed_seen` columns let you distinguish a *real* complete imbalance
from an allele that was never observable due to mapping loss.

### 4. Statistics respect the nested replication, and stay honest at small n

Replication is hierarchical: background → F1 plant → flower. Flowers from one
plant share a genome, so they are **not** independent biological replicates.

* **Collapse flowers** — each plant's flowers are summed into per-SNP (variable,
  fixed) counts. The F1 plant is the unit of inference.
* **Per-plant test** — each plant is tested on its own: a beta-binomial across
  that plant's informative SNPs when it has ≥3 (the overdispersion absorbs
  SNP/mapping noise), otherwise a pooled binomial. A multi-start optimiser and a
  convergence gate keep the fit stable; non-converged fits fall back to binomial.
* **Combine by intersection–union (max-p)** — the gene's p-value is the *largest*
  per-plant p, so a call requires *every* contributing plant to be individually
  significant in the same direction. This is valid and conservative at small n
  (it needs no across-plant variance estimate — a beta-binomial *across* only 2
  plants is unidentifiable), and it bakes the cross-background consistency in
  rather than treating it as a footnote.
* **Consistency** — `consistent_backgrounds` requires ≥2 backgrounds present and
  one shared direction; a gene resolved in a single background cannot pass.
* **Honesty flags** — `phase_concordant` (per-SNP directions agree within each
  plant), `fixed_allele_seen` / `possible_ref_bias` (separating real complete ASE
  from a fixed allele lost to mapping bias under a single-parent reference).
* **Descriptive only** — the pooled binomial and the per-plant logit t-test are
  reported but never drive calls.

This replaces an earlier across-plants beta-binomial that was uncalibrated at
n=2 plants (see AUDIT.md, finding C2).

## File-format contracts

All inter-stage files are TSV with a `#` comment block and one header line; the
readers/writers live in `tables.py`. Stages can therefore be run individually or
chained.

## Pure-Python footprint

Statistics use only `numpy`/`scipy`. The only external programs invoked (via
`subprocess`, in `external.py`) are `samtools` (always) and, optionally, GATK or
WASP. There is no R dependency.
