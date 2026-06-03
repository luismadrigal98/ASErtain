# ASErtain — design notes

ASErtain detects **cis-regulatory divergence** from **allele-specific expression
(ASE)** in F1 hybrids. It is built around a general two-parent / F1 design and is
not tied to any particular organism.

## Terminology (generic)

The toolkit deliberately uses neutral labels so any cross fits:

| Generic term        | Meaning                                                              |
|---------------------|---------------------------------------------------------------------|
| **variable species**| the parental lineage that may have several cross parents and carries the expression difference of interest |
| **fixed species**   | the other parental lineage, typically a single cross parent          |
| **background**      | an F1 grouping by which variable-species parent it descends from     |
| **diagnostic SNP**  | a site fixed for different alleles between the two parents, used to assign F1 reads to a parent of origin |

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
diagnose  VCF  -> diagnostic SNPs        (per-parent genotypes, not pooled freqs)
count     BAMs -> per-SNP allele counts  (mpileup; flag-driven ref-bias handling)
test      counts -> gene-level ASE       (replicate-aware statistics)
contrast  +parental DE -> cis/trans classes
report    -> HTML summary (+ optional volcano plot)
run       -> all of the above from one config
```

## Three design decisions that matter

### 1. Diagnostic-SNP selection uses *exact-parent* genotypes

Rather than pooling a whole species and calling fixation by frequency, each named
parent plant is genotyped individually. A site is kept only when:

1. every variable-species parent is homozygous and they **agree**,
2. the fixed-species parent(s) are homozygous, and
3. the two species are fixed for **different** alleles.

Sites where the variable-species parents disagree are not discarded — they are
emitted as **`background_specific`** and used only for the F1s descending from
the parent for which they are cleanly diagnostic. This preserves a robust
**`shared`** set for combined analysis while retaining per-background power.

Genotype calling tolerates RNA-seq noise: a parent is called homozygous when the
minor-allele fraction (from `AD`) is below `--maf-threshold`, not strictly zero.

### 2. Reference-mapping bias is selectable by flag

Reads carrying the reference allele map preferentially, biasing the allelic
ratio. Because the appropriate handling depends on what the reads were mapped to
(one parent, the other, or a third reference), every strategy is a flag:

| `--bias-mode` | behaviour |
|---------------|-----------|
| `none`        | null expectation = 0.5, no bookkeeping |
| `report`      | null = 0.5, but record `variable_is_ref` so systematic pull is detectable |
| `null-shift`  | per-SNP null from a balanced-control table (e.g. F1 gDNA reference-allele fraction) |
| `wasp`        | expect WASP-filtered BAMs (remap-and-filter) |
| `nmask`       | expect BAMs aligned to an N-masked reference |

The `test` stage consumes the per-record `null_p` written by `count`, so the
choice is transparent downstream.

### 3. Statistics respect the replicate structure

The unit of biological replication is the **F1 plant**, not the SNP or the read.

* **Primary test** — collapse each gene's diagnostic SNPs to one (variable,
  fixed) count per replicate, take the per-replicate logit allelic ratio, and run
  a one-sample t-test against the null. No pseudoreplication, correct error model.
* **Secondary test** — a beta-binomial LRT across replicates that uses read depth
  while modelling overdispersion.
* **Descriptive only** — the pooled binomial across all reads is reported but
  flagged anti-conservative.
* **Consistency** — a gene is called ASE only if the effect points the same way
  across all backgrounds (e.g. both variable-species parents), which guards
  against parent-specific artefacts.

## File-format contracts

All inter-stage files are TSV with a `#` comment block and one header line; the
readers/writers live in `tables.py`. Stages can therefore be run individually or
chained.

## Pure-Python footprint

Statistics use only `numpy`/`scipy`. The only external programs invoked (via
`subprocess`, in `external.py`) are `samtools` (always) and, optionally, GATK or
WASP. There is no R dependency.
