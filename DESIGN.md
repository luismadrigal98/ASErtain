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
count           BAMs -> allele counts (per-SNP pileup, or read-backed haplotype: --counter)
test            counts -> gene-level ASE    (nested flower->plant statistics, flower-normalised)
parental-de     parent RNA BAMs -> variable-vs-fixed DE (total = cis + trans)
contrast        gene ASE + parental DE -> cis/trans classes + ASE-vs-DE sanity check
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

### Flower normalisation — equalising technical-replicate contribution

Flowers are technical replicates nested in a plant, and they differ in
sequencing depth. If we just summed raw counts, a flower with 10× the reads
would set the plant's allelic ratio almost single-handedly — the imbalance call
would rest on one library, not on the plant. Before pooling, each flower is
therefore rescaled by a **per-plant size factor**:

```
size_factor(flower) = depth(flower) / geomean(depths of the plant's flowers)
normalised count    = raw count / size_factor(flower)
```

The depth is the flower's total allele-bearing reads summed across **all** of the
plant's SNPs (a stable, library-level quantity, not a noisy per-gene one). Every
flower of a plant is thereby brought to the same effective depth, so all flowers
contribute comparably. Because both alleles are divided by the *same* factor,
each flower's own allelic ratio is untouched — only its **weight** in the pooled
(k, n) changes. The rescaled counts are rounded back to integers for the
binomial / beta-binomial likelihoods (rounding is monotone, so k ≤ n holds).

This is `--flower-norm equalize` (default). `--flower-norm none` recovers the old
raw-summing behaviour. The displayed gene `log2_ratio` / `mean_*_ratio` are
computed from the *normalised* plant-summed counts so they match what the test
saw; `*_reads` stay raw (true read totals, for provenance).

*Honesty.* Equalising **weight** means each flower's ratio is trusted equally
regardless of depth — so a deep flower is down-weighted and a shallow one
up-weighted. That is the point (a flower is an observation, not a depth quota),
but it reweights statistical confidence rather than keeping it strictly
depth-proportional. The plant stays the unit of inference, and the across-plant
max-p rule plus the cross-background consistency requirement are the real guards
against any single noisy flower. In practice flowers of one plant have similar
depth, so the reweighting is mild and mainly caps the occasional outlier-deep
flower that would otherwise set the call on its own (see the worked example in
the test suite).

### Read-backed counting removes within-gene SNP non-independence

The per-SNP pileup counter treats each informative SNP as a separate
observation, but SNPs in one gene are **not** independent: a read/fragment
spanning several SNPs is counted once at *each*, so the same molecule is
double-counted and the per-plant depth — and p-value — is inflated. The
beta-binomial absorbs SNP-to-SNP *dispersion* (via ρ) but not this *read-level*
correlation.

`--counter haplotype` fixes it at the level of the **read**. For each fragment
overlapping a gene's informative SNPs, ASErtain reads the allele it carries at
every such SNP (walking the CIGAR, so introns and indels are handled), assigns
the *whole fragment* to the variable or fixed haplotype, and counts it **once**
per gene. The result is a single (K, N) of genuinely independent reads per
gene × plant, so the per-plant test is a clean **binomial over independent
reads** — no SNP pseudo-replication. A fragment carrying *both* a variable and a
fixed allele (sequencing error, mis-phased SNP, or recombinant) is conservatively
called **ambiguous** and excluded from K/N (kept in `other_count` for QC); paired
mates share a QNAME so their votes are pooled and the fragment counted once,
which also removes mate-overlap double-counting.

The ambiguous fraction is itself a phasing-quality check, surfaced per gene as
`ambiguous_fraction` and the `low_ambiguity` flag: a gene whose reads are more
than `--max-other-fraction` (default 0.10) ambiguous is flagged and **not
called**, because a high ambiguous rate signals a mis-phased SNP or a mapping
artefact rather than a clean allelic signal. (For `--counter pileup` the same
columns report the third-allele/error read fraction.)

It needs no reference FASTA (the read's own base is used, not an mpileup match
symbol) and works with the `nmask` / `wasp` de-biasing BAMs. It requires the
SNPs to be gene-annotated (it groups by gene). `--counter pileup` (default)
keeps the per-SNP behaviour and the across-SNP beta-binomial; choose `haplotype`
when reads span multiple SNPs per gene (the usual case) for the most defensible
per-plant statistics. The `test` stage is unchanged: haplotype counts arrive as
one record per gene, so each gene×plant is a single (K, N) and the binomial path
is taken; `n_snps` still reports the real number of SNPs phased into the reads.

## Parental differential expression and the ASE sanity check

The F1 allelic ratio is **cis** only. The parental expression difference is
**cis + trans** (the *total*). `asertain parental-de` estimates that total from
the parents' own RNA-seq:

* count reads per gene per parental library (`samtools view -c` over the gene
  interval — a gene-*region* proxy for an exon-union count, deliberately
  pure-`samtools`);
* normalise for sequencing depth by **library size** (`samtools idxstats` total
  mapped reads) — robust even for a handful of candidate genes, where a
  median-of-ratios size factor would be unstable;
* **collapse flowers to their genotype** (the biological unit) before testing,
  exactly as F1 flowers collapse into their plant on the ASE side — so each
  genotype counts once regardless of its flower count (amphorellae's 6 flowers
  do not outvote k2's 4 + k3's 3) and the test is across genotypes, not across
  pseudoreplicated flowers;
* test each gene variable-vs-fixed with a Welch t-test on log2(genotype mean
  + 1), oriented variable/fixed, BH-adjusted. A valid test needs ≥2 genotypes
  per lineage; with one genotype it falls back to a flower-level test flagged
  `method = welch_flower_pseudorep` (`n_variable`/`n_fixed` report the genotype
  counts, `n_*_flowers` the flower counts).

**Honesty.** Two limits are surfaced, not hidden. (a) Gene-region counts include
intronic overlap — fine for a direction/magnitude check, not a substitute for a
featureCounts/DESeq2 table (which you can pass to `contrast` directly). (b) If a
lineage has a single genotype sampled as several flowers (as the fixed lineage
often is), those flowers are technical replicates and the p-value is
anticonservative — ASErtain warns, and you should trust the fold-change
*direction* over the p-value.

**The sanity check.** `contrast` joins ASE (cis) with parental DE (total) and,
for every ASE candidate, asks whether the allelic shift points toward the parent
the DE says is more highly expressed:

| `sanity_check` | meaning |
|----------------|---------|
| `concordant`   | ASE cis direction agrees with the parental DE direction (expected for a cis-driven difference) |
| `discordant_compensatory` | significant ASE opposing the DE — opposing cis/trans (compensatory); biologically real but flagged for inspection, not silently trusted |
| `de_not_sig` / `no_parental_data` | DE absent or non-significant — cannot check |
| `not_ase`      | gene is not an ASE candidate |

This bakes the user's request — *an ASE candidate's shift should be toward the
significantly more-expressed parent* — into a column, while still surfacing the
genuinely interesting compensatory cases rather than discarding them.

## User labels in every output

The code works internally in the canonical roles `variable` / `fixed`, but the
tables and report show your names (`variable_label` / `fixed_label` in the
config). On write, the canonical column names and `direction` values are
rewritten to the labels (`variable_count` → `kunthii_count`,
`direction = kunthii`) and the labels are stamped into the file's `#` header; on
read, the header labels are mapped back to canonical so the stages still chain.
Files with no ASErtain label header (an external DESeq2 DE table) are read
verbatim. This lives in `labels.py` / `tables.py`.

## File-format contracts

All inter-stage files are TSV with a `#` comment block and one header line; the
readers/writers live in `tables.py`. Stages can therefore be run individually or
chained.

## Pure-Python footprint

Statistics use only `numpy`/`scipy`. The only external programs invoked (via
`subprocess`, in `external.py`) are `samtools` (always) and, optionally, GATK or
WASP. There is no R dependency.
