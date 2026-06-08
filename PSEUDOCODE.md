# ASErtain — pipeline pseudocode

Detailed, decision-by-decision pseudocode of the whole workflow, from raw BAMs to
ASE calls. Notation: `DECISION:` marks a branch point and why it exists; `WHY:`
gives the scientific rationale. Reflects the post-audit implementation.

---

## STAGE 0 — Preparation & variant calling  (examples/run_pkunthii_ase.sh)

```
INPUT : coordinate-sorted RNA-seq BAMs (one per flower), reference FASTA (= one
        parent's genome), Helixer GFF, candidate-loci list
OUTPUT: a multi-sample VCF (GT:AD:DP) over the candidate genes

0.1  index reference (faidx) and every flower BAM

0.2  MERGE flowers into one BAM per "genotyping unit":
        k2, k3            <- pool that kunthii parent's corolla flowers
        amphorellae       <- pool all amphorellae corolla flowers
        F1_k2a1, F1_k3a1  <- pool that F1 genome's flowers
     WHY: pooling deepens coverage so outbred parents are genotyped accurately.
     DECISION: F1 flowers may be pooled ONLY because they are one plant (one
               genome). Two distinct F1 plants of the same cross must NOT be
               pooled (siblings can differ) — this is why the config asks.

0.3  build candidate-gene BED:
        gene_id := text after "gene:" in each anthocyanin_loci line
        look gene_id up in the GFF "gene" features -> (chrom, start-1, end)

0.4  joint variant calling, RESTRICTED to the candidate BED:
        bcftools mpileup -a AD,DP --min-MQ 20 --min-BQ 20 --max-depth 100000
          DECISION (min-MQ 20): keep STAR-unique reads (MAPQ 255), drop
                                multimappers (MAPQ <= 3).
        | bcftools call -mv            # -v: keep only variant sites
        | bcftools norm
        | bcftools view -v snps -m2 -M2 # biallelic SNPs only
          WHY: phasing logic below assumes 2 alleles.

0.5  REHEADER the VCF samples to the clean names the config expects
     (k2 k3 amphorellae F1_k2a1 F1_k3a1), in bamlist order.
     WHY: bcftools names samples by the bamlist PATH (work/k2.bam); without this
          every genotype would read as "missing" and the run returns nothing.
```

---

## STAGE 1 — diagnose: informative-SNP discovery  (genotypes.py)

```
INPUT : VCF, CrossConfig, optional GeneIndex
OUTPUT: list of InformativeSNP (per-plant variable/fixed allele + tier)

1.0  GUARD: read VCF sample names.
     DECISION: if NONE of the config's vcf_sample names are in the VCF
               -> ABORT with both name lists (catches the reheader bug early).
               if SOME are missing -> WARN (those genotypes treated as missing).

1.1  FOR each VCF record (streamed):
        DECISION: skip if not biallelic SNP            (need exactly 2 alleles)
        DECISION: skip if QUAL < min_qual              (default 30)
        DECISION: skip if chrom_filter set and not matched

        # --- genotype every parent and the F1, once ---
        FOR each sample s in (all parents + each F1's vcf_sample):
            state[s] := CALL_STATE(s)        # see 1.2

        # --- resolve each F1 plant independently ---
        per_plant := {}
        FOR each F1 plant P:
            vparent := the parent of P whose lineage == "variable"   # by LINEAGE,
            fparent := the parent of P whose lineage == "fixed"      # not by sex
            allele  := INFORMATIVE_FOR_PLANT(state[vparent],
                                             state[fparent],
                                             state[P.vcf_sample])     # see 1.3
            IF allele is not None:
                per_plant[P] := allele

        DECISION: if per_plant is empty -> this SNP is uninformative, skip.

        # --- classify shared vs plant-specific ---
        DECISION:
          IF every F1 plant is informative AND all agree on the variable
             nucleotide AND all agree on the fixed nucleotide:
                classification := "shared"        (usable for all plants)
          ELSE:
                classification := "plant_specific"(usable only for listed plants)

        IF GeneIndex provided: annotate (chrom,pos) -> gene_id, location
        emit InformativeSNP
```

### 1.2  CALL_STATE(sample)  — genotype with RNA-noise tolerance

```
INPUT : the sample's GT / AD / DP at this site
OUTPUT: one of HOM_REF, HOM_ALT, HET, MISSING

DECISION order (prefer allelic depth, the most informative):
  IF AD present (ref_d, alt_d):
        total := ref_d + alt_d
        IF total < min_depth                    -> MISSING   (too little evidence)
        alt_frac := alt_d / total
        IF alt_frac        <= maf_threshold     -> HOM_REF   (default maf 0.10)
        IF (1 - alt_frac)  <= maf_threshold     -> HOM_ALT
        ELSE                                    -> HET
        WHY (maf, not 0): RNA reads carry mapping/sequencing noise; a true
             homozygote shows a few contaminating reads.
  ELSE IF DP present AND DP < min_depth          -> MISSING
  ELSE IF GT present: hom -> HOM_REF/HOM_ALT, het -> HET
  ELSE                                           -> MISSING
```

### 1.3  INFORMATIVE_FOR_PLANT(variable, fixed, f1)  — the phasing core

```
INPUT : genotype STATES of the variable-lineage parent, the fixed-lineage
        parent, and the F1 (each is an allele-index set {0}/{1}/{0,1} or None)
OUTPUT: PlantAllele(variable_nuc, fixed_nuc, tier)  or  None

KEY PRINCIPLE: take phase from the PARENTS (a genetic fact). Use the F1 genotype
only to OBSERVE which allele it carries — never require the F1 to look
heterozygous (RNA ASE can make a true het look homozygous; those are the
strongest-signal sites and must be kept).

DECISION 1 — both parents homozygous for different alleles ("both_hom"):
    IF variable in {hom} AND fixed in {hom} AND variable != fixed:
        RETURN variable=variable_allele, fixed=fixed_allele, tier="both_hom"
        WHY: inheritance fixes both contributed alleles; no F1 genotype needed,
             and an F1 that merely LOOKS homozygous (ASE) is NOT vetoed.
             >>> This is the tier that catches the kunthii(ref)/amphorellae(alt)
                 diagnostic SNPs even when the F1 is called 0/0 due to ASE. <<<

DECISION 2 — fixed-lineage parent homozygous ("phased" via F1):
    ELSE IF fixed in {hom}:
        fx := the fixed parent's allele
        IF f1 is None                          -> None   (cannot observe F1 allele)
        others := f1 - {fx}
        IF |others| != 1                       -> None   (F1 shows only fx -> no contrast)
        w := the single non-fx allele the F1 carries
        IF variable known AND w not in variable -> None  (inconsistent inheritance)
        RETURN variable=w, fixed=fx, tier="phased"
        WHY: the homozygous fixed parent pins its contribution; the F1's other
             allele must be the variable one — rescues HETEROZYGOUS variable
             parents, and an F1 expressing only the variable allele (extreme ASE)
             still resolves.

DECISION 3 — variable-lineage parent homozygous ("phased", symmetric):
    ELSE IF variable in {hom}:
        va := the variable parent's allele
        IF f1 is None                          -> None
        others := f1 - {va}
        IF |others| != 1                       -> None
        w := the single non-va allele the F1 carries
        IF fixed known AND w not in fixed       -> None
        RETURN variable=va, fixed=w, tier="phased"

DECISION 4 — otherwise (both parents heterozygous / ambiguous):
        RETURN None
        WHY: phase cannot be resolved from genotypes alone. (Read-backed phasing
             across the transcript is the principled future rescue.)
```

---

## STAGE 2 — count: allele-specific read counting  (counting.py)

```
INPUT : InformativeSNPs, CrossConfig, bias_mode, reference FASTA
OUTPUT: one row per (flower x informative SNP): variable_count, fixed_count, ...

2.0  GUARD: if reference FASTA is missing -> ABORT.
     WHY: mpileup '.'/',' (match-to-reference) symbols can only be resolved to a
          base with the reference; without it the reference allele is miscounted.

2.1  index SNPs by plant: snps_by_plant[P] = SNPs whose per_plant contains P

2.2  FOR each F1 plant P:
        usable := snps_by_plant[P]
        FOR each flower F of P:
            ensure F.bam is indexed
            FOR each SNP in usable:
                (var_nuc, fix_nuc) := SNP.per_plant[P]
                counts := COUNT_ALLELES(F.bam, SNP.pos, var_nuc, fix_nuc)  # 2.3
                DECISION: if (variable_count + fixed_count) < min_count_depth
                          -> skip this (flower, SNP)   (too few allele-bearing reads)
                          WHY: filter on reads that actually carry an allele, NOT
                               raw pileup depth (which includes intron skips).
                null_p := NULL_EXPECTATION(SNP, P, bias_mode)              # 2.4
                emit row {flower, plant, background, snp_id, var_nuc, fix_nuc,
                          variable_is_ref, tier, variable_count, fixed_count,
                          other_count, total_depth=usable_calls, null_p, gene_id}
```

### 2.3  COUNT_ALLELES(bam, pos, var_nuc, fix_nuc)

```
run: samtools mpileup -f reference -q min_mapq -Q min_baseq  at this 1bp region
parse the pileup base string into a list of base calls:
    '.'/','              -> the reference base
    'A/C/G/T'            -> that base
    '^x'                 -> read start, skip the mapping-quality char
    '$'                  -> read end
    '*'                  -> deletion placeholder (skipped)
    '<' / '>'            -> reference SKIP from a spliced RNA read (skipped)
    '+N.../-N...'        -> indel run (skipped)
variable_count := count of var_nuc ; fixed_count := count of fix_nuc
total_depth    := number of usable base calls (NOT raw mpileup depth)
```

### 2.4  NULL_EXPECTATION(SNP, plant, bias_mode)  — reference-bias handling

```
DECISION on --bias-mode:
  "none"  | "report" -> null_p := 0.5
        ("report" additionally records variable_is_ref so systematic pull toward
         the reference allele is visible downstream.)
  "null-shift"       -> look up a balanced-control table at (chrom,pos):
                          null_p := control_ref_fraction        if variable is ref
                                    1 - control_ref_fraction     otherwise
  "nmask" | "wasp"   -> null_p := 0.5
        (bias removed at the ALIGNMENT step; counting just proceeds. The N-masked
         reference / WASP SNP files are produced by `asertain mask-reference`.)
```

---

## STAGE 3 — test: nested gene-level ASE  (testing.py + stats.py)

```
INPUT : the allele-count rows, alpha, min_plants, ref_is_variable
OUTPUT: one row per gene with the ASE call

3.1  group count rows by gene_id (DECISION: drop "intergenic").

3.2  FOR each gene:
        null_p := allele-depth-weighted mean of the rows' null_p

        # --- collapse flowers -> plant, keep SNP structure ---
        FOR each plant P contributing to this gene:
            FOR each SNP:
                sum variable_count and (variable+fixed) ACROSS P's flowers
            -> P has a list of per-SNP (k, n) pairs

        # --- per-plant test (each plant on its own) ---  see 3.3
        FOR each plant P:
            plant_result[P] := PLANT_TEST(P's per-SNP pairs, null_p)
        DECISION: if no plant has a result -> skip gene.

        # --- combine plants by INTERSECTION-UNION (max-p) ---
        p_gene := MAX over plants of plant_result[P].p
        WHY: a gene is ASE only if EVERY contributing plant is individually
             significant. Max-p is valid & conservative at small n and needs no
             across-plant variance estimate (a beta-binomial across 2 plants is
             unidentifiable — audit finding C2).

        # --- direction & consistency ---
        per_plant_dir[P]   := sign of (P's pooled ratio  -  null_p)
        n_backgrounds      := number of distinct backgrounds among contributing plants
        DECISION consistent_backgrounds :=
                 (n_backgrounds >= 2) AND (all per_plant_dir agree, non-balanced)
                 WHY: a single background cannot be "consistent" (audit M4).

        # --- honesty flags ---
        phase_concordant := within EACH plant, the gene's per-SNP directions agree
                            (a mis-phased SNP would otherwise contaminate the sum)
        fixed_allele_seen := (total fixed reads > 0)
        possible_ref_bias := ref_is_variable AND NOT fixed_allele_seen
                 WHY: under a single-parent reference, a fixed allele lost to
                      mapping looks like complete ASE — flag, don't silently call.

        record p_primary=p_gene, log2_ratio (pooled), per-plant effects, flags

3.3  PLANT_TEST(per-SNP (k,n) pairs, null_p):
        DECISION:
          IF (#SNPs >= 3):
              fit a BETA-BINOMIAL across the SNPs (LRT mu vs null_p),
              multi-start optimiser.
              IF it CONVERGED -> use its p, log2
              ELSE            -> fall through to binomial
          pooled binomial:  p := binomial_test( sum k, sum n, null_p )
        WHY: >=3 SNPs lets the beta-binomial absorb SNP-to-SNP / mapping
             overdispersion (more conservative); too few -> binomial.

3.4  multiple testing & final call:
        q_value := Benjamini-Hochberg over all genes' p_primary
        DECISION ase_call := TRUE iff ALL of:
            q_value < alpha
            |log2_ratio| >= min_effect_log2
            consistent_backgrounds == TRUE      (>=2 backgrounds, same direction)
            n_plants >= min_plants              (default 2)
            phase_concordant == TRUE
        (possible_ref_bias is reported, NOT auto-vetoed: complete ASE can be real.)
```

---

## STAGE 4 — contrast: cis / trans decomposition  (contrast.py)

```
INPUT : gene-level ASE table (cis), a parental DE table (total), thresholds
OUTPUT: per-gene regulatory category

FOR each gene present in BOTH tables:
    cis    := F1 allelic log2 ratio (variable / fixed)          [from STAGE 3]
    total  := parental expression log2 fold change              [from your DE]
    trans  := total - cis
    cis_sig      := (gene's ASE q < alpha AND it was an ase_call)
    parental_sig := (DE padj < de_alpha)
    trans_sig    := |trans| >= trans_log2_threshold   # PROVISIONAL (see note)

    DECISION (McManus-style category):
       not parental_sig and not cis_sig and not trans_sig -> "conserved"
       cis_sig and not trans_sig                          -> "cis_only"
       trans_sig and not cis_sig                          -> "trans_only"
       cis_sig and trans_sig:
            same sign(cis), sign(trans):
                not parental_sig -> "compensatory"   (cis & trans cancel)
                else             -> "cis_plus_trans"  (same direction)
            opposite sign        -> "cis_x_trans_opposing"
       else                                               -> "ambiguous"

NOTE: trans_sig is currently a magnitude threshold, not a formal test, so the
      compensatory / cis_x_trans categories are PROVISIONAL. The principled
      version is a parent-vs-F1 interaction test on the allele counts.
```

---

## STAGE 5 — report  (report.py)

```
INPUT : gene-level ASE table
OUTPUT: self-contained HTML (+ volcano plot if matplotlib present)
  - headline counts (genes tested, ASE genes, variable- vs fixed-biased)
  - per-gene table (log2 ratio, q, direction, consistency, fixed_allele_seen, call)
  - DECISION: if matplotlib missing -> render table only, no plot (no failure)
```

---

## Orchestration — `asertain run`  (pipeline.py)

```
load CrossConfig
ref_is_variable := (config.reference.identity == "variable")
STAGE 1 diagnose  -> *.informative_snps.tsv / .bed
STAGE 2 count     -> *.allele_counts.tsv
STAGE 3 test      -> *.gene_ase.tsv         (ref_is_variable passed through)
IF parental DE provided: STAGE 4 contrast -> *.cis_trans.tsv
STAGE 5 report    -> *.report.html
each stage writes its own TSV, so any stage can be run / inspected alone.
```

---

## The decisions that matter most, in one place

1. **Phase from parents, not from F1 expression** (1.3) — keeps strong-ASE sites
   that RNA genotyping miscalls as homozygous.
2. **`both_hom` needs no F1 genotype** (1.3 D1) — recovers the classic diagnostic
   SNPs (one parent ref-homozygous, other alt-homozygous) under any F1 ASE.
3. **Variable/fixed by parent LINEAGE, not sex** (1.1) — reciprocal-cross safe.
4. **Filter on allele-bearing reads, not raw depth** (2.2) — RNA intron skips
   inflate raw mpileup depth.
5. **Per-plant test + intersection-union max-p** (3.2-3.3) — honest at n=2 plants;
   makes cross-background consistency a requirement, not a footnote.
6. **Flag, don't hide, reference-bias risk** (3.2) — `fixed_allele_seen` /
   `possible_ref_bias` separate real complete-ASE from mapping dropout.
```
