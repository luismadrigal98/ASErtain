"""Informative-SNP discovery for outbred parents, resolved per F1 individual.

With outbred (heterozygous) parents there is rarely a site "fixed between
species". What we need instead is *phase*: at a heterozygous F1 site, which
allele is maternal (variable lineage) and which is paternal (fixed lineage)?
That is resolved per F1 plant from its own genotype plus its two named parents:

    * both parents homozygous for different alleles  -> 'both_hom'  (F1 must be
      het; informative even without the F1 genotype)
    * one parent homozygous + F1 genotyped heterozygous -> 'phased'  (the
      homozygous parent fixes its contributed allele, so the other F1 allele is
      assigned to the opposite parent — this rescues sites where the *other*
      parent is heterozygous, common with outbreeding)
    * both parents heterozygous (or phase otherwise ambiguous) -> uninformative

Because each F1 has its own parents and its own genotype, the variable/fixed
allele identity is tracked **per F1 plant**, not per species. A site informative
and concordant for every F1 plant is 'shared'; otherwise it is 'plant_specific'
and only used for the plants it is valid for.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from .annotation import GeneIndex
from .config import CrossConfig, F1Plant
from .vcf import SampleCall, Variant, iter_variants

HOM_REF, HOM_ALT, HET, MISSING = "hom_ref", "hom_alt", "het", "missing"


@dataclass
class PlantAllele:
    variable: str   # nucleotide inherited from the variable-lineage (maternal) parent
    fixed: str      # nucleotide inherited from the fixed-lineage (paternal) parent
    tier: str       # 'both_hom' | 'phased'


@dataclass
class InformativeSNP:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    per_plant: Dict[str, PlantAllele]   # plant name -> resolved alleles
    classification: str                  # 'shared' | 'plant_specific'
    backgrounds: List[str]               # backgrounds of the informative plants
    gene_id: str = "intergenic"
    gene_name: str = "intergenic"
    location: str = "intergenic"

    def for_plant(self, plant: str) -> Optional[PlantAllele]:
        return self.per_plant.get(plant)

    def variable_is_ref(self, plant: str) -> bool:
        pa = self.per_plant.get(plant)
        return bool(pa and pa.variable == self.ref)


@dataclass
class DiagnoseStats:
    total: int = 0
    biallelic_snp: int = 0
    informative: int = 0
    shared: int = 0
    plant_specific: int = 0


# ---------------------------------------------------------------------------
# Genotype calling (AD-aware, tolerant of low-level sequencing noise)
# ---------------------------------------------------------------------------

def call_state(call: SampleCall, *, min_depth: int, maf_threshold: float) -> str:
    """Return HOM_REF / HOM_ALT / HET / MISSING for one sample."""
    if call.ad is not None and len(call.ad) >= 2:
        ref_d, alt_d = call.ad[0], call.ad[1]
        total = ref_d + alt_d
        if total < min_depth:
            return MISSING
        alt_frac = alt_d / total
        if alt_frac <= maf_threshold:
            return HOM_REF
        if (1 - alt_frac) <= maf_threshold:
            return HOM_ALT
        return HET
    if call.dp is not None and call.dp < min_depth:
        return MISSING
    if call.gt is None:
        return MISSING
    a, b = call.gt
    if a == b:
        return HOM_REF if a == "0" else HOM_ALT
    return HET


def _allele_set(state: str) -> Optional[Set[str]]:
    return {HOM_REF: {"0"}, HOM_ALT: {"1"}, HET: {"0", "1"}}.get(state)


def informative_for_plant(mother: Optional[Set[str]],
                          father: Optional[Set[str]],
                          f1: Optional[Set[str]],
                          ref: str, alt: str) -> Optional[PlantAllele]:
    """Resolve the maternal (variable) and paternal (fixed) allele for one F1.

    `mother`/`father`/`f1` are allele-index sets ({'0'}, {'1'}, {'0','1'}) or
    None when the genotype is missing. Returns a PlantAllele or None.
    """
    def nuc(idx: str) -> str:
        return ref if idx == "0" else alt

    hom = ({"0"}, {"1"})
    f1_het = f1 == {"0", "1"}

    # Strongest tier: both parents homozygous for different alleles.
    if mother in hom and father in hom and mother != father:
        if f1 is not None and not f1_het:
            return None                      # genotype conflicts with expectation
        return PlantAllele(variable=nuc(next(iter(mother))),
                           fixed=nuc(next(iter(father))), tier="both_hom")

    # Phased tiers require a heterozygous F1 genotype.
    if not f1_het:
        return None

    if father in hom:                        # paternal allele fixed -> maternal known
        pat = next(iter(father))
        mat = "1" if pat == "0" else "0"
        if mother is not None and mat not in mother:
            return None                      # mother cannot transmit the maternal allele
        return PlantAllele(variable=nuc(mat), fixed=nuc(pat), tier="phased")

    if mother in hom:                        # maternal allele fixed -> paternal known
        mat = next(iter(mother))
        pat = "1" if mat == "0" else "0"
        if father is not None and pat not in father:
            return None
        return PlantAllele(variable=nuc(mat), fixed=nuc(pat), tier="phased")

    return None                              # both heterozygous / ambiguous phase


# ---------------------------------------------------------------------------
# Site classification across all F1 plants
# ---------------------------------------------------------------------------

def classify_site(variant: Variant, cfg: CrossConfig, *,
                  min_depth: int, maf_threshold: float,
                  stats: DiagnoseStats) -> Optional[InformativeSNP]:
    alt = variant.alt[0]

    # Genotype every parent and F1 once.
    def state(sample: Optional[str]) -> Optional[Set[str]]:
        if not sample:
            return None
        return _allele_set(call_state(variant.call(sample),
                                      min_depth=min_depth,
                                      maf_threshold=maf_threshold))

    parent_states = {p.name: state(p.vcf_sample) for p in cfg.parents}

    per_plant: Dict[str, PlantAllele] = {}
    for pl in cfg.f1_plants:
        pa = informative_for_plant(
            parent_states.get(pl.mother),
            parent_states.get(pl.father),
            state(pl.vcf_sample),
            variant.ref, alt)
        if pa is not None:
            per_plant[pl.name] = pa
    if not per_plant:
        return None

    # 'shared' iff every F1 plant is informative and concordant on the alleles.
    all_plants = {pl.name for pl in cfg.f1_plants}
    vari(0)  # placeholder removed below
    return _build(variant, alt, per_plant, all_plants, cfg, stats)


def _build(variant, alt, per_plant, all_plants, cfg, stats) -> InformativeSNP:
    variable_nucs = {pa.variable for pa in per_plant.values()}
    fixed_nucs = {pa.fixed for pa in per_plant.values()}
    shared = (set(per_plant) == all_plants
              and len(variable_nucs) == 1 and len(fixed_nucs) == 1)
    if shared:
        stats.shared += 1
        classification = "shared"
    else:
        stats.plant_specific += 1
        classification = "plant_specific"
    stats.informative += 1

    bg_by_plant = {pl.name: pl.bg for pl in cfg.f1_plants}
    backgrounds = sorted({bg_by_plant[name] for name in per_plant})

    return InformativeSNP(
        chrom=variant.chrom, pos=variant.pos, ref=variant.ref, alt=alt,
        qual=variant.qual, per_plant=per_plant,
        classification=classification, backgrounds=backgrounds)


def find_informative_snps(cfg: CrossConfig, vcf_path: str, *,
                          min_depth: int = 8,
                          maf_threshold: float = 0.10,
                          min_qual: float = 30.0,
                          snps_only: bool = True,
                          chrom_filter: Optional[str] = None,
                          gene_index: Optional[GeneIndex] = None,
                          ) -> Tuple[List[InformativeSNP], DiagnoseStats]:
    out: List[InformativeSNP] = []
    stats = DiagnoseStats()
    window = cfg.annotation_window
    for v in iter_variants(vcf_path, snps_only=snps_only,
                           chrom_filter=chrom_filter, min_qual=min_qual):
        stats.total += 1
        if not v.is_biallelic_snp:
            continue
        stats.biallelic_snp += 1
        snp = classify_site(v, cfg, min_depth=min_depth,
                            maf_threshold=maf_threshold, stats=stats)
        if snp is None:
            continue
        if gene_index is not None:
            hit = gene_index.annotate(snp.chrom, snp.pos, window=window)
            snp.gene_id, snp.gene_name, snp.location = (
                hit.gene_id, hit.gene_name, hit.location)
        out.append(snp)
    return out, stats
