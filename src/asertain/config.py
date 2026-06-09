"""Cross-design configuration: pedigree + nested replication.

The data model reflects two realities of real hybrid crosses:

* **Outbred parents.** Diagnostic sites are not "fixed between species" but
  *informative for a specific F1*: phase (which allele is maternal vs paternal)
  is resolved from that F1's own genotype together with its two named parents.
  So the config names, for every F1, its exact mother and father.

* **Nested replication.** RNA samples (flowers) are grouped under the F1 plant
  they came from. The plant is the biological replicate; flowers are technical /
  observational sub-samples nested within it. The config encodes that hierarchy
  so the statistics never mistake flowers for independent replicates.

Terminology (generic — assign real labels in the config):
    variable lineage  the parental lineage carrying the expression difference of
                      interest, often with several cross parents
    fixed lineage     the other parental lineage
    background        an F1 grouping, by default its mother (variable-lineage
                      parent), used for the cross-background consistency check

Two schemas are accepted:
    * new  : `parents:` + `f1_plants:` (with `flowers:`)  ← full pedigree/nesting
    * legacy: `variable_parents:`/`fixed_parents:`/`f1:`   ← auto-adapted, one
              flower per plant, mother = background, father = first fixed parent
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

VARIABLE, FIXED = "variable", "fixed"


@dataclass
class Parent:
    name: str                 # short label, e.g. 'k1'
    vcf_sample: str           # the column name in the VCF
    lineage: str              # 'variable' or 'fixed'
    # Optional per-sample RNA libraries for this parent (one BAM per flower /
    # library). Only needed for the parental differential-expression stage; the
    # genotyping stages use the VCF column, not these BAMs.
    flowers: List["Flower"] = field(default_factory=list)


@dataclass
class Flower:
    name: str                 # RNA sample / library label
    bam: str                  # path to sorted, indexed BAM
    vcf_sample: Optional[str] = None   # optional, rarely needed at flower level


@dataclass
class F1Plant:
    name: str                 # biological replicate label
    mother: str               # variable-lineage Parent.name
    father: str               # fixed-lineage Parent.name
    flowers: List[Flower]
    vcf_sample: Optional[str] = None   # F1 plant's own VCF column (for phasing)
    background: Optional[str] = None    # defaults to `mother`

    @property
    def bg(self) -> str:
        return self.background or self.mother


@dataclass
class Reference:
    fasta: Optional[str] = None
    # which biological entity the reference equals, for bias bookkeeping:
    # 'variable' | 'fixed' | 'third_species' | a specific parent name | 'unknown'
    identity: str = "unknown"


@dataclass
class CrossConfig:
    project: str
    reference: Reference
    variable_label: str
    fixed_label: str
    parents: List[Parent]
    f1_plants: List[F1Plant]
    gtf: Optional[str] = None
    annotation_window: int = 500

    # -- parent lookups ----------------------------------------------------
    def parent(self, name: str) -> Parent:
        for p in self.parents:
            if p.name == name:
                return p
        raise KeyError(f"parent '{name}' not defined")

    @property
    def variable_parents(self) -> List[Parent]:
        return [p for p in self.parents if p.lineage == VARIABLE]

    @property
    def fixed_parents(self) -> List[Parent]:
        return [p for p in self.parents if p.lineage == FIXED]

    @property
    def parental_flowers(self) -> List["Flower"]:
        """All parental RNA libraries (flowers), across every parent."""
        return [fl for p in self.parents for fl in p.flowers]

    def has_parental_expression(self) -> bool:
        """True if at least one parent of each lineage declares RNA libraries,
        so a variable-vs-fixed parental DE can be computed."""
        var = any(p.flowers for p in self.variable_parents)
        fix = any(p.flowers for p in self.fixed_parents)
        return var and fix

    def lineage_of_parent_flower(self, flower_name: str) -> Optional[str]:
        for p in self.parents:
            if any(fl.name == flower_name for fl in p.flowers):
                return p.lineage
        return None

    # -- F1 lookups --------------------------------------------------------
    @property
    def flowers(self) -> List[Flower]:
        return [fl for pl in self.f1_plants for fl in pl.flowers]

    def plant_of_flower(self, flower_name: str) -> F1Plant:
        for pl in self.f1_plants:
            if any(fl.name == flower_name for fl in pl.flowers):
                return pl
        raise KeyError(flower_name)

    def backgrounds_present(self) -> List[str]:
        return sorted({pl.bg for pl in self.f1_plants})

    def plants_in_background(self, bg: str) -> List[F1Plant]:
        return [pl for pl in self.f1_plants if pl.bg == bg]

    def variable_parent_of(self, plant: F1Plant) -> str:
        """The variable-lineage parent of this F1 (by lineage, not by sex)."""
        for n in (plant.mother, plant.father):
            if self.parent(n).lineage == VARIABLE:
                return n
        raise KeyError(f"F1 plant '{plant.name}' has no variable-lineage parent")

    def fixed_parent_of(self, plant: F1Plant) -> str:
        """The fixed-lineage parent of this F1 (by lineage, not by sex)."""
        for n in (plant.mother, plant.father):
            if self.parent(n).lineage == FIXED:
                return n
        raise KeyError(f"F1 plant '{plant.name}' has no fixed-lineage parent")


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _load_raw(path: str) -> dict:
    with open(path) as fh:
        if path.endswith((".yaml", ".yml")):
            import yaml
            return yaml.safe_load(fh)
        return json.load(fh)


def load_config(path: str) -> CrossConfig:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    raw = _load_raw(path)
    ref = Reference(**(raw.get("reference") or {}))
    common = dict(
        project=raw.get("project", "asertain"),
        reference=ref,
        variable_label=raw.get("variable_label", "variable"),
        fixed_label=raw.get("fixed_label", "fixed"),
        gtf=raw.get("gtf"),
        annotation_window=int(raw.get("annotation_window", 500)),
    )

    if "f1_plants" in raw:
        cfg = _load_new(raw, common)
    elif "f1" in raw:
        cfg = _load_legacy(raw, common)
    else:
        raise ValueError("Config must contain either 'f1_plants' (new schema) "
                         "or 'f1' (legacy schema).")
    _validate(cfg)
    return cfg


def _load_new(raw: dict, common: dict) -> CrossConfig:
    parents = [Parent(name=p["name"], vcf_sample=p["vcf_sample"],
                      lineage=p.get("lineage", VARIABLE),
                      flowers=[Flower(name=fl["name"], bam=fl["bam"],
                                      vcf_sample=fl.get("vcf_sample"))
                               for fl in p.get("flowers", [])])
               for p in raw.get("parents", [])]
    plants: List[F1Plant] = []
    for pl in raw["f1_plants"]:
        flowers = [Flower(name=fl["name"], bam=fl["bam"],
                          vcf_sample=fl.get("vcf_sample"))
                   for fl in pl.get("flowers", [])]
        plants.append(F1Plant(
            name=pl["name"], mother=pl["mother"], father=pl["father"],
            flowers=flowers, vcf_sample=pl.get("vcf_sample"),
            background=pl.get("background")))
    return CrossConfig(parents=parents, f1_plants=plants, **common)


def _load_legacy(raw: dict, common: dict) -> CrossConfig:
    """Adapt the original schema: each F1 becomes a one-flower plant."""
    parents: List[Parent] = []
    for p in raw.get("variable_parents", []):
        parents.append(Parent(p["name"], p["vcf_sample"], VARIABLE))
    for p in raw.get("fixed_parents", []):
        parents.append(Parent(p["name"], p["vcf_sample"], FIXED))
    father = next((p.name for p in parents if p.lineage == FIXED), None)
    plants: List[F1Plant] = []
    for r in raw["f1"]:
        bg = r["background"]
        plants.append(F1Plant(
            name=r["name"], mother=bg, father=father,
            flowers=[Flower(name=r["name"], bam=r["bam"],
                            vcf_sample=r.get("vcf_sample"))],
            vcf_sample=r.get("vcf_sample"), background=bg))
    return CrossConfig(parents=parents, f1_plants=plants, **common)


def _validate(cfg: CrossConfig) -> None:
    errs: List[str] = []
    if not cfg.variable_parents:
        errs.append("no variable-lineage parents defined")
    if not cfg.fixed_parents:
        errs.append("no fixed-lineage parents defined")
    if not cfg.f1_plants:
        errs.append("no F1 plants defined")
    names = {p.name for p in cfg.parents}
    lineage = {p.name: p.lineage for p in cfg.parents}
    for pl in cfg.f1_plants:
        ok = True
        if pl.mother not in names:
            errs.append(f"F1 plant '{pl.name}': mother '{pl.mother}' is not a parent")
            ok = False
        if pl.father not in names:
            errs.append(f"F1 plant '{pl.name}': father '{pl.father}' is not a parent")
            ok = False
        # The two parents must be one variable + one fixed lineage, in either
        # sex role — otherwise variable/fixed alleles would be assigned wrong
        # (audit M2). This makes reciprocal crosses safe and forbids same-lineage
        # 'F1's being analysed as inter-lineage hybrids.
        if ok:
            ls = sorted({lineage[pl.mother], lineage[pl.father]})
            if ls != [FIXED, VARIABLE]:
                errs.append(
                    f"F1 plant '{pl.name}': parents must be one '{VARIABLE}' and "
                    f"one '{FIXED}' lineage; got mother={lineage[pl.mother]}, "
                    f"father={lineage[pl.father]}")
        if not pl.flowers:
            errs.append(f"F1 plant '{pl.name}' has no flowers")
    if errs:
        raise ValueError("Invalid config:\n  - " + "\n  - ".join(errs))
