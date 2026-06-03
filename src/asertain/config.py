"""Cross-design configuration: who is a parent, who is an F1, how F1s group.

The whole point of ASErtain's design is that your *cross structure* drives SNP
selection and statistics. A single config file (YAML or JSON) names the exact
parent plants and the F1 replicates, and groups F1s by which parent-of-the-
variable-species (here: kunthii) they descend from. Everything downstream reads
this object, so the CLI stays terse.

Terminology (generic, with the Penstemon mapping in parentheses):
    variable_species  the up-/down-regulated, multi-parent species (kunthii)
    fixed_species     the single-parent reference species (amphorellae)
    backgrounds       F1 groupings by variable-species parent (k1, k2, ...)
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class Parent:
    name: str                 # short label, e.g. 'k1'
    vcf_sample: str           # the column name in the VCF
    species: str              # 'variable' or 'fixed'


@dataclass
class F1Replicate:
    name: str                 # short label, e.g. 'f1_k1_r1'
    vcf_sample: Optional[str] # VCF column (may be None if F1s not genotyped)
    bam: str                  # path to sorted, indexed BAM
    background: str           # which variable-species parent (e.g. 'k1')


@dataclass
class Reference:
    fasta: Optional[str] = None
    # which biological entity the reference equals, for ref-bias bookkeeping:
    # 'variable' | 'fixed' | 'third_species' | 'unknown'
    identity: str = "unknown"


@dataclass
class CrossConfig:
    project: str
    reference: Reference
    variable_label: str               # display name, e.g. 'kunthii'
    fixed_label: str                  # display name, e.g. 'amphorellae'
    variable_parents: List[Parent]    # the exact variable-species plants (k1, k2)
    fixed_parents: List[Parent]       # the exact fixed-species plant(s) (amphorellae)
    f1: List[F1Replicate]
    gtf: Optional[str] = None
    annotation_window: int = 500
    # convenience: derived at load
    backgrounds: Dict[str, List[str]] = field(default_factory=dict)

    # -- lookups -----------------------------------------------------------
    @property
    def variable_vcf_samples(self) -> List[str]:
        return [p.vcf_sample for p in self.variable_parents]

    @property
    def fixed_vcf_samples(self) -> List[str]:
        return [p.vcf_sample for p in self.fixed_parents]

    def f1_by_name(self, name: str) -> F1Replicate:
        for r in self.f1:
            if r.name == name:
                return r
        raise KeyError(name)

    def backgrounds_present(self) -> List[str]:
        return sorted({r.background for r in self.f1})

    def reference_is(self, species: str) -> Optional[bool]:
        """True/False if reference identity is known to (mis)match `species`
        ('variable' or 'fixed'); None if unknown."""
        if self.reference.identity in ("unknown", ""):
            return None
        if self.reference.identity == "third_species":
            return False
        return self.reference.identity == species


def _load_raw(path: str) -> dict:
    with open(path) as fh:
        if path.endswith((".yaml", ".yml")):
            import yaml  # optional dependency; ubiquitous in bioinfo envs
            return yaml.safe_load(fh)
        return json.load(fh)


def load_config(path: str) -> CrossConfig:
    """Parse a YAML/JSON cross-design file into a CrossConfig."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    raw = _load_raw(path)

    ref = Reference(**(raw.get("reference") or {}))

    def _parents(key: str, species: str) -> List[Parent]:
        return [Parent(name=p["name"], vcf_sample=p["vcf_sample"],
                       species=species)
                for p in raw.get(key, [])]

    f1 = [F1Replicate(name=r["name"], vcf_sample=r.get("vcf_sample"),
                      bam=r["bam"], background=r["background"])
          for r in raw.get("f1", [])]

    cfg = CrossConfig(
        project=raw.get("project", "asertain"),
        reference=ref,
        variable_label=raw.get("variable_label", "variable"),
        fixed_label=raw.get("fixed_label", "fixed"),
        variable_parents=_parents("variable_parents", "variable"),
        fixed_parents=_parents("fixed_parents", "fixed"),
        f1=f1,
        gtf=raw.get("gtf"),
        annotation_window=int(raw.get("annotation_window", 500)),
    )
    cfg.backgrounds = {bg: [r.name for r in cfg.f1 if r.background == bg]
                       for bg in cfg.backgrounds_present()}
    _validate(cfg)
    return cfg


def _validate(cfg: CrossConfig) -> None:
    errs: List[str] = []
    if len(cfg.variable_parents) < 1:
        errs.append("at least one variable-species parent is required")
    if len(cfg.fixed_parents) < 1:
        errs.append("at least one fixed-species parent is required")
    if not cfg.f1:
        errs.append("no F1 replicates defined")
    for r in cfg.f1:
        if r.background not in {p.name for p in cfg.variable_parents}:
            errs.append(
                f"F1 '{r.name}' background '{r.background}' does not match any "
                f"variable parent name {[p.name for p in cfg.variable_parents]}")
    if errs:
        raise ValueError("Invalid config:\n  - " + "\n  - ".join(errs))
