"""User-facing labels for the two lineage roles, applied across every output.

Internally ASErtain always works in terms of the canonical roles ``variable``
and ``fixed`` (column names like ``variable_count``, ``fixed_reads`` and the
``direction`` values ``variable`` / ``fixed``). That keeps the inter-stage code
stable. But scientists read the tables, and "variable / fixed" is opaque — they
want to see *kunthii* and *amphorellae*.

This module is the single place that maps between the two:

    * on WRITE  — canonical column names and ``direction`` values are rewritten
      to the user's labels, and the labels are stamped into the file's ``#``
      header so the meaning travels with the data.
    * on READ   — the header labels are parsed and the columns/values are mapped
      back to canonical, so downstream stages keep seeing ``variable``/``fixed``.

Because the labels live in the file header, a stage run on its own (e.g.
``asertain test`` on a hand-made counts file) recovers them without needing the
config. Files that carry no label header (an external DESeq2 DE table, say) are
left exactly as-is.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

# Columns whose *values* (not just their names) are role tokens and must be
# relabelled too. Keep this list in sync with any new direction-valued column.
DIRECTION_COLS = ("direction",)


@dataclass(frozen=True)
class Labels:
    """Display names for the variable- and fixed-lineage roles."""
    variable: str = "variable"
    fixed: str = "fixed"

    @property
    def is_default(self) -> bool:
        return self.variable == "variable" and self.fixed == "fixed"

    @classmethod
    def from_config(cls, cfg) -> "Labels":
        return cls(variable=cfg.variable_label, fixed=cfg.fixed_label)

    # -- name mapping (column headers) -------------------------------------
    def to_display(self, name: str) -> str:
        """canonical column name -> labelled column name."""
        if self.is_default:
            return name
        return name.replace("variable", self.variable).replace("fixed", self.fixed)

    def to_canonical(self, name: str) -> str:
        """labelled column name -> canonical column name."""
        if self.is_default:
            return name
        out = name
        # Replace the longer label first so that a label which is a substring of
        # the other (pathological, but cheap to guard) cannot be clobbered.
        for canon, lab in sorted((("variable", self.variable), ("fixed", self.fixed)),
                                 key=lambda kv: -len(kv[1])):
            out = out.replace(lab, canon)
        return out

    # -- value mapping (direction cells) -----------------------------------
    def value_to_display(self, value: str) -> str:
        if self.is_default:
            return value
        if value == "variable":
            return self.variable
        if value == "fixed":
            return self.fixed
        return value

    def value_to_canonical(self, value: str) -> str:
        if self.is_default:
            return value
        if value == self.variable:
            return "variable"
        if value == self.fixed:
            return "fixed"
        return value


def format_comment(labels: Labels) -> str:
    """The two header lines that stamp the labels into a written file."""
    return (f"# variable_label: {labels.variable}\n"
            f"# fixed_label: {labels.fixed}\n")


def parse_comment(path: str) -> Optional[Labels]:
    """Recover the labels from a file's ``#`` header, or None if absent."""
    var = fix = None
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if "variable_label:" in line:
                var = line.split(":", 1)[1].strip()
            elif "fixed_label:" in line:
                fix = line.split(":", 1)[1].strip()
    if var is None and fix is None:
        return None
    return Labels(variable=var or "variable", fixed=fix or "fixed")


def display_columns(cols: List[str], labels: Labels) -> List[str]:
    return [labels.to_display(c) for c in cols]
