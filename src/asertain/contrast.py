"""cis / trans decomposition — the biological pay-off of the F1 design.

In an F1 both alleles share one nucleus, so the F1 allelic ratio isolates the
*cis* component of regulatory divergence, while the parental expression
difference reflects *cis + trans*. Therefore:

    cis    = log2 F1 allelic ratio        (variable / fixed)
    total  = log2 parental expression ratio (variable / fixed)   [from your DE]
    trans  = total - cis

Genes are placed into the McManus et al. (2010) / Landry framework categories.
This module is a working scaffold: the category logic is implemented, but the
formal *trans* significance test (a parent-vs-F1 interaction on allele counts,
or comparing F1 allelic counts to parental allele-equivalent counts) is left as
a clearly marked extension — by default we use a sign/҂threshold approximation
and flag it.
"""
from __future__ import annotations

from typing import Dict, List, Optional

CONTRAST_COLS = [
    "gene_id", "gene_name",
    "cis_log2", "parental_log2", "trans_log2",
    "cis_sig", "parental_sig", "trans_sig", "category",
    "de_concordant", "sanity_check",
]


def de_sanity_check(cis_log2: Optional[float], parental_log2: Optional[float],
                    *, cis_sig: bool, parental_sig: bool,
                    tol: float = 0.0) -> tuple:
    """Sanity-check an ASE candidate against the parental DE direction.

    The biological expectation: if a gene is a real ASE candidate, the F1's
    allelic shift should point toward the parent that the DE analysis says is
    more highly expressed (i.e. the *cis* effect should agree in sign with the
    parental expression difference). A significant ASE call pointing the OTHER
    way is not a bug — it is the signature of opposing cis/trans (compensatory)
    regulation — but it should be surfaced, not assumed away.

    Returns (de_concordant, verdict):
      de_concordant : True / False / "NA"
      verdict       : 'concordant' | 'discordant_compensatory'
                      | 'de_not_sig' | 'no_parental_data' | 'not_ase'
    """
    if not cis_sig:
        return "NA", "not_ase"
    if parental_log2 is None:
        return "NA", "no_parental_data"
    if not parental_sig:
        return "NA", "de_not_sig"
    cis_dir = (cis_log2 or 0.0) > tol
    de_dir = parental_log2 > tol
    if cis_dir == de_dir:
        return True, "concordant"
    return False, "discordant_compensatory"


def classify(cis_log2: float, parental_log2: float, *,
             cis_sig: bool, parental_sig: bool, trans_sig: bool,
             tol: float = 0.0) -> str:
    """Return a McManus-style regulatory category."""
    trans_log2 = parental_log2 - cis_log2

    if not parental_sig and not cis_sig and not trans_sig:
        return "conserved"
    if cis_sig and not trans_sig:
        return "cis_only"
    if trans_sig and not cis_sig:
        return "trans_only"
    if cis_sig and trans_sig:
        same_dir = (cis_log2 > tol) == (trans_log2 > tol)
        if not parental_sig:
            # cis and trans cancel out in the parents
            return "compensatory"
        return "cis_plus_trans" if same_dir else "cis_x_trans_opposing"
    return "ambiguous"


def run_contrast(gene_ase: List[Dict], parental_de: List[Dict], *,
                 de_gene_col: str = "gene_id",
                 de_log2_col: str = "log2FoldChange",
                 de_padj_col: str = "padj",
                 ase_alpha: float = 0.05,
                 de_alpha: float = 0.05,
                 trans_log2_threshold: float = 1.0) -> List[Dict]:
    """Join gene-level ASE (cis) with parental DE (total) and classify.

    `parental_de` rows must give a log2 fold change oriented as
    variable/fixed (variable-species over fixed-species) and an adjusted p-value.
    """
    de_index: Dict[str, Dict] = {row[de_gene_col]: row for row in parental_de}
    out: List[Dict] = []

    for g in gene_ase:
        gid = g["gene_id"]
        de = de_index.get(gid)
        cis = _as_float(g.get("log2_ratio"))
        if cis is None:
            continue
        cis_sig = _as_float(g.get("q_value"), 1.0) < ase_alpha and g.get("ase_call") in (True, "True")

        if de is None:
            parental_log2, parental_sig = None, False
        else:
            parental_log2 = _as_float(de.get(de_log2_col))
            parental_sig = _as_float(de.get(de_padj_col), 1.0) < de_alpha

        if parental_log2 is None:
            # No parental DE for this gene: report cis only.
            concordant, verdict = de_sanity_check(
                cis, None, cis_sig=cis_sig, parental_sig=False)
            out.append({
                "gene_id": gid, "gene_name": g.get("gene_name", gid),
                "cis_log2": round(cis, 4), "parental_log2": "NA",
                "trans_log2": "NA", "cis_sig": cis_sig,
                "parental_sig": False, "trans_sig": "NA",
                "category": "cis_only" if cis_sig else "no_parental_data",
                "de_concordant": concordant, "sanity_check": verdict,
            })
            continue

        trans = parental_log2 - cis
        # Approximate trans significance: magnitude threshold (FLAGGED).
        trans_sig = abs(trans) >= trans_log2_threshold
        category = classify(cis, parental_log2,
                            cis_sig=cis_sig, parental_sig=parental_sig,
                            trans_sig=trans_sig)
        concordant, verdict = de_sanity_check(
            cis, parental_log2, cis_sig=cis_sig, parental_sig=parental_sig)
        out.append({
            "gene_id": gid, "gene_name": g.get("gene_name", gid),
            "cis_log2": round(cis, 4),
            "parental_log2": round(parental_log2, 4),
            "trans_log2": round(trans, 4),
            "cis_sig": cis_sig, "parental_sig": parental_sig,
            "trans_sig": trans_sig, "category": category,
            "de_concordant": concordant, "sanity_check": verdict,
        })
    return out


def _as_float(v, default: Optional[float] = None) -> Optional[float]:
    try:
        if v is None or v == "NA" or v == "":
            return default
        return float(v)
    except (TypeError, ValueError):
        return default
