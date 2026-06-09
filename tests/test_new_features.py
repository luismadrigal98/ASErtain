"""Proofreading tests for the four added features (pure-Python paths only).

Run:  python3 -m pytest tests/ -q     (or)     python3 tests/test_new_features.py

These cover everything that does NOT need samtools/BAMs:
  * Task 1  flower-contribution normalisation (equalize vs none)
  * Task 3  per gene×SNP counts table
  * Task 2  parental-DE statistics + cis/trans DE-concordance sanity check
  * Task 4  label round-trip through the TSV read/write layer
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from asertain.labels import Labels, parse_comment
from asertain import tables, testing, contrast, expression


# ---------------------------------------------------------------------------
# Helpers: build a synthetic allele-count record set
# ---------------------------------------------------------------------------

def _rec(flower, plant, snp, v, f, *, bg="bg1", gene="g1", tier="both_hom"):
    return {
        "flower": flower, "plant": plant, "background": bg,
        "chrom": "chr1", "pos": int(snp.split(":")[1]), "snp_id": snp,
        "variable_allele": "A", "fixed_allele": "T", "variable_is_ref": "True",
        "tier": tier, "variable_count": v, "fixed_count": f, "other_count": 0,
        "total_depth": v + f, "null_p": 0.5, "gene_id": gene, "gene_name": gene,
    }


# ---------------------------------------------------------------------------
# Task 1 — flower normalisation
# ---------------------------------------------------------------------------

def test_flower_size_factors_equalize():
    # Plant P1: shallow flower (10x) variable-biased; deep flower (1000x) balanced.
    recs = [
        _rec("f_shallow", "P1", "chr1:100", 9, 1),       # depth 10, ratio 0.9
        _rec("f_deep", "P1", "chr1:100", 500, 500),       # depth 1000, ratio 0.5
    ]
    sf = testing.flower_size_factors(recs, mode="equalize")
    # Deep flower gets a >1 factor, shallow <1; geometric mean of factors ~1.
    assert sf[("P1", "f_deep")] > 1.0
    assert sf[("P1", "f_shallow")] < 1.0

    raw = testing._plant_snp_counts(recs, None)["P1"]
    norm = testing._plant_snp_counts(recs, sf)["P1"]
    raw_ratio = raw[0][0] / raw[0][1]
    norm_ratio = norm[0][0] / norm[0][1]
    # Raw pooling is dominated by the deep balanced flower (~0.5); equalising
    # pulls the pooled ratio toward the average of the two flowers.
    assert abs(raw_ratio - 0.5) < 0.05
    assert norm_ratio > raw_ratio + 0.1


def test_flower_norm_none_is_raw_sum():
    recs = [_rec("f1", "P1", "chr1:100", 9, 1),
            _rec("f2", "P1", "chr1:100", 500, 500)]
    none = testing._plant_snp_counts(recs, testing.flower_size_factors(recs, mode="none"))
    assert none["P1"] == [(509, 1010)]


def test_single_flower_factor_is_one():
    recs = [_rec("only", "P1", "chr1:100", 9, 1)]
    sf = testing.flower_size_factors(recs, mode="equalize")
    assert sf[("P1", "only")] == 1.0


def test_flower_norm_recovers_ase_masked_by_a_deep_flower():
    # Each plant: two shallow strongly-biased flowers + one huge balanced flower.
    # Raw summing is dominated by the deep balanced flower (no call); equalising
    # recovers the real imbalance and makes the call.
    recs = []
    for pl, bg in [("P1", "k2"), ("P2", "k3")]:
        recs.append(_rec(pl + "s1", pl, "chr1:100", 27, 3, bg=bg))
        recs.append(_rec(pl + "s2", pl, "chr1:160", 24, 3, bg=bg))
        recs.append(_rec(pl + "deep", pl, "chr1:100", 2000, 2000, bg=bg))
    raw = testing.test_genes(recs, flower_norm="none")[0]
    eq = testing.test_genes(recs, flower_norm="equalize")[0]
    assert raw["ase_call"] is False and abs(raw["log2_ratio"]) < 0.2
    assert eq["ase_call"] is True and eq["log2_ratio"] > 1.0
    assert eq["direction"] == "variable"


# ---------------------------------------------------------------------------
# Task 3 — per gene×SNP counts table
# ---------------------------------------------------------------------------

def test_snp_gene_summary_collapses_plants_and_flowers():
    recs = [
        _rec("f1", "P1", "chr1:100", 8, 2),
        _rec("f2", "P1", "chr1:100", 6, 4),
        _rec("f1b", "P2", "chr1:100", 5, 5, bg="bg2"),
        _rec("f1", "P1", "chr1:200", 1, 9),
    ]
    rows = testing.snp_gene_summary(recs)
    by_snp = {r["snp_id"]: r for r in rows}
    s1 = by_snp["chr1:100"]
    assert s1["variable_count"] == 8 + 6 + 5
    assert s1["fixed_count"] == 2 + 4 + 5
    assert s1["n_plants"] == 2
    assert s1["n_flowers"] == 3
    assert "P1=14/20" in s1["per_plant_counts"]   # variable=8+6, total=10+10
    assert "P2=5/10" in s1["per_plant_counts"]
    assert set(testing.SNP_GENE_COLS).issuperset(s1.keys())


# ---------------------------------------------------------------------------
# Task 2 — parental DE + DE-concordance sanity check
# ---------------------------------------------------------------------------

def test_parental_de_direction_and_columns():
    # g1: variable parents much higher than fixed.  g2: roughly equal.
    counts = {
        "g1": {"v1": 1000, "v2": 1100, "f1": 100, "f2": 90},
        "g2": {"v1": 200, "v2": 210, "f1": 195, "f2": 205},
    }
    lineage = {"v1": "variable", "v2": "variable", "f1": "fixed", "f2": "fixed"}
    names = {"g1": "g1", "g2": "g2"}
    # Equal library sizes -> library-size normalisation is identity, so the
    # fold change reflects the raw biology (robust for a tiny gene panel).
    libs = {"v1": 1_000_000, "v2": 1_000_000, "f1": 1_000_000, "f2": 1_000_000}
    de = expression.differential_expression(counts, lineage, names, library_sizes=libs)
    by_gene = {r["gene_id"]: r for r in de}
    assert by_gene["g1"]["log2FoldChange"] > 1.5
    assert by_gene["g1"]["higher_in"] == "variable"
    assert by_gene["g1"]["pvalue"] != "NA"
    assert abs(by_gene["g2"]["log2FoldChange"]) < 0.3
    assert set(expression.DE_COLS).issuperset(by_gene["g1"].keys())


def test_parental_de_library_size_correction():
    # Same true expression, but the fixed library is 4x deeper. Without depth
    # normalisation g1 would look fixed-biased; library sizes must correct it.
    counts = {"g1": {"v1": 100, "v2": 110, "f1": 400, "f2": 440}}
    lineage = {"v1": "variable", "v2": "variable", "f1": "fixed", "f2": "fixed"}
    libs = {"v1": 1_000_000, "v2": 1_000_000, "f1": 4_000_000, "f2": 4_000_000}
    de = expression.differential_expression(counts, lineage, {"g1": "g1"},
                                            library_sizes=libs)
    assert abs(de[0]["log2FoldChange"]) < 0.2   # depth-corrected -> balanced


def test_size_factors_fallback_when_no_full_genes():
    import numpy as np
    # Every gene has a zero somewhere -> median-of-ratios unusable -> library size.
    mat = np.array([[10.0, 0.0], [0.0, 20.0]])
    sf = expression.size_factors(mat)
    assert sf.shape == (2,)
    assert np.all(sf > 0)


def test_de_sanity_check_concordant_and_discordant():
    # ASE variable-biased (cis>0), DE variable-higher (parental>0) -> concordant.
    conc, verdict = contrast.de_sanity_check(0.8, 1.2, cis_sig=True, parental_sig=True)
    assert conc is True and verdict == "concordant"
    # ASE variable-biased but DE fixed-higher -> discordant (compensatory).
    disc, verdict = contrast.de_sanity_check(0.8, -1.2, cis_sig=True, parental_sig=True)
    assert disc is False and verdict == "discordant_compensatory"
    # Not an ASE candidate.
    na, verdict = contrast.de_sanity_check(0.8, 1.2, cis_sig=False, parental_sig=True)
    assert na == "NA" and verdict == "not_ase"
    # DE not significant.
    na, verdict = contrast.de_sanity_check(0.8, 0.1, cis_sig=True, parental_sig=False)
    assert na == "NA" and verdict == "de_not_sig"


def test_run_contrast_adds_sanity_columns():
    genes = [{"gene_id": "g1", "gene_name": "g1", "log2_ratio": "1.0",
              "q_value": "0.001", "ase_call": "True"}]
    de = [{"gene_id": "g1", "log2FoldChange": "2.0", "padj": "0.001"}]
    rows = contrast.run_contrast(genes, de)
    assert rows[0]["sanity_check"] == "concordant"
    assert rows[0]["de_concordant"] is True
    assert "de_concordant" in contrast.CONTRAST_COLS


# ---------------------------------------------------------------------------
# Task 4 — label round-trip
# ---------------------------------------------------------------------------

def test_label_name_and_value_mapping():
    lab = Labels("kunthii", "amphorellae")
    assert lab.to_display("variable_count") == "kunthii_count"
    assert lab.to_display("n_plants_fixed_seen") == "n_plants_amphorellae_seen"
    assert lab.to_canonical("kunthii_count") == "variable_count"
    assert lab.to_canonical(lab.to_display("mean_variable_ratio")) == "mean_variable_ratio"
    assert lab.value_to_display("variable") == "kunthii"
    assert lab.value_to_canonical("amphorellae") == "fixed"


def test_allele_counts_label_roundtrip(tmp_path=None):
    lab = Labels("kunthii", "amphorellae")
    recs = [_rec("f1", "P1", "chr1:100", 8, 2)]
    d = tempfile.mkdtemp()
    path = os.path.join(d, "ac.tsv")
    tables.write_allele_counts(recs, path, bias_mode="report", labels=lab)
    with open(path) as fh:
        head = fh.read()
    # Header carries the labels; columns are relabelled.
    assert "variable_label: kunthii" in head
    assert "kunthii_count" in head and "amphorellae_count" in head
    assert "variable_count" not in head.split("\n")[4]  # header line relabelled
    # Reader maps everything back to canonical.
    back = tables.read_allele_counts(path)
    assert back[0]["variable_count"] == 8
    assert back[0]["fixed_count"] == 2
    assert parse_comment(path) == lab


def test_gene_table_direction_value_roundtrip():
    lab = Labels("kunthii", "amphorellae")
    rows = [{"gene_id": "g1", "direction": "variable", "log2_ratio": 1.0}]
    d = tempfile.mkdtemp()
    path = os.path.join(d, "g.tsv")
    tables.write_table(rows, ["gene_id", "direction", "log2_ratio"], path,
                       comment="x", labels=lab)
    with open(path) as fh:
        body = fh.read()
    assert "kunthii" in body  # direction VALUE relabelled
    back = tables.read_table(path)
    assert back[0]["direction"] == "variable"   # canonicalised on read


def test_external_table_without_header_untouched():
    # A DESeq2-style table with no ASErtain label header is read verbatim.
    d = tempfile.mkdtemp()
    path = os.path.join(d, "de.tsv")
    with open(path, "w") as fh:
        fh.write("gene_id\tlog2FoldChange\tpadj\n")
        fh.write("g1\t2.0\t0.01\n")
    rows = tables.read_table(path)
    assert rows[0]["gene_id"] == "g1" and rows[0]["log2FoldChange"] == "2.0"


def _run_all():
    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
        passed += 1
    print(f"\n{passed}/{len(fns)} tests passed")


if __name__ == "__main__":
    _run_all()
