"""Proofreading tests for the four added features (pure-Python paths only).

Run:  python3 -m pytest tests/ -q     (or)     python3 tests/test_new_features.py

These cover everything that does NOT need samtools/BAMs:
  * Task 1  flower-contribution normalisation (equalize vs none)
  * Task 3  per gene×SNP counts table
  * Task 2  parental-DE statistics + cis/trans DE-concordance sanity check
  * Task 4  label round-trip through the TSV read/write layer
"""
import math
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from asertain.labels import Labels, parse_comment
from asertain import tables, testing, contrast, expression, haplotype


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
    # Two genotypes per lineage (so a valid across-genotype test exists).
    counts = {
        "g1": {"v1a": 1000, "v2a": 1100, "f1a": 100, "f2a": 90},
        "g2": {"v1a": 200, "v2a": 210, "f1a": 195, "f2a": 205},
    }
    lineage = {"v1a": "variable", "v2a": "variable", "f1a": "fixed", "f2a": "fixed"}
    geno = {"v1a": "V1", "v2a": "V2", "f1a": "F1", "f2a": "F2"}   # 4 genotypes
    names = {"g1": "g1", "g2": "g2"}
    libs = {s: 1_000_000 for s in lineage}
    de = expression.differential_expression(counts, lineage, names,
                                            library_sizes=libs, sample_genotype=geno)
    by_gene = {r["gene_id"]: r for r in de}
    assert by_gene["g1"]["log2FoldChange"] > 1.5
    assert by_gene["g1"]["higher_in"] == "variable"
    assert by_gene["g1"]["pvalue"] != "NA"
    assert by_gene["g1"]["method"] == "welch_genotype"
    assert by_gene["g1"]["n_variable"] == 2 and by_gene["g1"]["n_fixed"] == 2
    assert abs(by_gene["g2"]["log2FoldChange"]) < 0.3
    assert set(expression.DE_COLS).issuperset(by_gene["g1"].keys())


def test_parental_de_unequal_flower_counts_do_not_bias_fc():
    # Variable lineage: k2 (4 flowers) + k3 (3 flowers). Fixed: amph (6 flowers).
    # Each genotype's flowers carry the SAME per-genotype level, so collapsing to
    # genotype means must give exactly the genotype-level fold change regardless
    # of how many flowers each genotype has.
    counts = {"g": {}}
    lineage, geno, libs = {}, {}, {}
    def add(sample, gt, lin, level):
        counts["g"][sample] = level
        lineage[sample] = lin; geno[sample] = gt; libs[sample] = 1_000_000
    for i in range(4): add(f"k2_{i}", "k2", "variable", 300)   # k2 level 300
    for i in range(3): add(f"k3_{i}", "k3", "variable", 100)   # k3 level 100
    for i in range(6): add(f"a_{i}", "amph", "fixed", 50)      # amph level 50
    de = expression.differential_expression(counts, lineage, {"g": "g"},
                                            library_sizes=libs, sample_genotype=geno)[0]
    # Genotype means: variable = mean(300,100)=200, fixed = 50.
    # log2((200+1)/(50+1)) ~ log2(3.94) ~ 1.98 -- NOT pulled by amph's 6 flowers.
    assert abs(de["log2FoldChange"] - math.log2(201/51)) < 1e-3   # value is 4-dp rounded
    assert de["n_variable"] == 2 and de["n_fixed"] == 1
    assert de["n_variable_flowers"] == 7 and de["n_fixed_flowers"] == 6
    # Fixed has 1 genotype -> no valid across-genotype test -> flagged fallback.
    assert de["method"] == "welch_flower_pseudorep"


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


# ---------------------------------------------------------------------------
# Read-backed haplotype counting (CIGAR walk + fragment assignment)
# ---------------------------------------------------------------------------

def test_read_bases_at_simple_match():
    # 10M read starting at ref 100 -> covers 100..109. Base at 105 is seq[5].
    seq = "ACGTACGTAC"
    out = haplotype.read_bases_at(100, "10M", seq, "I" * 10, [100, 105, 109, 200])
    assert out[100][0] == "A"
    assert out[105][0] == "C"   # seq[5]
    assert out[109][0] == "C"   # seq[9]
    assert 200 not in out       # not covered


def test_read_bases_at_with_intron_and_indel():
    # 5M 100N 5M: covers ref 100..104 and 205..209 (intron skips 105..204).
    seq = "AAAAACCCCC"
    out = haplotype.read_bases_at(100, "5M100N5M", seq, "I" * 10, [102, 150, 205])
    assert out[102][0] == "A"
    assert 150 not in out                 # falls in the intron (N) -> absent
    assert out[205][0] == "C"             # seq[5], first base after the skip
    # Soft clip consumes query but not ref: 3S7M at ref 100 -> ref 100 = seq[3].
    out2 = haplotype.read_bases_at(100, "3S7M", "TTTACGTAC", "I" * 9, [100])
    assert out2[100][0] == "A"            # seq[3], past the 3 soft-clipped bases


def test_fragment_assignment_logic(tmp_path=None):
    import subprocess, tempfile, os, shutil
    if shutil.which("samtools") is None:
        print("    (skipped: samtools not on PATH)")
        return
    d = tempfile.mkdtemp()
    sam = os.path.join(d, "x.sam")
    # SNP at pos 105 (var=A, fix=G) and pos 115 (var=A, fix=G).
    # frag1: both mates/reads variable -> 1 variable. frag2: fixed at both -> fixed.
    # frag3: variable at 105 but fixed at 115 -> ambiguous (excluded).
    # frag4: single read, variable at 105 only -> variable.
    with open(sam, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:300\n")
        def rd(q, pos, seq):
            fh.write(f"{q}\t0\tchr1\t{pos}\t60\t{len(seq)}M\t*\t0\t0\t{seq}\t{'I'*len(seq)}\n")
        # read spanning 100..119 (20M); index of 105 = 5, 115 = 15.
        rd("frag1", 100, "AAAAA" + "A" + "AAAAAAAAA" + "A" + "AAAA")  # A at 105 and 115
        rd("frag2", 100, "GGGGG" + "G" + "GGGGGGGGG" + "G" + "GGGG")  # G at 105 and 115
        rd("frag3", 100, "AAAAA" + "A" + "AAAAAAAAA" + "G" + "AAAA")  # A@105, G@115
        rd("frag4", 100, "AAAAA" + "A" + "AAAA")                       # covers 100..109, A@105
    bam = os.path.join(d, "x.bam")
    subprocess.run(f"samtools view -bS {sam} | samtools sort -o {bam} -",
                   shell=True, check=True)
    subprocess.run(["samtools", "index", bam], check=True)
    snps = [(105, "A", "G"), (115, "A", "G")]
    var, fix, amb, n = haplotype.count_gene_haplotypes(
        bam, "chr1", 105, 115, snps, min_mapq=20, min_baseq=20)
    assert var == 2          # frag1 + frag4
    assert fix == 1          # frag2
    assert amb == 1          # frag3 (both haplotypes)
    assert n == 4


def test_haplotype_counts_each_read_once_endtoend():
    """A fragment spanning 3 SNPs must be counted once, not 3x; the per-plant
    test is then a binomial over independent reads, and n_snps reports the real
    phased-SNP count."""
    import subprocess, tempfile, os, shutil
    if shutil.which("samtools") is None:
        print("    (skipped: samtools not on PATH)")
        return
    from asertain.config import CrossConfig, Parent, Flower, Reference, F1Plant
    from asertain.genotypes import InformativeSNP, PlantAllele
    d = tempfile.mkdtemp()

    def make_bam(name, n_var, n_fix):
        sam = os.path.join(d, name + ".sam")
        with open(sam, "w") as fh:
            fh.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:300\n")
            rid = 0
            for hap, n in [("A", n_var), ("G", n_fix)]:
                for _ in range(n):
                    s = list("A" * 30)
                    for off in (5, 15, 25):       # SNPs at 105/115/125
                        s[off] = hap
                    fh.write(f"r{name}_{rid}\t0\tchr1\t100\t60\t30M\t*\t0\t0\t"
                             f"{''.join(s)}\t{'I'*30}\n"); rid += 1
        bam = os.path.join(d, name + ".bam")
        subprocess.run(f"samtools view -bS {sam} | samtools sort -o {bam} -",
                       shell=True, check=True)
        subprocess.run(["samtools", "index", bam], check=True)
        return bam

    bams = {n: make_bam(n, v, f) for n, v, f in
            [("p1f1", 80, 20), ("p1f2", 78, 22), ("p2f1", 82, 18), ("p2f2", 79, 21)]}
    parents = [Parent("V", "V", "variable"), Parent("F", "F", "fixed")]
    plants = [
        F1Plant("P1", "V", "F", [Flower("p1f1", bams["p1f1"]), Flower("p1f2", bams["p1f2"])],
                vcf_sample="P1", background="bgV"),
        F1Plant("P2", "V", "F", [Flower("p2f1", bams["p2f1"]), Flower("p2f2", bams["p2f2"])],
                vcf_sample="P2", background="bgF"),
    ]
    cfg = CrossConfig("t", Reference(None, "unknown"), "v", "f", parents, plants)

    def snp(pos):
        pa = {pl.name: PlantAllele("A", "G", "both_hom") for pl in plants}
        return InformativeSNP("chr1", pos, "A", "G", 100.0, pa, "shared",
                              ["bgV", "bgF"], gene_id="G1", gene_name="G1", location="genic")
    recs = haplotype.count_flowers_haplotype(cfg, [snp(105), snp(115), snp(125)],
                                             min_depth=10, progress=False)
    assert len(recs) == 4
    # Each fragment spans all 3 SNPs but is counted ONCE: ~80 variable, not 240.
    by_fl = {r["flower"]: r for r in recs}
    assert by_fl["p1f1"]["variable_count"] == 80 and by_fl["p1f1"]["fixed_count"] == 20
    assert by_fl["p1f1"]["n_hap_snps"] == 3
    genes = testing.test_genes(recs, alpha=0.05)
    g = genes[0]
    assert g["n_snps"] == 3                 # real phased-SNP count, not the pseudo-SNP
    assert g["method"] == "binomial"        # one independent (K,N) per plant
    assert g["ase_call"] is True and g["direction"] == "variable"


def test_ambiguous_fraction_flags_and_blocks_call():
    # A clean, strongly-biased gene with FEW ambiguous reads -> called.
    clean = [_rec(f"f{p}", p, "chr1:100", 90, 10, bg=b) for p, b in
             [("P1", "b1"), ("P2", "b2")]]
    for r in clean:
        r["other_count"] = 2          # ~2% ambiguous
        r["total_depth"] = r["variable_count"] + r["fixed_count"] + 2
    g = testing.test_genes(clean, alpha=0.05)[0]
    assert g["low_ambiguity"] is True and g["ase_call"] is True
    assert g["ambiguous_fraction"] < 0.05 and g["other_reads"] == 4

    # Same signal but 40% of reads ambiguous (mis-phasing) -> flagged, not called.
    noisy = [_rec(f"f{p}", p, "chr1:100", 90, 10, bg=b) for p, b in
             [("P1", "b1"), ("P2", "b2")]]
    for r in noisy:
        r["other_count"] = 67         # 67/(100+67) ~ 0.40
        r["total_depth"] = r["variable_count"] + r["fixed_count"] + 67
    g2 = testing.test_genes(noisy, alpha=0.05)[0]
    assert g2["ambiguous_fraction"] > 0.30
    assert g2["low_ambiguity"] is False
    assert g2["ase_call"] is False        # blocked by the QC flag despite signal
    # A looser ceiling lets it through again.
    g3 = testing.test_genes(noisy, alpha=0.05, max_other_fraction=0.5)[0]
    assert g3["low_ambiguity"] is True and g3["ase_call"] is True


# ---------------------------------------------------------------------------
# Audit fixes — robustness & extensibility
# ---------------------------------------------------------------------------

def _write_vcf(lines):
    d = tempfile.mkdtemp()
    path = os.path.join(d, "v.vcf")
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for ln in lines:
            fh.write(ln + "\n")
    return path


def test_qual_dot_is_not_dropped():
    # QUAL='.' must pass the filter (missing != zero); a numeric QUAL below the
    # threshold must still be dropped.
    from asertain import vcf
    path = _write_vcf([
        "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",     # missing QUAL -> keep
        "chr1\t200\t.\tA\tT\t5\tPASS\t.\tGT\t0/1",     # QUAL 5 -> drop at min 30
        "chr1\t300\t.\tA\tT\t99\tPASS\t.\tGT\t0/1",    # QUAL 99 -> keep
    ])
    keep = list(vcf.iter_variants(path, min_qual=30.0))
    positions = {v.pos for v in keep}
    assert positions == {100, 300}
    missing = next(v for v in keep if v.pos == 100)
    assert missing.qual != missing.qual          # NaN preserved as "missing"


def test_polyploid_and_haploid_gt_calls():
    from asertain.vcf import _parse_gt, SampleCall
    from asertain.genotypes import call_state, HOM_REF, HOM_ALT, HET, MISSING
    assert _parse_gt("0/0/1/1") == ("0", "0", "1", "1")   # tetraploid het
    assert _parse_gt("0") == ("0",)                        # haploid
    assert _parse_gt("0/.") is None                        # half-missing
    sc = lambda gt: SampleCall(gt, None, 50)
    assert call_state(sc(_parse_gt("0/0/1")), min_depth=8, maf_threshold=0.1) == HET
    assert call_state(sc(_parse_gt("1/1/1/1")), min_depth=8, maf_threshold=0.1) == HOM_ALT
    assert call_state(sc(_parse_gt("0/0/0")), min_depth=8, maf_threshold=0.1) == HOM_REF
    assert call_state(sc(None), min_depth=8, maf_threshold=0.1) == MISSING


def test_reference_lineage_resolution():
    from asertain.config import CrossConfig, Parent, F1Plant, Flower, Reference
    parents = [Parent("k2", "k2", "variable"), Parent("amph", "amph", "fixed")]
    pl = [F1Plant("P1", "k2", "amph", [Flower("f", "b")], vcf_sample="P1")]
    cfg = CrossConfig("t", Reference("r.fa", "amph"), "kun", "amph", parents, pl)
    assert cfg.reference_lineage() == "fixed"          # parent name -> its lineage
    cfg.reference.identity = "variable"
    assert cfg.reference_lineage() == "variable"
    cfg.reference.identity = "third_species"
    assert cfg.reference_lineage() is None             # no single-parent bias


def test_possible_ref_bias_is_symmetric():
    # Variable allele never seen (all reads fixed). Under a FIXED reference this
    # is the suspect mapping artefact; under a VARIABLE reference it is not.
    recs = [_rec(f"f{p}", p, "chr1:100", 0, 80, bg=b)
            for p, b in [("P1", "b1"), ("P2", "b2")]]
    g_fix = testing.test_genes(recs, ref_lineage="fixed")[0]
    g_var = testing.test_genes(recs, ref_lineage="variable")[0]
    assert g_fix["variable_allele_seen"] is False
    assert g_fix["possible_ref_bias"] is True
    assert g_var["possible_ref_bias"] is False


def test_null_shift_heterogeneous_nulls_not_pooled():
    # Each SNP sits exactly at its OWN (shifted) null -> no real ASE. The old
    # pool-then-test-vs-averaged-null path could manufacture signal; the per-SNP
    # shift test must report no imbalance.
    def r(snp, v, f, null, p, bg):
        d = _rec(f"fl{p}", p, snp, v, f, bg=bg, tier="phased")
        d["null_p"] = null
        return d
    recs = []
    for p, bg in [("P1", "b1"), ("P2", "b2")]:
        recs += [r("chr1:100", 30, 70, 0.3, p, bg),
                 r("chr1:200", 280, 120, 0.7, p, bg)]
    g = testing.test_genes(recs, alpha=0.05)[0]
    assert g["method"] in ("binom_shift", "betabinom_shift")
    assert g["p_primary"] > 0.5 and g["ase_call"] is False
    # A genuine shift relative to each SNP's null is detected and called.
    recs2 = []
    for p, bg in [("P1", "b1"), ("P2", "b2")]:
        recs2 += [r("chr1:100", 45, 55, 0.3, p, bg),
                  r("chr1:200", 340, 60, 0.7, p, bg)]
    g2 = testing.test_genes(recs2, alpha=0.05)[0]
    assert g2["p_primary"] < 1e-3 and g2["direction"] == "variable"
    assert g2["phase_concordant"] is True and g2["ase_call"] is True


def test_stouffer_combine_scales_power():
    # Six plants all variable-biased, one weak/shallow. max-p is killed by the
    # weak plant; Stouffer aggregates the consistent shift and calls.
    recs = [_rec(f"fl{p}", p, "chr1:100", v, f, bg=b) for p, b, v, f in [
        ("P1", "b1", 60, 40), ("P2", "b2", 58, 42), ("P3", "b3", 62, 38),
        ("P4", "b4", 59, 41), ("P5", "b5", 61, 39), ("P6", "b6", 9, 6)]]
    mp = testing.test_genes(recs, alpha=0.05, combine="maxp")[0]
    so = testing.test_genes(recs, alpha=0.05, combine="stouffer")[0]
    assert mp["ase_call"] is False
    assert so["p_primary"] < mp["p_primary"] and so["ase_call"] is True


def test_annotation_prefers_genic_over_window():
    from asertain.annotation import GeneIndex, Gene
    # geneA at 1000-2000 (listed first), geneB at 2100-3000. A SNP at 2150 is
    # genic in B but within A's downstream 500-window. Must annotate as B/genic.
    idx = GeneIndex({"chr1": [
        Gene("A", "A", "protein", 1000, 2000, "+"),
        Gene("B", "B", "protein", 2100, 3000, "+")]})
    hit = idx.annotate("chr1", 2150, window=500)
    assert hit.gene_id == "B" and hit.location == "genic"
    # A SNP at 2050 is in no gene but in A's downstream window -> A/downstream.
    hit2 = idx.annotate("chr1", 2050, window=500)
    assert hit2.gene_id == "A" and hit2.location == "downstream"


def test_mpileup_uses_uncapped_depth_and_matched_flags():
    from asertain import external
    # The pileup command must disable the depth cap (-d 0) and use the same read
    # exclusion mask as the haplotype counter's samtools_view, so both counters
    # see the same reads.
    assert external._COUNT_EXCLUDE_FLAGS == 3844
    import inspect
    src = inspect.getsource(external.samtools_mpileup)
    assert '"-d"' in src and "_COUNT_EXCLUDE_FLAGS" in src and '"--ff"' in src
    assert '"-B"' in src                          # BAQ disabled for ASE counting


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
