"""Thin wrappers around external bioinformatics tools invoked via subprocess.

Everything that shells out (samtools, bcftools, GATK, WASP, bgzip/tabix) lives
here so the rest of the package stays pure-Python and testable. Tool paths are
configurable so the same code runs on a laptop or an HPC module environment.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from typing import List, Optional, Sequence


class ExternalToolError(RuntimeError):
    """Raised when an external command fails."""


def which(tool: str) -> Optional[str]:
    """Return the resolved path of `tool` or None if not on PATH."""
    return shutil.which(tool)


def check_tool(tool_path: str, version_flag: str = "--version") -> bool:
    """Return True if `tool_path` is runnable; print a one-line version banner."""
    try:
        res = subprocess.run(
            [tool_path, version_flag],
            capture_output=True, text=True, check=False,
            errors="replace",
        )
        if res.returncode == 0:
            first = (res.stdout or res.stderr).strip().split("\n")[0]
            print(f"  ✓ {os.path.basename(tool_path)}: {first}")
            return True
    except FileNotFoundError:
        pass
    print(f"  ✗ {os.path.basename(tool_path)} not found ({tool_path})")
    return False


def require_tools(**tools: str) -> None:
    """Validate a set of {label: path} tools, raising if any is missing.

    Example: require_tools(samtools="samtools", bcftools="bcftools")
    """
    missing = [f"{label} ({path})" for label, path in tools.items()
               if which(path) is None and not os.path.exists(path)]
    if missing:
        raise ExternalToolError(
            "Required external tool(s) not found: " + ", ".join(missing)
        )


def run(cmd: Sequence[str], *, capture: bool = True,
        check: bool = True, text: bool = True) -> subprocess.CompletedProcess:
    """Run a command, raising ExternalToolError with stderr on failure."""
    try:
        return subprocess.run(
            list(cmd), capture_output=capture, text=text, check=check,
            errors="replace" if text else None,
        )
    except subprocess.CalledProcessError as exc:
        raise ExternalToolError(
            f"Command failed (exit {exc.returncode}): {' '.join(map(str, cmd))}\n"
            f"{exc.stderr or ''}"
        ) from exc


# ---------------------------------------------------------------------------
# samtools
# ---------------------------------------------------------------------------

# Read flags excluded by both counters, so --counter pileup and --counter
# haplotype see the SAME reads. 3844 = UNMAP(0x4) + SECONDARY(0x100) +
# QCFAIL(0x200) + DUPLICATE(0x400) + SUPPLEMENTARY(0x800). samtools mpileup's
# own default --ff omits SUPPLEMENTARY, so we set it explicitly to match
# samtools_view's -F 3844.
_COUNT_EXCLUDE_FLAGS = 3844


def samtools_mpileup(bam: str, region: str, *,
                     min_mapq: int, min_baseq: int,
                     reference: Optional[str] = None,
                     max_depth: int = 0,
                     samtools: str = "samtools") -> str:
    """Return the raw mpileup stdout for a single region (e.g. 'chr1:100-100').

    `max_depth=0` disables samtools' per-file depth cap. The samtools default
    (8000) silently downsamples deep loci — and the downsampling is positional,
    not allele-randomised — so it can bias the allelic ratio at exactly the
    highest-expression genes. We disable it by default for unbiased counting;
    pass a positive value to re-impose a cap.

    BAQ (`-B` disables it) is intentionally OFF: with a reference, samtools
    recomputes base qualities (BAQ) by default, which downweights bases near
    indels asymmetrically and can drop legitimate alt-allele reads — a known
    reference-bias confounder for allele-specific counting. Standard ASE counters
    (GATK ASEReadCounter, phASER) likewise disable it.
    """
    cmd = [samtools, "mpileup", "-r", region,
           "-q", str(min_mapq), "-Q", str(min_baseq),
           "-d", str(max_depth),
           "--ff", str(_COUNT_EXCLUDE_FLAGS),
           "-B",                          # disable BAQ
           "--no-output-ins", "--no-output-del"]
    if reference:
        cmd += ["-f", reference]
    cmd.append(bam)
    res = run(cmd)
    return res.stdout.strip()


def samtools_count(bam: str, region: str, *,
                   min_mapq: int = 20,
                   exclude_flags: int = 0x900,
                   samtools: str = "samtools") -> int:
    """Count primary alignments overlapping `region` (e.g. 'chr1:100-2000').

    Used by the parental-expression stage as a lightweight gene read count.
    `exclude_flags` defaults to 0x900 = secondary (0x100) + supplementary
    (0x800), so each read is counted once. This is a gene-*region* count (it
    includes intronic overlap), a deliberate pure-`samtools` proxy for a proper
    exon-union count — adequate for the cis/trans sanity check, documented as
    approximate.
    """
    cmd = [samtools, "view", "-c", "-q", str(min_mapq),
           "-F", str(exclude_flags), bam, region]
    res = run(cmd)
    out = res.stdout.strip()
    return int(out) if out else 0


def samtools_view(bam: str, region: str, *,
                  min_mapq: int = 20,
                  exclude_flags: int = _COUNT_EXCLUDE_FLAGS,
                  samtools: str = "samtools") -> str:
    """Return raw SAM alignment lines overlapping `region` (no header).

    Used by the read-backed haplotype counter. `exclude_flags` defaults to 3844
    = unmapped(0x4) + secondary(0x100) + qcfail(0x200) + duplicate(0x400) +
    supplementary(0x800), so only primary, mapped, non-duplicate alignments are
    returned and each fragment's two mates are the only repeats (deduplicated by
    QNAME downstream).
    """
    cmd = [samtools, "view", "-q", str(min_mapq),
           "-F", str(exclude_flags), bam, region]
    return run(cmd).stdout


def samtools_total_mapped(bam: str, samtools: str = "samtools") -> int:
    """Total mapped reads in a BAM via `idxstats` (one fast call).

    Used as the library size for depth normalisation in the parental-DE stage.
    Robust for small candidate-gene panels where a median-of-ratios size factor
    (which needs many genes) would be unstable.
    """
    res = run([samtools, "idxstats", bam])
    total = 0
    for line in res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3:
            try:
                total += int(parts[2])     # column 3 = mapped read-segments
            except ValueError:
                continue
    return total


def ensure_bam_index(bam: str, samtools: str = "samtools") -> None:
    """Index a BAM if no index exists alongside it.

    Checks the two standard layouts (foo.bam.bai / foo.bam.csi and the
    splitext form foo.bai) explicitly rather than via string slicing.
    """
    stem, _ = os.path.splitext(bam)
    candidates = (bam + ".bai", bam + ".csi", stem + ".bai", stem + ".csi")
    if any(os.path.exists(c) for c in candidates):
        return
    run([samtools, "index", bam])


# ---------------------------------------------------------------------------
# Secondary alignment handling
# ---------------------------------------------------------------------------

def count_secondary_alignments(bam: str, samtools: str = "samtools") -> int:
    """Return the number of secondary alignments in a BAM (FLAG 0x100).

    Uses ``samtools view -c -f 0x100`` — a single pass over the BAM index for
    coordinate-sorted files. Fast even on large BAMs.
    """
    res = run([samtools, "view", "-c", "-f", "0x100", bam])
    out = res.stdout.strip()
    return int(out) if out else 0


def filter_secondary_alignments(bam: str, out_bam: str,
                                samtools: str = "samtools") -> str:
    """Write a new BAM excluding secondary alignments (FLAG 0x100); index it.

    Returns ``out_bam``. The output BAM is indexed immediately so it is ready
    for random-access queries without a separate ``ensure_bam_index`` call.
    """
    run([samtools, "view", "-F", "0x100", "-b", "-o", out_bam, bam])
    run([samtools, "index", out_bam])
    return out_bam


def prepare_bam(bam: str, *, filter_secondary: bool = False,
                samtools: str = "samtools") -> str:
    """Check for secondary alignments and optionally remove them.

    Returns the BAM path to use for downstream analysis — either ``bam``
    unchanged or a filtered copy with a ``.no_secondary.bam`` suffix written
    alongside the original.

    When ``filter_secondary`` is False (the default) secondary alignments are
    still excluded at runtime via the ``-F 3844`` flag in
    :func:`samtools_mpileup` / :func:`samtools_view`; this function only warns
    so the analyst is aware they exist.  Pass ``filter_secondary=True`` (CLI
    flag ``--filter-secondary``) to strip them from the BAM before analysis —
    useful for GATK or other tools that do not expose a per-flag exclude filter.

    The filtered copy is reused on repeated runs if it already exists alongside
    the original.
    """
    n = count_secondary_alignments(bam, samtools=samtools)
    if n == 0:
        return bam
    stem, _ = os.path.splitext(bam)
    filtered = stem + ".no_secondary.bam"
    if filter_secondary:
        if not os.path.exists(filtered):
            print(f"  [secondary] {n:,} secondary alignments in "
                  f"{os.path.basename(bam)} — writing filtered copy "
                  f"{os.path.basename(filtered)}", flush=True)
            filter_secondary_alignments(bam, filtered, samtools=samtools)
        else:
            print(f"  [secondary] {os.path.basename(bam)} has {n:,} secondary "
                  f"alignments — reusing {os.path.basename(filtered)}", flush=True)
        return filtered
    print(f"  WARNING: {os.path.basename(bam)} contains {n:,} secondary "
          f"alignments (excluded at runtime by -F flags; use "
          f"--filter-secondary to strip them from the BAM)", flush=True)
    return bam


# ---------------------------------------------------------------------------
# GATK ASEReadCounter (optional counting backend)
# ---------------------------------------------------------------------------

def gatk_ase_read_counter(bam: str, sites_vcf: str, reference: str,
                          output_tsv: str, *,
                          min_mapq: int, min_baseq: int,
                          min_depth: int = 1,
                          gatk: str = "gatk") -> str:
    """Run GATK ASEReadCounter for a single BAM against a sites VCF.

    Returns the path to the written TSV. Caller parses it (chr/position/
    refCount/altCount/...). Used only when --counter gatk is selected.
    """
    cmd = [gatk, "ASEReadCounter",
           "-R", reference, "-I", bam, "-V", sites_vcf,
           "-O", output_tsv,
           "--min-mapping-quality", str(min_mapq),
           "--min-base-quality", str(min_baseq),
           "--min-depth-of-non-filtered-base", str(min_depth)]
    run(cmd)
    return output_tsv


# ---------------------------------------------------------------------------
# WASP (optional reference-bias remap filtering)
# ---------------------------------------------------------------------------

def wasp_filter(bam: str, snp_dir: str, output_bam: str, *,
                wasp_dir: str, aligner_cmd: str,
                samtools: str = "samtools") -> str:
    """Placeholder wrapper for the WASP remap-and-filter workflow.

    WASP (van de Geijn et al. 2015) removes reads whose mapping changes when
    their overlapping SNP allele is swapped — the gold standard for removing
    reference-mapping bias. A full implementation chains:
        find_intersecting_snps.py -> re-align -> filter_remapped_reads.py
    which requires the project's chosen aligner. This stub documents the
    contract and raises until wired to your alignment command.
    """
    raise NotImplementedError(
        "WASP filtering is not yet wired. Provide --wasp-dir and the project "
        "aligner command, then implement the remap chain in external.wasp_filter. "
        "Until then use --bias-mode report or nmask."
    )
