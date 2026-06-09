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

def samtools_mpileup(bam: str, region: str, *,
                     min_mapq: int, min_baseq: int,
                     reference: Optional[str] = None,
                     samtools: str = "samtools") -> str:
    """Return the raw mpileup stdout for a single region (e.g. 'chr1:100-100')."""
    cmd = [samtools, "mpileup", "-r", region,
           "-q", str(min_mapq), "-Q", str(min_baseq),
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
                  exclude_flags: int = 3844,
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
