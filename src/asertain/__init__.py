"""ASErtain — allele-specific expression discovery & validation for F1 hybrids.

A staged toolkit for detecting cis-regulatory divergence from allele-specific
expression (ASE) in F1 hybrids, designed around the Penstemon kunthii ×
amphorellae cross but generic to any two-parent / hybrid design.

Pipeline stages (exposed as `asertain <subcommand>`):

    diagnose   identify diagnostic SNPs from exact-parent genotypes
    count      count allele-specific reads in F1 BAMs
    test       replicate-aware ASE statistics, gene-level calls
    contrast   cis/trans decomposition against parental DE
    report     plots + HTML summary
    run        orchestrate the whole pipeline from one config file

See DESIGN.md for the scientific rationale and file-format contracts.
"""

__version__ = "0.1.0"
__author__ = "Luis Javier Madrigal-Roca"
