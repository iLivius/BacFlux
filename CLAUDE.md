# BacFlux — Claude Code project notes

Snakemake workflow for bacterial short-read WGS. MIT licensed, archived on Zenodo.
Run: `snakemake --sdm conda`. Conda-per-rule is the dependency model —
every new tool must be bioconda-installable.

## Current work
Mobilome / AMR-mobility module. Full design spec: @docs/mobilome_module_SPEC.md
Read §2 (conceptual model), §3.3 (rejected tools, with reasons) and §11
(licensing constraints) before proposing any change to this module.

## Hard rules
- Never vendor licence-encumbered databases; config URL + README notice only.
- Never copy CC BY-NC-SA code into this repo (see spec §11).
- Every filtering decision gets an audit TSV with a reason column.
