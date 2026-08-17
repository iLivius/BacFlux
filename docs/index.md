# Overview

BacFlux is a [Snakemake](https://snakemake.github.io/) workflow for bacterial
whole-genome sequencing. It carries one isolate from raw reads — or from a genome you
already assembled — through assembly, decontamination, quality control, taxonomic
placement, functional annotation, and screening for resistance genes, virulence
factors, plasmids and prophages. It can also answer the question that usually comes
next: for each resistance gene found, is it sitting inside something that can move to
another bacterium? A single key in the config file, `mode`, decides which of four
front ends runs. This site documents version **v2.0.0**.

BacFlux belongs to the [BioFlux](https://github.com/stars/iLivius/lists/bioflux)
family of pipelines.

## The four modes

|                          | `illumina` | `nanopore` | `hybrid` | `contigs` |
|--------------------------|------------|------------|----------|-----------|
| **Input**                | Illumina paired-end reads | long reads from Oxford Nanopore Technologies (ONT) | Illumina and ONT reads from the same isolate | a finished assembly (FASTA) |
| **Read QC**              | PhiX removal (`bowtie2`), then `fastp` | `NanoPlot` either side of `Filtlong` | both legs | none |
| **Assembly**             | `SPAdes` | `Flye`, reoriented by `dnaapler`, polished by `Medaka`\* | `Flye` as in `nanopore`, then corrected with the short reads by `Polypolish`; a `SPAdes` draft is kept as the comparator | no assembler |
| **Contig filter**        | ≥500 bp and ≥2× coverage, read off SPAdes' own headers | none | the same cutoffs, on the SPAdes draft only | the same cutoffs when the headers are SPAdes-style; otherwise headers are trimmed and nothing is removed |
| **Delivered genome**     | the screened SPAdes contigs | the polished Flye assembly | the Polypolish-corrected ONT assembly | the screened input contigs |
| **Needs reads later on** | yes — the CARD read leg, and IS copy number if the mobilome module is on | no | yes, as `illumina` | no |

\* optional: `parameters.{nanopore,hybrid}.medaka_model: false` skips the polishing
step. `hybrid` is the expensive mode — it assembles twice, and `Snippy` then counts how
far the two assemblies still differ.

Each front end ends by writing one file,
`02.assembly/{sample}/contigs_final.fasta`. Two shared stages sit either side of it:
decontamination, which *is* the step that writes that file in `illumina` and
`contigs` mode and runs just before the final polish in the two long-read modes; and
assembly QC, which scores whatever came out (in `hybrid`, both genomes). From
taxonomy onwards every rule reads `contigs_final.fasta` and nothing else — taxonomy,
annotation, AMR, plasmids, prophages, the optional mobilome stage and the report
exist in one copy each, so a result means the same thing whichever way the isolate
was sequenced. See [Modes](modes/index.md) for the front ends and
[Choosing a mode](getting-started/choosing-a-mode.md) for which one a dataset belongs
in.

## One command, one config, one report

All four modes are launched the same way, from the repository root:

```bash
snakemake --sdm conda --cores 16 --configfile config/config.yaml
```

`--sdm conda` (short for `--software-deployment-method`) tells Snakemake to build
each rule's own conda environment on first use, so only Snakemake itself has to be
installed by hand. `--cores` is the single CPU knob: every CPU-bound rule declares
its own `threads:`, so `--cores` already sets the ceiling and `--jobs` must not be
added alongside it. The reasoning is on [Running BacFlux](reference/running.md).

The config file is required on the command line; BacFlux deliberately ships no
default. Mode, input and output paths, database locations and every tunable
parameter live in the file you pass, and a missing key stops the run at parse time —
in seconds, before any job starts — rather than quietly falling back to a placeholder
path. Every key is documented in [Configuration](reference/configuration.md).

Every mode ends with one [MultiQC](https://multiqc.info/) page,
`09.report/multiqc_report.html`, so a finished run can be read end to end from a
single HTML file.

!!! note "What has to be prepared in advance"

    BacFlux installs its own software, but five reference databases are too large or
    too version-sensitive to fetch automatically and must be on disk before the first
    run: Bakta (database v6.0), NCBI core nt, eggNOG, GTDB-Tk (GTDB R232) and Platon.
    Several others are downloaded for you on first use, or can be pointed at a copy
    you already hold. See [Reference databases](getting-started/databases.md).

## Where AMR mobility fits

Stage `08.mobilome` is optional and off by default (`mobilome.run: false`). Turned
on, it runs AMRFinderPlus over the annotated genome and, for every resistance gene it
calls, reports the mobile-element context around that gene and a mobility tier, on a
ladder from *chromosomal, no mobile-element context — intrinsic candidate* up to
*inside an ICE or on a conjugative plasmid — predicted self-transmissible*. That
distinction, intrinsic versus acquired, is the one regulators ask about, and it is
not something a gene list alone can answer.

Every call is a prediction and carries its own confidence tier and contig-edge flags,
because on a fragmented assembly the structure you are looking for is often the thing
that broke the assembly. Start at [the mobilome module](mobilome/index.md), and read
[Draft assemblies](mobilome/draft-assemblies.md) before drawing conclusions from a
short-read run.

## Where to go next

| Page                                                  | Read it for                                                                       |
|-------------------------------------------------------|-----------------------------------------------------------------------------------|
| [Installation](getting-started/installation.md)       | cloning the repository and the Snakemake launcher environment                     |
| [Reference databases](getting-started/databases.md)   | what to download, what BacFlux fetches, and what you can point at a copy you hold |
| [Quick start](getting-started/quick-start.md)         | the shortest path from a clone to a finished run                                  |
| [Choosing a mode](getting-started/choosing-a-mode.md) | which mode a given dataset belongs in, and what a wrong one does                  |
| [Coming from v1](getting-started/from-v1.md)          | the `mode:` key, the config that has to be rebuilt, and the renumbered output     |
| [Modes](modes/index.md)                               | the four front ends, one page each                                                |
| [Antimicrobial resistance](analysis/amr.md)           | the ABRicate and CARD legs, and what each is and is not evidence of               |
| [The mobilome module](mobilome/index.md)              | AMR mobility: the ladder, the switch, the output, and its limits                  |
| [Configuration](reference/configuration.md)           | every config key, which modes read it, its default and its trade-off              |
| [Output files](reference/output.md)                   | what each numbered stage holds and which files are the deliverables               |
| [Running BacFlux](reference/running.md)               | dry runs, CPU budget, restarting a run                                            |
| [Troubleshooting](troubleshooting.md)                 | the failures that actually happen, and how to get past them                       |
| [Rationale](about/rationale.md)                       | why one workflow with four entry points                                           |
| [Citation and references](about/citation.md)          | how to cite BacFlux and the tools that do the work                                |

BacFlux is released under the MIT License; the tools and databases it calls carry
their own terms, listed on [Licensing](about/licensing.md).
