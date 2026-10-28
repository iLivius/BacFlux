# Overview

BacFlux is a [Snakemake](https://snakemake.github.io/) workflow for bacterial
isolate genome analysis. It accepts Illumina reads, Oxford Nanopore reads, both read
types from the same isolate, or an existing assembly. The workflow combines
assembly processing, decontamination, quality control, taxonomic classification,
functional annotation, and screening for resistance genes, virulence factors,
plasmids and prophages.

BacFlux was developed to support the genomic characterisation of bacterial strains
considered for biocontrol, biostimulant and biofertiliser applications. Its annotation
and screening results help assess candidate strains, including the presence of
resistance-associated and virulence-associated genes. These predictions support
strain selection alongside phenotypic evidence and application-specific assessment.

An optional mobilome module assesses the genomic context and potential mobility
of detected antimicrobial resistance genes. The configuration key `mode` selects
the input mode. This site documents version **v2.1.0**.

BacFlux belongs to the [BioFlux](https://github.com/stars/iLivius/lists/bioflux)
family of pipelines.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **IS** | insertion sequence — the smallest mobile element, carrying only the genes it needs to move |
    | **AMR** | antimicrobial resistance |
    | **GTDB** | Genome Taxonomy Database, the reference used for taxonomic placement |
    | **ICE** | integrative and conjugative element — integrates into the chromosome and carries its own conjugation machinery |
    | **IME** | integrative mobilisable element — integrates like an ICE but has no conjugation machinery of its own, so it needs a helper element to transfer |
    | **MGE** | mobile genetic element — DNA that can move within a genome or between bacteria: insertion sequences, transposons, integrons, plasmids, ICEs |

## The workflow at a glance

[![Metro map of BacFlux v2.1.0. Blue (illumina), purple (hybrid), green (nanopore) and yellow (contigs) lines run through read QC, assembly, decontamination and polishing to contigs_final.fasta, then through assembly QC, taxonomy, annotation, antimicrobial resistance, plasmids and prophages to the MultiQC report. A dashed red line, the optional mobilome module, reads the annotation, the taxonomy and the plasmid calls and ends in the mobility-ladder table.](assets/bacflux_metro.png#only-light)](assets/bacflux_metro.png)
[![Metro map of BacFlux v2.1.0. Blue (illumina), purple (hybrid), green (nanopore) and yellow (contigs) lines run through read QC, assembly, decontamination and polishing to contigs_final.fasta, then through assembly QC, taxonomy, annotation, antimicrobial resistance, plasmids and prophages to the MultiQC report. A dashed red line, the optional mobilome module, reads the annotation, the taxonomy and the plasmid calls and ends in the mobility-ladder table.](assets/bacflux_metro_dark.png#only-dark)](assets/bacflux_metro_dark.png)

Each coloured line is one mode. The four lines take their own routes to
`contigs_final.fasta`, the genome each mode delivers, and from there run through the
same stages. In `hybrid` the purple line has two arms: the Illumina arm assembles a
SPAdes draft and decontaminates it, and the Illumina reads that map back to the kept
contigs are what Polypolish uses to correct the long-read assembly. An open ring marks
a step the config switches on or off; the prophage stage runs VirSorter2 by default,
or geNomad instead. The dashed red line is the optional mobilome module. Select the
map to open it at full size.

## The four modes

|                          | `illumina` | `nanopore` | `hybrid` | `contigs` |
|--------------------------|------------|------------|----------|-----------|
| **Input**                | Illumina paired-end reads | long reads from Oxford Nanopore Technologies (ONT) | Illumina and ONT reads from the same isolate | a finished assembly (FASTA) |
| **Read QC**              | PhiX removal (`bowtie2`), then `fastp` | `NanoPlot` before and after `Filtlong` | both of the above | none |
| **Assembly**             | `SPAdes` | `Flye`, reoriented by `dnaapler`, polished by `Medaka`\* | `Flye` as in `nanopore`, then corrected with the short reads by `Polypolish`; a `SPAdes` draft is kept as the comparator | no assembler |
| **Contig filter**        | ≥500 bp and ≥2× coverage, read off SPAdes' own headers | none | the same cutoffs, on the SPAdes draft only | the same cutoffs when the headers are SPAdes-style; otherwise headers are trimmed and nothing is removed |
| **Delivered genome**     | the screened SPAdes contigs | the polished Flye assembly | the Polypolish-corrected ONT assembly | the screened input contigs |
| **Needs reads later on** | yes — the CARD read analysis, and IS copy number when the mobilome module and its ISOSDB layer are on | no | yes, as `illumina` | no |

\* Medaka polishing is optional. Set
`parameters.{nanopore,hybrid}.medaka_model: false` to skip it. Hybrid mode requires
more computation because it assembles both read types and compares the assemblies
with Snippy.

Every mode produces `02.assembly/{sample}/contigs_final.fasta` for the shared analysis
stages. The primary decontamination screen filters the draft assembly; in hybrid mode
this is the Illumina draft, while a second screen of the final long-read assembly writes
a verdict per contig and removes nothing. Assembly QC evaluates the final genome and, in
hybrid mode, the Illumina comparator. See
[Decontamination](analysis/decontamination.md).

The shared downstream stages use the final assembly together with the annotation,
taxonomy, reads or other inputs required by each analysis. This common processing
improves consistency across datasets. Differences in read quality and assembly
fragmentation still affect detection and interpretation. See [Modes](modes/index.md)
and [Choosing a mode](getting-started/choosing-a-mode.md).

## One command, one config, one report

All four modes are launched the same way, from the repository root:

```bash
snakemake --sdm conda --cores 16 --configfile config/config.yaml
```

`--sdm conda` enables the conda environments specified by the workflow. Install
Snakemake first; the analysis environments are then created as required.
`--cores` sets the CPU budget for local execution. Use it without a separate
`--jobs` setting for a local run. See [Running BacFlux](reference/running.md).

Pass the configuration file explicitly with `--configfile`. BacFlux supplies a template,
but does not load it automatically. Use a copy to set the mode, input and output paths,
database locations and analysis parameters. Required missing keys stop planning; an
existing key with an invalid or empty path may fail later. See
[Configuration](reference/configuration.md) and
[Quick start](getting-started/quick-start.md).

Every mode produces a [MultiQC](https://multiqc.info/) overview at
`09.report/multiqc_report.html`. Use it for the batch-level QC summary and read the
stage-specific tables and reports for annotation, AMR and mobilome interpretation.

!!! note "What has to be prepared in advance"

    BacFlux installs its own software, but five reference databases are too large or
    too version-sensitive to fetch automatically and must be on disk before the first
    run: Bakta (database v6.0), NCBI core nt, eggNOG, GTDB-Tk (GTDB R232) and Platon.
    Several others are downloaded for you on first use, or can be pointed at a copy
    you already hold. See [Reference databases](getting-started/databases.md).

## Where AMR mobility fits

Stage `08.mobilome` is disabled by default (`mobilome.run: false`). When enabled,
it uses AMRFinderPlus and mobile-element evidence to assign a context and mobility
tier to each assessed resistance gene. The reported categories range from
`intrinsic_candidate` to `predicted_self_transmissible`.

These categories describe genomic evidence and predicted transfer potential: the
evidence for the intrinsic-versus-acquired question a strain assessment asks of each
resistance gene, which a gene list alone cannot answer. They do not demonstrate
transfer or establish that resistance is intrinsic or acquired.
Interpret the tiers together with confidence fields, contig-edge flags and assembly
quality. Start with [the mobilome module](mobilome/index.md) and
[Draft assemblies](mobilome/draft-assemblies.md).

## Where to go next

| Page                                                  | Read it for                                                                       |
|-------------------------------------------------------|-----------------------------------------------------------------------------------|
| [Installation](getting-started/installation.md)       | cloning the repository and the Snakemake launcher environment                     |
| [Reference databases](getting-started/databases.md)   | what to download, what BacFlux fetches, and what you can point at a copy you hold |
| [Quick start](getting-started/quick-start.md)         | the shortest path from a clone to a finished run                                  |
| [Choosing a mode](getting-started/choosing-a-mode.md) | which mode a given dataset belongs in, and what a wrong one does                  |
| [Coming from v1](getting-started/from-v1.md)          | the `mode:` key, the config that has to be rebuilt, and the renumbered output     |
| [Modes](modes/index.md)                               | input processing and assembly for each mode                                                |
| [Antimicrobial resistance](analysis/amr.md)           | ABRicate screening, CARD read mapping and interpretation limits               |
| [The mobilome module](mobilome/index.md)              | mobility tiers, activation, outputs and interpretation limits                  |
| [Configuration](reference/configuration.md)           | every config key, which modes read it, its default and its trade-off              |
| [Output files](reference/output.md)                   | what each numbered stage holds and which files are the deliverables               |
| [Running BacFlux](reference/running.md)               | dry runs, CPU budget, restarting a run                                            |
| [Troubleshooting](troubleshooting.md)                 | common failures, and how to get past them                                         |
| [Rationale](about/rationale.md)                       | why one workflow with four entry points                                           |
| [Citation and references](about/citation.md)          | how to cite BacFlux and the tools that do the work                                |

BacFlux is released under the MIT License; the tools and databases it calls carry
their own terms, listed on [Licensing](about/licensing.md).
