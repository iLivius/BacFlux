# Workflow architecture

BacFlux is one workflow with four ways in. Whichever way you enter, you leave through
the same nine numbered stages, and the `mode` key in the config decides nothing except
how a set of contigs gets produced in the first place.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **ONT** | Oxford Nanopore Technologies — the long-read sequencing platform |
    | **AMR** | antimicrobial resistance |
    | **EFSA** | European Food Safety Authority, whose thresholds the ABRicate leg applies |

## The shape of a run

The diagram is in two labelled halves. Everything in the first box changes with the
mode you chose; everything in the second box is the same whichever mode ran. That is
why there is one workflow instead of four.

```mermaid
flowchart TD
    subgraph P1["PART 1 &nbsp; what changes with the mode you chose"]
        direction TB
        I[illumina<br/>paired short reads]
        N[nanopore<br/>ONT long reads]
        H[hybrid<br/>short reads and long reads]
        C[contigs<br/>an assembly you already have]

        IF[fastp trimming, PhiX removal<br/>assembly with SPAdes]
        NF[Filtlong, assembly with Flye<br/>Medaka polishing, dnaapler]
        HF[both front ends, then Polypolish<br/>and a Snippy comparison]
        CF[length filter only]

        I --> IF
        N --> NF
        H --> HF
        C --> CF

        DRAFT([a draft assembly])
        IF --> DRAFT
        NF --> DRAFT
        HF --> DRAFT
        CF --> DRAFT
    end

    subgraph P2["PART 2 &nbsp; the same whichever mode ran"]
        direction TB
        DECON[02.assembly<br/>decontamination: BLAST, BlobTools, genus selector]
        FINAL([contigs_final.fasta<br/>every stage below reads this one file])
        DECON --> FINAL

        TAX[03.taxonomy<br/>GTDB-Tk]
        ANN[04.annotation<br/>Bakta, eggNOG, antiSMASH, dbCAN]
        AMR[05.amr<br/>ABRicate, and BBMap against CARD]
        PLA[06.plasmids<br/>Platon, optionally geNomad]
        PHA[07.phages<br/>VirSorter2 or geNomad, then CheckV]
        REP[09.report<br/>MultiQC]

        FINAL --> TAX
        FINAL --> ANN
        FINAL --> AMR
        FINAL --> PLA
        FINAL --> PHA

        MOB[08.mobilome<br/>OFF unless you switch it on<br/>ISEScan, CONJScan, AMRFinderPlus]
        LADDER([mobility table<br/>one row per AMR gene, tier 1 to 6])

        AMR --> MOB
        PLA --> MOB
        MOB --> LADDER

        TAX --> REP
        ANN --> REP
        PHA --> REP
    end

    DRAFT --> DECON
```

Only `08.mobilome` is optional. Everything else runs on every sample, in every mode.

!!! note "Why the entry points converge so early"

    The four modes exist only because a genome can arrive as reads or as contigs, and
    reads can be short, long or both. Once there is a decontaminated assembly, nothing
    downstream cares how it was made — which is why Part 1 is short and Part 2 is
    everything else. The four v1 workflows shared roughly two-thirds of their code and
    had to be kept in step by hand; here that two-thirds exists once.

    One thing does stay mode-aware: a short-read assembly is more fragmented, so every
    mobilome call carries contig-edge flags and anything spanning contigs is capped at
    low confidence — see [the mobilome overview](mobilome/index.md).

## What runs where

| stage | runs in | optional? |
|---|---|---|
| `01.reads` | every mode that has reads | — |
| `02.assembly` | all four | — |
| `03.taxonomy` | all four | — |
| `04.annotation` | all four | — |
| `05.amr` | all four; the CARD read leg needs short reads | — |
| `06.plasmids` | all four | Platon always; the geNomad second opinion is opt-in |
| `07.phages` | all four | always runs; VirSorter2 by default, geNomad optional |
| `08.mobilome` | all four | **off by default** |
| `09.report` | all four | — |

## The three AMR legs, and why there are three

The one place BacFlux deliberately does the same job three times. Each leg fails in a
different way, so the three together say more than any one of them.

```mermaid
flowchart LR
    R([trimmed reads])
    F([contigs_final.fasta])

    A[ABRicate on the contigs<br/>8 databases, EFSA 80% and 70%]
    B[BBMap, reads against CARD<br/>at least 70% of the gene covered]
    P[AMRFinderPlus on the proteins<br/>curated cutoffs, one per gene]

    F --> A
    R --> B
    F --> P

    A --> AW[blind spot: a gene the<br/>assembler collapsed]
    B --> BW[blind spot: no gene boundaries,<br/>and it needs short reads]
    P --> PW[blind spot: only with the<br/>mobilome module switched on]

    A --> S([AMR summary matrix])
    B --> CR([CARD report, grouped by<br/>CARD's own resistance mechanism])
    P --> M([mobility ladder, tier 1 to 6])
```

Read the detail on [the AMR page](analysis/amr.md).

## Where the files land

```text
output_dir/
├── 01.reads/       trimmed, PhiX-free, filtlong-filtered reads   (temporary)
├── 02.assembly/    draft and decontaminated contigs, CheckM, QUAST
├── 03.taxonomy/    GTDB-Tk classification
├── 04.annotation/  Bakta, eggNOG, antiSMASH, dbCAN
├── 05.amr/         ABRicate tables, CARD read mapping
├── 06.plasmids/    Platon calls, optional geNomad concordance
├── 07.phages/      prophage predictions, CheckV quality
├── 08.mobilome/    mobility table                    (only when switched on)
├── 09.report/      MultiQC
└── logs/           one log per rule, per sample
```

Every path and every file is listed on [the output reference page](reference/output.md).
