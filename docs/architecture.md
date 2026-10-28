# Workflow architecture

BacFlux supports four input modes within one Snakemake workflow. The mode selects
input processing and assembly; shared rules then analyse the final assembly.
Read-dependent analyses run only when the required inputs are available.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **ONT** | Oxford Nanopore Technologies — the long-read sequencing platform |
    | **AMR** | antimicrobial resistance |
    | **EFSA** | European Food Safety Authority, whose thresholds the ABRicate analysis applies |

## The shape of a run

Two diagrams show the shape of a run: the first what changes with the mode, the second
the analyses every mode shares. They are schematics rather than the execution DAG: the
position of decontamination and the availability of individual rules depend on the
selected mode. For the mode-specific sequence, see [Modes](modes/index.md) and
[Decontamination](analysis/decontamination.md).

**Part 1 — what changes with the mode you chose.**

```mermaid
---
config:
  themeVariables:
    fontSize: 14px
  flowchart:
    nodeSpacing: 35
    rankSpacing: 35
---
flowchart LR
    I[illumina<br/>paired short reads] --> IF[PhiX removal, fastp<br/>SPAdes assembly]
    N[nanopore<br/>ONT long reads] --> NF[Filtlong, Flye,<br/>dnaapler, Medaka]
    H[hybrid<br/>short and long reads] --> HF[both front ends,<br/>Polypolish, Snippy]
    C[contigs<br/>an existing assembly] --> CF[length filter only]
    IF --> DRAFT([a draft<br/>assembly])
    NF --> DRAFT
    HF --> DRAFT
    CF --> DRAFT
```

**Part 2 — the same whichever mode ran.**

```mermaid
---
config:
  themeVariables:
    fontSize: 14px
  flowchart:
    nodeSpacing: 35
    rankSpacing: 35
---
flowchart LR
    DECON[02.assembly<br/>decontamination<br/>of the draft] --> FINAL([contigs_final.fasta])
    FINAL --> QC[02.assembly/eval<br/>QUAST, CheckM,<br/>Qualimap]
    FINAL --> TAX[03.taxonomy<br/>GTDB-Tk]
    FINAL --> ANN[04.annotation<br/>Bakta, eggNOG,<br/>antiSMASH, dbCAN]
    FINAL --> AMR[05.amr<br/>ABRicate,<br/>CARD read mapping]
    FINAL --> PLA[06.plasmids<br/>Platon,<br/>optionally geNomad]
    FINAL --> PHA[07.phages<br/>VirSorter2 or<br/>geNomad, CheckV]
    ANN --> MOB[08.mobilome<br/>optional<br/>ISEScan, CONJScan,<br/>AMRFinderPlus<br/>→ mobility tiers]
    TAX --> MOB
    PLA --> MOB
    QC --> REP[09.report<br/>MultiQC]
    TAX --> REP
    ANN --> REP
```

The numbered output stages are shared across modes, but their individual rules
can be conditional. `08.mobilome` and its optional layers are disabled by default.
Medaka polishing runs by default in the long-read modes; `medaka_model: false` skips
it. Read QC and CARD read mapping require the corresponding reads.

!!! note "Why the entry points converge so early"

    The modes differ in what arrives — short reads, long reads, both, or contigs.
    Once a final assembly exists, the downstream rules are the same in every mode,
    and steps that need reads run only where those reads exist. Shared
    implementation also removes the maintenance differences between the earlier
    separate workflows.

    The mobilome module applies the same contig-edge flags and low-confidence caps
    in every mode; a fragmented assembly, most often a short-read one, simply
    triggers them more often — see [the mobilome overview](mobilome/index.md).

## What runs where

| stage | runs in | optional? |
|---|---|---|
| `01.reads` | every mode that has reads | — |
| `02.assembly` | all four | — |
| `03.taxonomy` | all four | — |
| `04.annotation` | all four | — |
| `05.amr` | all four; the CARD read analysis needs short reads | — |
| `06.plasmids` | all four | Platon always; the geNomad second opinion is opt-in |
| `07.phages` | all four | always runs; VirSorter2 by default, geNomad optional |
| `08.mobilome` | all four | **off by default** |
| `09.report` | all four | — |

## Complementary AMR analyses { #the-three-amr-legs-and-why-there-are-three }

The three analyses use different inputs and detection methods. ABRicate screens
the contigs, BBMap maps short reads to CARD, and AMRFinderPlus provides curated
calls for the optional mobilome analysis. Their outputs are complementary and
have different interpretation limits.

```mermaid
---
config:
  themeVariables:
    fontSize: 14px
  flowchart:
    nodeSpacing: 35
    rankSpacing: 35
---
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
├── 01.reads/       read QC reports and the PhiX genome and index (kept);
│                   trimmed, PhiX-free, filtlong-filtered reads   (temporary)
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
