# Modes

BacFlux has four front ends and one analysis. `mode` in the config picks the front
end; nothing on the command line changes with it.

A front end is everything up to and including the assembly: read QC, the assembler,
the contig filter and the contamination screen. It ends by writing one file,
`02.assembly/{sample}/contigs_final.fasta`. From `03.taxonomy` onwards every rule
reads that file and nothing else about the assembly, so the whole downstream half
exists in one copy rather than four.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **ONT** | Oxford Nanopore Technologies — the long-read sequencing platform |
    | **QC** | quality control |
    | **AMR** | antimicrobial resistance |
    | **IS** | insertion sequence — the smallest mobile element, carrying only the genes it needs to move |

| | `illumina` | `nanopore` | `hybrid` | `contigs` |
|---|---|---|---|---|
| **Input** | paired-end reads | ONT long reads | both, same isolate | a FASTA |
| **Read QC** | PhiX removal (Bowtie2), fastp | NanoPlot, Filtlong | both | — |
| **Assembler** | SPAdes `--isolate` | Flye | Flye, plus SPAdes for the comparator | — |
| **After assembly** | length and coverage filter | dnaapler, then Medaka\* | dnaapler, Medaka\*, then Polypolish | header-aware filter |
| **Screened** | the SPAdes draft | the reoriented Flye assembly | the Illumina draft; the ONT leg is screened at read level | the supplied contigs |
| **Delivered genome** | the screened contigs | the polished assembly | the Polypolish-corrected ONT assembly | the screened contigs |
| **Front-end rules** | 6 | 7 | 16 | 1 |

\* optional: `parameters.{nanopore,hybrid}.medaka_model: false` skips the polishing
step, and the rule is then not defined at all.

## What is identical in all four

Taxonomy, annotation, AMR and virulence screening, plasmid and prophage detection,
the optional mobilome module and the MultiQC report are one set of rules, run over
`contigs_final.fasta`. So are the contamination screen and assembly QC — those two
sit inside the front end, at the point each mode needs them, but they are the same
code everywhere. Only the mapper that builds the coverage track differs: Bowtie2
where there are short reads, minimap2 otherwise.

That is what makes a cross-technology comparison meaningful: when two answers come
out of two sets of rules, every difference between them is confounded by the code.
Out of one set, the difference is the data. See [Rationale](../about/rationale.md).

## What the mode still decides downstream

A few shared steps need something a mode does not have. They are gated on the
capability rather than on the mode name — `HAS_SHORT_READS`, `HAS_LONG_READS`,
`HAS_READS` in `workflow/rules/shared/00_common.smk` — and when the capability is
absent the rule is never defined, so nothing waits for it.

| Shared step | Needs | `illumina` | `nanopore` | `hybrid` | `contigs` |
|---|---|:-:|:-:|:-:|:-:|
| CARD read mapping (`05.amr/mapping`) | short reads | ✓ | — | ✓ | — |
| Qualimap mapping evaluation | any reads | ✓ | ✓ | ✓ | — |
| Coverage that carries information for BlobTools | any reads | ✓ | ✓ | ✓ | — |
| Circular-replicon table handed to Bakta | Flye's topology calls | — | ✓ | ✓ | — |
| Second BLAST of the final assembly, for the plasmid stage | long reads | — | ✓ | ✓ | — |
| IS copy number from reads (mobilome, optional) | short reads | ✓ | — | ✓ | — |
| Two genomes through QC and taxonomy, plus Snippy | both technologies | — | — | ✓ | — |

## The four front ends

| Page | What it covers |
|---|---|
| [Illumina](illumina.md) | PhiX removal, fastp, SPAdes, and the length and coverage cutoffs |
| [Nanopore](nanopore.md) | NanoPlot, Filtlong, Flye, dnaapler, Medaka — and the filter settings that decide whether a small plasmid survives |
| [Hybrid](hybrid.md) | how the Illumina leg becomes the ONT leg's contamination filter, Polypolish, and the Snippy comparison |
| [Pre-assembled contigs](contigs.md) | what runs when the input is a finished genome, and why the coverage filter is conditional |

[Choosing a mode](../getting-started/choosing-a-mode.md) answers the other question:
which one a given dataset belongs in.
