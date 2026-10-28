# Rationale

Characterising a microbial strain requires several linked analyses: assembly and
quality assessment, taxonomic identification, functional annotation, and screening
for resistance-associated genes, virulence factors and mobile elements. BacFlux was
developed to support this work for candidate strains in biocontrol, biostimulant and
biofertiliser research. Linking the tools in one workflow reduces manual processing
and makes the analysis easier to repeat and inspect.

BacFlux writes these analyses down once, as a [Snakemake](https://snakemake.github.io/)
workflow. Every step declares its inputs and its outputs, and every step that calls an
external tool declares its own conda environment, so a run can be resumed after an
interruption, repeated on another machine, or shared with a colleague. Reproduction also
requires the same inputs, software and reference data; the configuration alone is not
sufficient.

BacFlux v2 consolidates four related workflows while retaining their input modes.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **AMR** | antimicrobial resistance |

## Four ways in, one analysis

BacFlux began as one workflow, for Illumina paired-end reads. Three more followed, each
because a project needed it — FastaFlux for genomes that arrived already assembled,
BacFluxL when Oxford Nanopore Technologies (ONT) sequencing became routine, and
BacFluxL+ for isolates sequenced both ways. FastaFlux was a second Snakefile alongside
the first; the two long-read ones had repositories of their own. All four were released,
and work published with them stays citable — BacFluxL and BacFluxL+ are retired in
favour of this repository and their DOIs stay valid. See
[Coming from v1](../getting-started/from-v1.md).

They also had almost everything in common. The downstream analyses share most of their
processing once a final assembly exists, although available reads, assembly quality and
mode-specific metadata still matter. Taxonomic placement, annotation, AMR and virulence
screening, plasmid and prophage detection are the same tools, in the same order, at the
same thresholds, whatever went in at the top; so is the decontamination screen, apart
from the mapper that builds its coverage track — bowtie2 for short reads, minimap2 for
long ones. Keeping that in four Snakefiles meant four places to make one improvement,
and four places to look when a result seemed odd.

Version 2.0.0 keeps all four ways in and gives them one analysis. One key in the config
file, `mode`, chooses the front end; nothing else on the command line changes.

| `mode:` | Input | Front end | Front-end rules |
|---|---|---|--:|
| `illumina` | Illumina paired-end reads | PhiX removal, fastp, SPAdes | 6 |
| `nanopore` | Oxford Nanopore reads | NanoPlot, Filtlong, Flye, dnaapler, Medaka | 7 |
| `hybrid` | both | both front ends, then Polypolish correction of the ONT assembly with the short reads, and a Snippy comparison of the two | 16 |
| `contigs` | a finished assembly (FASTA) | contig filtering only | 1 |

Those four front ends are 30 rules between them. The analysis they feed is 62 rules in
`workflow/rules/shared/` — 22 of them the optional mobilome module — written once
instead of four times, plus six `_local` variants that link a database copy you
already hold. [Choosing a mode](../getting-started/choosing-a-mode.md) covers
which one a given dataset belongs in.

## Why the shared rules matter { #what-writing-a-step-once-buys }

**Shared processing reduces one source of variation.** The same annotation, taxonomy
and screening rules are used across modes wherever the required inputs are available.
This reduces differences caused by maintaining separate copies of the analysis code.
It does not remove effects of sequencing technology, filtering or assembly quality.
In hybrid mode, QUAST, CheckM and GTDB-Tk assess the final assembly and the Illumina
comparator, while Snippy compares four stages of the ONT assembly with the Illumina
assembly.

**A fix reaches every mode at once.** A re-run in July found that the CheckV database
ships without a DIAMOND index, so the index in any shared copy was built by whichever
DIAMOND version that site happened to have, and the pinned one could not read it. The
fix — build a local view of the read-only database from symlinks plus an index this
workflow builds itself — was written in one rule, and all four modes had it. Two of the
three bugs found that week stayed hidden in the `illumina` and `nanopore` test runs: a
sample with no viral contigs takes the empty-input shortcut and never touches the
database. Only the `hybrid` run went down that path — and fixing it there fixed it
everywhere.

**A new capability is built once.** The mobilome module is 22 rules and eleven Python
helpers. Across the four it would have been copied four times and maintained in
parallel. That is what settled the order of work: merge first, then write the module.
See [the mobilome module](../mobilome/index.md).

## What the four modes actually share { #how-the-four-modes-actually-share }

**A common final assembly.** Each mode writes
`02.assembly/{sample}/contigs_final.fasta`. Shared stages use this assembly together
with the annotation, taxonomy, reads and metadata they require. Long-read topology
information is supplied to Bakta separately. Hybrid mode also retains an Illumina
assembly for comparison.

**Fixed stage numbers.** All technology-specific work lives under `01.reads` and
`02.assembly`, so every shared stage lands on the same number in every mode: taxonomy is
always `03.taxonomy`, AMR always `05.amr`. Previously each workflow had a different
number of front-end stages, so the shared stages were numbered differently in each — the
single largest reason "the same" rule was not literally the same file. See
[Output files](../reference/output.md).

**Capability flags rather than mode names.** Where a shared step genuinely depends on the
technology, it asks what the run *has*, not which mode it *is*: `HAS_SHORT_READS`,
`HAS_LONG_READS`, `HAS_READS`, defined in `workflow/rules/shared/00_common.smk`. The CARD
read-mapping analysis is gated on `HAS_SHORT_READS` because it needs reads — not because of a
hard-coded list of mode names that a fifth mode would have to be added to.

**One environment set.** 33 files in `workflow/envs/`: one per tool, with the version
pinned in the file, plus `download.yaml` (wget, tar, gzip and bzip2 for the download
rules) and `python.yaml` (Python for the helper scripts). A few tool environments carry
a second pinned package the tool needs at run time (`diamond` alongside Bakta,
`samtools` alongside bowtie2). The merge turned up
exactly one version drift across the four, fastp 1.0.1 against 1.1.0, which says the
copies had been kept in step by hand — at a cost that does not show up anywhere in the
results.

!!! note "Where sharing stops"

    `hybrid` runs both technologies, so it keeps its own copies of the Illumina rules and
    of Flye and dnaapler, each marked *keep in sync* where it sits. Its Medaka rule is a
    copy too, but a deliberately different one, and the comment above it lists the three
    places it departs from `nanopore`'s. Meanwhile two long-read rules that would
    otherwise be byte-identical in `nanopore` and `hybrid` sit in `rules/shared/` behind
    a `HAS_LONG_READS` gate instead. Code was shared where it is genuinely the same, and
    left alone where forcing it together would have made it harder to read.

## Validation of the consolidation { #the-merge-was-checked-not-assumed }

The consolidation was evaluated by rerunning each read-based mode on its previous
test dataset and comparing the delivered assemblies byte for byte. The table records
those tests; it does not establish equivalence for every possible input.

| Mode | Compared against | Genome chain |
|---|---|---|
| `illumina` | BacFlux 1.3.1 | draft and decontaminated assembly identical |
| `nanopore` | BacFluxL | Flye assembly, dnaapler reorientation, decontamination, Medaka consensus and delivered genome identical |
| `hybrid` | BacFluxL+ | SPAdes draft, Flye assembly, Medaka consensus and the Polypolish-corrected genome identical |
| `contigs` | — | no earlier run to compare against; exercised on its own code with assemblies from both SPAdes and Flye |

One downstream difference is deliberate and appears in the two long-read modes: Bakta
gains or loses a single feature in each, because BacFlux now hands it a `--replicons`
table declaring a closed chromosome circular. That makes Bakta's gene caller run in
closed mode, which shifts a couple of marginal start codons anywhere in the genome. This
change reflects different topology information supplied to the gene caller; the
comparison alone does not establish which individual feature call is correct.

The helper scripts that carry real logic — the decontamination selector, the Medaka
model check, the replicon table, CARD reporting, plasmid concordance and the mobilome
helpers — are plain Python and are tested without Snakemake:

```bash
pytest workflow/scripts        # three tests skip unless BACFLUX_TEST_* variables name real GTDB files
```

## What did not change

The analysis steps, the thresholds and the decontamination policy came across as they
were, and the byte comparisons above are the evidence. The merge only added: QUAST now
runs in `nanopore` mode, which never had a QUAST rule, and in `hybrid` it evaluates the
delivered genome as well as the Illumina draft. Neither touches an existing output. What
moved is where the output lands, and how the run is configured.

## Where to go next

| Page | Read it for |
|---|---|
| [Coming from v1](../getting-started/from-v1.md) | the `mode` key, the config to rebuild, the renumbered output, the retired repositories |
| [Choosing a mode](../getting-started/choosing-a-mode.md) | which mode a dataset belongs in, and what a wrong one does |
| [Modes](../modes/index.md) | what each front end runs, step by step |
| [Output files](../reference/output.md) | what each numbered stage holds, and which files are the deliverables |
| [The mobilome module](../mobilome/index.md) | the part that is new rather than consolidated |
