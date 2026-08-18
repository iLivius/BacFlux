# Rationale

Getting a bacterial isolate from raw sequencing data to a statement you can defend —
what species it is, how complete the genome is, what the organism can do, what
resistance it carries and whether that resistance can move — takes a dozen tools.
Each has its own input format, its own reference database and its own defaults. Run
by hand, the chain is slow, easy to get subtly wrong, and hard to repeat a year later.

BacFlux is that chain written down once, as a [Snakemake](https://snakemake.github.io/)
workflow. Every step declares its inputs and its outputs, and every step that calls an
external tool declares its own conda environment, so a run can be resumed after an
interruption, repeated on another machine, or handed to a colleague along with one
config file and nothing else.

That much is true of any workflow manager. What is specific to this repository is that
there is now a single workflow where there used to be four.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **AMR** | antimicrobial resistance |

## Four ways in, one analysis

BacFlux began as one workflow, for Illumina paired-end reads. Three more followed, each
because a project needed it — FastaFlux for genomes that arrived already assembled,
BacFluxL when Oxford Nanopore Technologies (ONT) sequencing became routine here, and
BacFluxL+ for isolates sequenced both ways. FastaFlux was a second Snakefile alongside
the first; the two long-read ones had repositories of their own. All four were released,
and work published with them stays citable — BacFluxL and BacFluxL+ are retired in
favour of this repository and their DOIs stay valid. See
[Coming from v1](../getting-started/from-v1.md).

They also had almost everything in common. Once an assembly exists, nothing downstream
of it cares how it was produced. Taxonomic placement, annotation, AMR and virulence
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
instead of four times. [Choosing a mode](../getting-started/choosing-a-mode.md) covers
which one a given dataset belongs in.

## What writing a step once buys

**Results become comparable across sequencing strategies.** The comparison that comes up
most often here is between technologies: the same isolate assembled from short reads and
from long reads, or one collaborator's Illumina genome against another's ONT genome. When
those two answers come out of two sets of rules, every difference between them is
confounded by the code, and you cannot separate a real biological difference from a
threshold that was changed in one copy and not the other. Out of one set, the difference
is the data. Hybrid mode is the sharpest case: QUAST, CheckM and GTDB-Tk read both
assemblies of the isolate, and Snippy counts the variants between the Illumina assembly
and each stage of the ONT one. Those numbers mean something because the same code
produced both sides.

**A fix reaches every mode at once.** A re-run in July found that the CheckV database
ships without a DIAMOND index, so the index in any shared copy was built by whichever
DIAMOND version that site happened to have, and the pinned one could not read it. The
fix — build a local view of the read-only database from symlinks plus an index this
workflow builds itself — was written in one rule, and all four modes had it. Two of the
three bugs found that week were invisible in the `illumina` and `nanopore` runs, because
those samples call zero viral contigs and take the empty-input shortcut before ever
touching the database. Only the `hybrid` sample, which carries a prophage, went down
that path — and fixing it there fixed it everywhere.

**A new capability is built once.** The mobilome module is 22 rules and eleven Python
helpers. Across the four it would have been copied four times and maintained in
parallel. That is what settled the order of work: merge first, then write the module.
See [the mobilome module](../mobilome/index.md).

## How the four modes actually share

**One hand-off file.** Whatever a front end ran, it finishes by writing
`02.assembly/{sample}/contigs_final.fasta`, and every rule from `03.taxonomy` onwards
reads that file and nothing else about the assembly. Exactly one shared rule feeding
those stages looks further back — the one that reads Flye's own contig summary to tell
Bakta which contigs are circular chromosomes — and it is gated so that it defines
nothing in the two modes without long reads.

**Fixed stage numbers.** All technology-specific work lives under `01.reads` and
`02.assembly`, so every shared stage lands on the same number in every mode: taxonomy is
always `03.taxonomy`, AMR always `05.amr`. Previously each workflow had a different
number of front-end stages, so the shared stages were numbered differently in each — the
single largest reason "the same" rule was not literally the same file. See
[Output files](../reference/output.md).

**Capability flags rather than mode names.** Where a shared step genuinely depends on the
technology, it asks what the run *has*, not which mode it *is*: `HAS_SHORT_READS`,
`HAS_LONG_READS`, `HAS_READS`, defined in `workflow/rules/shared/00_common.smk`. The CARD
read-mapping leg is gated on `HAS_SHORT_READS` because it needs reads — not because of a
hard-coded list of mode names that a fifth mode would have to be added to.

**One environment set.** 31 files in `workflow/envs/`, one per tool, with the version
pinned in the file — plus, in a few of them, a second pinned package the tool needs at
run time (`diamond` alongside Bakta, `samtools` alongside bowtie2). The merge turned up
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

## The merge was checked, not assumed

"Consolidation changed the plumbing, not the biology" is a testable claim, so it was
tested. Each read-based mode was re-run on the isolate its predecessor had been run on —
an *Arthrobacter* — and the delivered genome compared byte for byte against the earlier
result.

| Mode | Compared against | Genome chain |
|---|---|---|
| `illumina` | BacFlux 1.3.1 | draft and decontaminated assembly identical |
| `nanopore` | BacFluxL | Flye assembly, dnaapler reorientation, decontamination, Medaka consensus and delivered genome identical |
| `hybrid` | BacFluxL+ | SPAdes draft, Flye assembly, Medaka consensus and the Polypolish-corrected genome identical |
| `contigs` | — | no earlier run to compare against; exercised on its own code with six genomes, three of SPAdes and three of Flye origin |

One downstream difference is deliberate and appears in the two long-read modes: Bakta
gains or loses a single feature (4,999 → 4,998 in `nanopore`, 4,773 → 4,774 in `hybrid`)
because BacFlux now hands it a `--replicons` table declaring a closed chromosome
circular. That makes Bakta's gene caller run in closed mode, which shifts a couple of
marginal start codons anywhere in the genome. It is 0.02% of the features, and it is the
better call.

The helper scripts that carry real logic — the Medaka model check, the replicon table,
CARD reporting, plasmid concordance and the mobilome helpers — are plain Python and are
tested without Snakemake:

```bash
pytest workflow/scripts        # 457 passed, 2 skipped
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
