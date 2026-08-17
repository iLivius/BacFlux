# Choosing a mode

BacFlux is one workflow with four front ends. The `mode` key picks the front end;
everything downstream of the assembly is the same code in all four. That is what
makes two isolates comparable when one was sequenced on Illumina and the other on
Oxford Nanopore Technologies (ONT from here on).

The choice follows from what is on disk, not from preference.

| `mode` | You have | File names BacFlux looks for | Input key |
|---|---|---|---|
| `illumina` | Illumina paired-end reads | `{sample}_R1.<ext>` and `{sample}_R2.<ext>` | `input.illumina_dir` |
| `nanopore` | ONT long reads | `{sample}_ont.<ext>` | `input.nanopore_dir` |
| `hybrid` | Both, for the same isolates | both of the above | both directories |
| `contigs` | A finished assembly and no reads | `{sample}.<fasta\|fa\|fna>` | `input.contigs_dir` |

`<ext>` for reads is one of `fastq`, `fq`, `fastq.gz`, `fq.gz`, and every sample in
one batch must share the same one — a directory holding both `.fastq` and
`.fastq.gz` stops the run and lists what it found. The same applies to the FASTA
extension in `contigs` mode. Sample names may contain underscores and dots;
spaces and the characters `* # @ % ^ / ! ? & : ; | < >` are rejected.

## The switch

Two keys in `config/config.yaml` pick the front end: `mode`, and the input
directory or directories that mode reads. The input directories it does not read
can stay empty — the output directory and the database paths under `directories:`
are needed whichever mode you pick.

```yaml
mode: illumina                 # illumina | nanopore | hybrid | contigs

input:
  illumina_dir: "/data/run1/illumina"
  nanopore_dir: ""
  contigs_dir:  ""
```

Both are checked before any job starts, so a mistake costs you seconds rather than
hours. Write `mode: ilumina` — one `l` short — and the run stops with:

```text
config.mode must be one of illumina | nanopore | hybrid | contigs (got: 'ilumina')
```

The message quotes back exactly what it found, so the typo is visible rather than
merely reported.

and a valid run echoes the mode and every sample it found, so the batch can be
confirmed from the top of the log:

```text
Mode: illumina — Genomic analysis of bacterial Illumina reads.
Sample isolate_01 will be processed.
Sample isolate_02 will be processed.
```

The mode can be overridden for a single run with `--config mode=contigs`, but the
matching `input.*_dir` almost always has to change with it, so editing the file is
usually simpler.

## What each front end does

| | `illumina` | `nanopore` | `hybrid` | `contigs` |
|---|---|---|---|---|
| Read QC | PhiX removal with Bowtie2, then fastp | NanoPlot before and after, Filtlong | both of those | none |
| Assembler | SPAdes `--isolate`, k chosen from read length | Flye | Flye, over the ONT reads | none |
| After assembly | keep contigs ≥ 500 bp and ≥ 2× coverage | dnaapler reorientation, then Medaka (*optional* — `parameters.<mode>.medaka_model: false` skips it) | dnaapler and Medaka as in `nanopore`, then Polypolish correction with the Illumina reads | SPAdes-style headers: the same 500 bp / 2× filter, headers left as they are. Any other header: trimmed to its first whitespace token, nothing filtered |
| Contamination screen | on the SPAdes draft, and its output *is* the delivered genome | on the reoriented assembly, before Medaka | on the **Illumina** draft only; the ONT leg is screened at the read level instead | on the supplied contigs, and its output *is* the delivered genome |
| Delivered genome | the screened SPAdes contigs | the Flye assembly, reoriented, screened and polished | the ONT assembly corrected by Polypolish | the screened input contigs |
| Also produced | — | — | the Illumina draft, kept as a comparator, and a Snippy comparison of the two | — |

`hybrid` assembles twice — SPAdes over the short reads and Flye over the long ones
— so it is the most expensive mode to run. The Illumina draft is not a by-product:
it is what the contamination screen screens, and the reads that map to the clean
version of it become Filtlong's reference for scoring the ONT reads. The ONT leg
has no contamination screen of its own.

## Everything after the assembly is shared

Each front end ends by writing one file, `contigs_final.fasta`, and from there the
run is identical in all four modes: assembly QC, completeness and contamination,
taxonomic placement, annotation, AMR and virulence screening, plasmids, prophages,
the optional mobilome module, and one MultiQC report. See
[Output files](../reference/output.md) for what each stage writes.

Decontamination is the one shared step that is not downstream of that file. It is
the same screening code in every mode, but it sits inside the front end, at the
different points the table above gives — which is why in `illumina` and `contigs`
its output *is* `contigs_final.fasta`. See
[Decontamination](../analysis/decontamination.md).

## What the mode still decides downstream

A few shared steps need something a given mode does not have. They are gated on
the capability — short reads, long reads, any reads — rather than on the mode name,
and when the capability is absent the rule is not defined at all, so it does not
appear in the DAG and nothing waits for it.

| Shared step | Needs | `illumina` | `nanopore` | `hybrid` | `contigs` |
|---|---|:-:|:-:|:-:|:-:|
| CARD read mapping (`05.amr/mapping`) | short reads | ✓ | — | ✓ | — |
| Qualimap mapping evaluation | any reads | ✓ | ✓ | ✓ | — |
| Coverage that carries information for BlobTools | any reads | ✓ | ✓ | ✓ | — |
| Circular-replicon table handed to Bakta | Flye's topology calls | — | ✓ | ✓ | — |
| Second BLAST of the final assembly, for the plasmid stage | long reads | — | ✓ | ✓ | — |
| IS copy number from reads (mobilome, *optional*) | short reads | ✓ | — | ✓ | — |
| Two genomes through QC and taxonomy, plus the Snippy comparison | both technologies | — | — | ✓ | — |

The same table appears on [Modes](../modes/index.md), which is where it is maintained.
Four consequences worth knowing before you choose:

- **`contigs` mode cannot separate contaminants by coverage.** BlobTools needs a
  BAM, and with no reads the contigs are mapped against themselves, which gives a
  near-uniform depth that means nothing. Only the BLAST taxonomy leg is doing real
  work in that mode. See [Decontamination](../analysis/decontamination.md).
- **`nanopore` and `contigs` lose the read-based AMR leg.** Mapping reads to CARD
  side-steps the assembly, so a resistance gene present in several copies, or
  sitting on a repeat that broke the contig, can be visible in the reads and absent
  from the assembly. Without reads there is no second chance; ABRicate on the
  contigs still runs in every mode. See
  [Antimicrobial resistance](../analysis/amr.md).
- **Only the long-read modes know that a contig is circular.** Flye reports it per
  contig, and Bakta is told through `--replicons`, which lets the gene caller call
  a gene that runs across the origin — often *dnaA* itself. In `illumina` and
  `contigs` mode Bakta treats every sequence as linear.
- **The mobilome module runs in all four modes, and is believable to different
  degrees in each.** On a fragmented assembly the located insertion-sequence count
  is a floor rather than a count, because the elements being hunted are a leading
  cause of the contig breaks. Every call therefore carries its distance to the
  contig end, and anything spanning contigs is capped at low confidence. In the
  short-read modes, and only when an ISOSDB source is configured, a read-based
  copy-number leg turns that warning into a number.
  See [Draft assemblies](../mobilome/draft-assemblies.md).

## What a wrong mode does

Usually it stops within seconds, because the four filename conventions do not
overlap. Pointing `illumina` at a directory of `{sample}_ont.fastq.gz`:

```text
No suitable illumina input files found.
```

and `hybrid` at two directories that do not describe the same isolates:

```text
Missing ONT files for the following samples:
  - isolate_03
```

An unset or non-existent input directory is caught in the same pass:

```text
[BacFlux] mode=nanopore requires 'input.nanopore_dir' to be set in the config.
```

!!! warning "The two cases that are not loud"

    **`illumina` on a dataset that also has long reads.** It runs to the end and
    quietly ignores `nanopore_dir`. There is nothing to detect: the Illumina half
    is a perfectly valid `illumina` run. If you meant to use both, the mode is
    `hybrid`.

    **`contigs` on an assembly whose headers are not SPAdes-style.** No length
    filter and no coverage filter are applied — only the headers are trimmed to
    their first whitespace token, which every later step keys on. That is
    deliberate, because the length and coverage of each contig are read out of the
    SPAdes header and inventing them for another assembler would silently delete
    real sequence. But it means the only contigs that can leave the assembly are
    the ones the contamination screen drops. The rule logs which branch it took.

## Picking one

- **Reads from one instrument** — use that instrument's mode. Nothing about
  `hybrid` improves a dataset that is missing half of what it needs.
- **Both technologies for the same isolates** — `hybrid`, and give it two
  directories in which every isolate appears in both. A sample present in only one
  of them stops the run rather than being silently dropped, so filter the
  directories to the complete isolates first.
- **Closed replicons, plasmid topology, or an assembly that survives its own
  repeats** — you need long reads, so `nanopore` or `hybrid`. This is also the only
  route to the circular-replicon information Bakta uses.
- **A genome from a collaborator, a database, or an earlier assembly** — `contigs`.
  It runs every analysis stage; it just cannot re-examine evidence that only exists
  in reads.
- **Something already analysed with a workflow older than v2.0.0** — these four
  modes replace the separate repositories that used to do these jobs. See
  [Coming from v1](from-v1.md).

## Next

- [Quick start](quick-start.md) — a first run from clone to finished report.
- [Reference databases](databases.md) — what must be on disk before the first launch.
- [Illumina](../modes/illumina.md), [Nanopore](../modes/nanopore.md),
  [Hybrid](../modes/hybrid.md) and
  [Pre-assembled contigs](../modes/contigs.md) — each front end step by step.
- [Configuration](../reference/configuration.md) — every key, which modes read it,
  and its default.
