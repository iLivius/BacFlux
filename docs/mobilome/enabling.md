# Turning it on

```yaml
mobilome:
  run: true
```

That is the whole switch. **No database path has to be configured**, in any mode.

Everything else in the `mobilome:` block is either a threshold with a sensible default
or one of the four optional layers, each of which stays off until you give it a source
of its own ([Optional layers](optional-layers.md)).

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **AMR** | antimicrobial resistance |
    | **IS** | insertion sequence — the smallest mobile element, carrying only the genes it needs to move itself |
    | **ICE** | integrative and conjugative element — sits in the chromosome and carries its own conjugation machinery |
    | **IME** | integrative mobilisable element — sits in the chromosome but needs a helper element to move |
    | ***att* site** | the short direct repeat left at each end of an element where it integrated |
    | **HMM** | hidden Markov model — a statistical profile of a gene or protein family |

## What `run: true` adds

Ten rules on the default configuration, of which three are per-sample tool runs at up to
8 threads.

| Added | What it does |
|---|---|
| **AMRFinderPlus** | the AMR calls the module reasons about, with coordinates, a `Method` column and a partial-at-contig-end flag |
| **ISEScan 1.7.3** | insertion sequences on the delivered genome |
| **MacSyFinder 2.1.6** running the **CONJscan** models | relaxase, coupling protein and mating-pair apparatus — the machinery behind an ICE, an IME or a conjugative plasmid |
| four helper steps | replicon calls from Platon, contig lengths, the ICE/IME caller with its *att*-site search, and the co-localisation that writes the deliverable |

Two of those are new conda environments, both from bioconda. AMRFinderPlus is not new:
it is already installed by Bakta, and reads its database from
`{bakta_db}/amrfinderplus-db/latest` — which is why the module needs no download of its
own.

Stage 08 **consumes** earlier stages rather than repeating them. It reads the delivered
genome, the Bakta proteins and GFF3, the GTDB-Tk species call, and Platon's replicon
calls from `06.plasmids`. Nothing in it re-detects a plasmid or a prophage.

!!! note "The first run needs network access"

    The CONJscan model package is fetched at run time by `msf_data install`. It is a
    *versioned* package, so unlike the optional layers it needs no `PROVENANCE.txt` —
    the version is already recoverable. The rule writes the models' licence notice to
    the top of its own log.

## What runs in which mode

The ICE/IME caller, the composite-transposon test and the *att*-site search run in
**all four modes**, `illumina` included. The constraint on them is assembly contiguity,
not the sequencer — and every call carries the flags that say how contiguous the
assembly under it was ([Draft assemblies](draft-assemblies.md)).

The one mode-restricted leg is the optional **ISOSDB copy-number** layer, which maps
reads and therefore needs them: `illumina` and `hybrid` only.

## `--organism`, and why the taxonomy stage matters here

AMRFinderPlus detects curated **point mutations** only when it is told which organism it
is looking at. Those mutations are the intrinsic, chromosomal, non-transferable
determinants that tier 1 rests on, and a homology screen cannot see them
([The mobility ladder](mobility-ladder.md)).

So the module maps the GTDB-Tk call onto an AMRFinderPlus organism name and passes
`--organism` when a confident mapping exists. There is no official convention for that
translation — GTDB's own FAQ states there is no direct mapping to NCBI taxa — so BacFlux
uses a hand-checked, project-local table. It is documented and audited, not
authoritative:

- every decision, **including every refusal**, is written to
  `{sample}_amrfinder_organism_audit.tsv` with its reason;
- a genome with no confident mapping gets **no** `--organism` rather than a guess;
- an empty mutations file therefore means *not assessed*, never *none found*.

Re-check the table when you change the GTDB release.

!!! note "The other AMR legs are untouched"

    ABRicate and the CARD read screen keep running exactly as before. The three legs are
    complementary: reads are immune to assembly collapse, ABRicate brings breadth across
    eight databases, and AMRFinderPlus brings the structured output this module needs.
    See [Antimicrobial resistance](../analysis/amr.md).

## What it costs

**Time.** Three extra tool runs per sample, each capped at 8 threads, plus a one-off
model download. On a typical bacterial isolate none of them is close to eggNOG or GTDB-Tk
in cost.

**Disk.** The model package and the per-sample `isescan/` and `conjscan/` trees, kept as
raw evidence.

**A licence decision.** Turning the module on fetches the CONJscan models, which Institut
Pasteur / CNRS license under CC BY-NC-SA 4.0. BacFlux never ships them — you download
them, under your own agreement with the licensor, exactly as with
`bakta_db`. Selecting the module prints the notice at parse time. Licence types for
everything the module can fetch are on [Licensing](../about/licensing.md).

## The rest of the block

```yaml
mobilome:
  run: false
  max_composite_span_bp: 20000       # two IS copies further apart are not one composite
  coverage_profile: 0.5              # how much of an HMM profile a hit must cover
  contig_boundary_bp: 100            # how close to a contig end counts as "at the edge"
  require_trna_boundary_for_high: false
  icescan:   { run: false, ... }     # ─┐ the four optional layers, all off by
  tncentral: { url: "", ... }        #  │ default: ICEscan by its own run flag,
  iceberg:   { urls: [], ... }       #  │ the other three because they have
  isosdb:    { fasta_url: "", ... }  # ─┘ nowhere yet to fetch from
```

ICEscan is the only one of the four that ships a working URL, because it has a single
canonical source — the bundle EBI publishes for ICEfinder2 — so `run: true` is all it
needs. For the other three you paste a URL in yourself, and the commented `# e.g.` line
beside each key is the address that was used to validate this module.

Two of them need explaining before a first run.

- **`coverage_profile`** (default 0.5, MacSyFinder's own default) is the fraction of a
  profile HMM an alignment must cover before the hit is kept. It is not identity and not
  an E-value: a hit can be statistically beyond doubt and still be dropped by this rule
  alone. Every validated result in these docs was measured at 0.5. What lowering it buys
  and costs, over all three benchmark sets, is on [Tuning](tuning.md).
- **`require_trna_boundary_for_high`** (default false) demands a tRNA-anchored *att* site
  before any element may be called high confidence. Sensible on closed long-read
  assemblies, where an unresolved boundary is a warning sign; on a draft it is the norm,
  which is why it is off by default.

Every key, in the file's own order, is in
[Configuration](../reference/configuration.md).

## Confirming it is on

The parse-time header names the module before any job starts, alongside the mode and the
phage caller. If it is not in the header, it is not running — check that you edited the
config you passed on the command line.

Output lands in `08.mobilome/`: one directory per sample, plus the shared model and
database directories at the top of the stage. See
[Reading the output](output.md).
