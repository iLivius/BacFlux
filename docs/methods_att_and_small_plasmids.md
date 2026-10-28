# Method reference: att-site detection and small-plasmid recovery

This page explains two limitations relevant to interpreting mobile-element context:
uncertain element boundaries and incomplete small-plasmid recovery. It describes the
methods used in BacFlux, the supporting observations and the limits of those observations.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | ***att* site** | the short repeated sequence left at each end of an element when it integrates into a genome. A supported pair can help locate an element's boundaries |
    | **attL / attR** | the left and right copies of that repeat |
    | **ICE** | integrative and conjugative element — a mobile element that sits in the chromosome and can move itself to another cell |
    | **direct repeat** | two copies of a sequence in the same orientation |
    | **HMM** | hidden Markov model, a statistical profile used to recognise a protein family |
    | **IS** | insertion sequence, the smallest kind of mobile element |
    | **MEM** | maximal exact match — a shared stretch that cannot be extended further at either end |
    | **ONT** | Oxford Nanopore Technologies, the long-read sequencing platform |
    | **DAG** | the graph of rules Snakemake builds to decide what still has to run |

The evidence comes from the primary sources, the tools' own source code, and two
failures on clinical *Klebsiella pneumoniae* hybrid data — isolates TUM24772
(PRJNA1168299) and K3 (PRJNA1291976).

---

## Part 1 — How att sites are actually detected

### The failure that prompted this

BacFlux's att search took a **25 bp probe** from the 3′ end of a tRNA and looked for
a second copy bracketing the candidate element, allowing one mismatch. On two clinical
ICE*Kp* elements it found no pair.

The reason is measured, not inferred: **the ICE*Kp* direct repeat has a 17 bp core**
(`CCAGTCAGAGGAGCCAA`; Lam *et al.* 2018), of which only 16 bp lie inside the tRNA-Asn
gene. A probe cut from the last 25 bp of the tRNA therefore extends 9 bp beyond the
sequence the two copies share. With a variable-length search, the pair recovered on
these isolates is a 23 bp exact repeat: the 17 bp core followed by `ATTTAA`.

> Lam MMC *et al.* (2018) *Genetic diversity, mobilisation and spread of the
> yersiniabactin-encoding mobile element ICEKp in Klebsiella pneumoniae
> populations.* Microbial Genomics.
> <https://pmc.ncbi.nlm.nih.gov/articles/PMC6202445/>

### What the reference tools do — every one searches a range of lengths

| tool | att method | length rule | source |
|---|---|---|---|
| **ICEfinder** (Ou lab, the reference ICE tool) | ARAGORN locates tRNA/tmRNA 3′ termini; **Vmatch** finds the direct repeats marking the tRNA-distal boundary | maximal exact repeats — variable by construction | [ICEberg 2.0 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6323972/) |
| **ICEfinder2** | same, in code | **`vmatch -l 15`** — floor 15 bp, length read from the output | [`script/single.py` L286](https://github.com/EBI-Metagenomics/icefinder2/blob/main/script/single.py) |
| **icefinder-opt** (maintained 2025 fork) | unchanged | still `vmatch -l 15` | [repo](https://github.com/guogenglin/icefinder-opt) |
| **DEPhT** | **BLASTN of the left flank against the right flank**, scoring *all* resulting pairs | whatever the local alignment returns | [DEPhT paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9303363/) |
| **DBSCAN-SWA** | same flank-vs-flank BLAST pattern | explicit **12 bp** minimum, ranked by bitscore | [`bin/dbscan-swa.py`](https://github.com/HIT-ImmunologyLab/DBSCAN-SWA/blob/master/bin/dbscan-swa.py) |
| **Islander / TIGER** | sidesteps repeats entirely — integration splits a tDNA and the island restores it, so it hunts the **displaced tDNA fragment** by sensitive BLASTN | whatever the alignment returns | [Islander paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4383910/) |
| **IntegronFinder** | covariance model (Infernal) for *attC* | bounded variable length, defaults 40–200 bp | [NAR 2016](https://academic.oup.com/nar/article/44/10/4539/2516972) |
| IslandViewer components (IslandPath-DIMOB, SIGI-HMM), alien_hunter, PAIDB | **do not detect att sites at all** — compositional / mobility-gene methods | n/a | [IslandViewer 4](https://academic.oup.com/nar/article/45/W1/W30/3787837) |

IslandViewer integrates Islander *specifically because* its own compositional
predictors have poor boundary accuracy — so they offer no precedent for
boundary calling.

### Why repeat length is allowed to vary { #why-a-fixed-length-cannot-work }

- Bacterial ICE direct repeats span roughly **10–60 bp**
  ([FEMS Microbiol Rev](https://academic.oup.com/femsre/article/41/4/512/3089980)),
  so any single value sits arbitrarily inside the real distribution.
- Integrase family predicts the range: **tyrosine**-integrase elements at tRNA
  sites (which ICE*Kp* is) sit at **20–45 bp** cores, **serine**-integrase
  elements go as low as **3–12 bp** ([DEPhT](https://pmc.ncbi.nlm.nih.gov/articles/PMC9303363/)).

### Maximal repeats and local alignment { #stepping-the-probe-down-no-precedent-and-weaker-anyway }

A maximal-repeat search avoids testing a sequence of fixed probe lengths: it
returns matches that cannot be extended further. Local alignment is a different
approach and can allow mismatches or gaps. An exact-repeat search, including
BacFlux's, does not capture those imperfect matches. A BLAST bit score measures
alignment support rather than the biological validity of an integration boundary.

### Boundary-control measurements { #nobody-has-published-a-false-positive-rate-for-boundary-calling }

No published false-positive rate for *att*-site boundary calling was found when this
method was developed (July 2026): benchmarks of genomic-island tools score island calls,
not boundaries, and the Islander paper reports only "the few false positives", without a
number. BacFlux therefore measured its own rate, on 300 randomly placed 15 kb non-ICE
spans of a clinical *K. pneumoniae* chromosome, with the real IS mask applied: **22%**
of them returned a confident de novo "boundary", falling to **16%** once the
repeat-family guard was added, and to **~1%** reaching the AMR table once de novo
repeats stopped being applied at all. Quote **16%** as the de novo false-boundary rate
in this specific control experiment — 22% is the *before* figure, and the two are easy
to confuse because both appear in the code comments.

This is also why de novo repeats are **reported but not applied**: at a
one-in-six error rate a de novo pair is a reasonable lead for a human to follow
and an unacceptable basis for silently redefining an element, since widening an
interval turns every gene inside it into predicted cargo. tRNA-anchored pairs
*are* applied, because they start from a position integrases are known to target.

### What BacFlux does as a result

**Exact maximal repeats between the two flanks, computed in standard-library
Python** — `vmatch -l N` semantics, which is what ICEfinder and ICEfinder2 use.
The search is in `workflow/scripts/80_mobilome/att_search.py`
(`find_maximal_repeats` / `search_maximal_repeat`).

It calls no external program. `att_search.py` is standard-library only, and the rule
that runs it, `conjscan_ice`, uses `workflow/envs/python.yaml`, a pinned Python with
nothing else in it. Shelling out to `blastn` would mean adding BLAST to that
environment. A BLAST search would be more sensitive — it absorbs mismatches and gaps,
which an exact repeat search cannot — so this is a deliberate trade of sensitivity
for a rule that needs nothing beyond Python, not an equivalence.

How it works, and where each choice comes from:

- **Exact maximal repeats**: index every *k*-mer of the left flank, find matches
  in the right flank, then extend each match outwards for as long as the two
  sequences agree. One pass returns the longest repeat, so nothing has to "step
  the probe down".
- **A 50 kb flank window each side** (`DEFAULT_FLANK_WINDOW_BP`), raised from
  30 kb. Measured on SPI-7 (`AL513382`): the correct attR sits in a tRNA-Phe
  5,800 bp *outside* a 30 kb window, so the element was reported 50 kb short for
  want of anywhere to look. Widening to 50 kb recovers it, and 80 kb, 120 kb and
  200 kb change no element's answer on the benchmark — so this is where the curve
  flattens, not simply a bigger round number. The cost is more candidate repeats
  to rank, which is why the ranking rule below matters more at 50 kb than at 30.
- **Floor of 15 bp**, matching ICEfinder2 and icefinder-opt; DBSCAN-SWA uses 12,
  so 15 is the more conservative of the two published choices. The floor actually
  applied is the **larger** of 15 and a chance-match threshold computed from the
  two window sizes: at 15 bp, two 50 kb flanks are expected to share **~2.3**
  repeats by luck alone, so a flat 15 would report noise on any large window. At
  the shipped 50 kb window this resolves to **18 bp** — as it also did at 30 kb,
  so widening the window did not move the floor.
  Note this threshold assumes DNA is a random string of four equally likely
  letters, which real chromosomes are not; it is a floor against random
  background, not a filter against repetitive sequence. The 16% figure above is
  what happens when it is asked to do more than that.
- **One search, not two modes.** tRNA proximity is a **ranking term inside the
  single search**, not a separate search: every maximal repeat that brackets the
  machinery is found first, and only then are they ranked. This follows
  ICEfinder's division of labour — the tRNA *locates* a candidate site, the
  repeat search *delimits* it — and DEPhT's treatment of integrase proximity as
  one term among several rather than as a gate.
- **Anchoring outranks length.** Candidates are ordered by whether *at least one*
  copy sits inside an annotated tRNA **first**, and by repeat length only within
  each group. The two are not comparable quantities: length says how unlikely the
  match is by chance, a tRNA 3′ end says the match is where integration actually
  happens. Measured on SPI-7 (`AL513382`), a 51 bp repeat in ordinary sequence
  beat the real 24 bp *att* pair at tRNA-Phe and the element came out 50 kb short;
  under this ordering the tRNA-anchored pair wins and the call lands on the
  curated interval.
- **Both copies inside tRNAs of the *same amino acid* is rejected outright** —
  that is two paralogous tRNA genes, not an integration scar. Two copies in
  tRNAs of *different* amino acids is allowed, and the distinction is not
  academic: ICE*Ec2* (`GU725392`) has its real 22 bp *att* pair in tRNA-Phe at
  one end and tRNA-Ser at the other, and the earlier blanket rule threw that
  boundary away, reporting the element 37 kb short. This is why the ordering
  above tests "at least one copy in a tRNA" rather than "exactly one": by the
  time ranking happens, the paralogue case has already been removed.
- The label written to `boundary_method` is therefore an **outcome**, not a mode:
  `tRNA` when the winning repeat has a copy in a tRNA, `denovo` when it does not,
  `none` when nothing survived.

The implementation uses Python's standard library and requires no external repeat finder.

### Why not Vmatch

Vmatch is packaged in Bioconda (2.3.1, linux-64 and osx-64). Its recipe declares
`license: Unknown / OTHER`, and the upstream site could not be reached when this was
checked in July 2026, so its terms could not be established. BacFlux therefore does
not install it (see [Licensing](about/licensing.md)).

A separate check implemented the same intended exact-repeat search using a
prefix-doubling suffix array, Kasai LCP array and cross-flank MEM enumeration. It
replayed actual caller windows after IS masking, at the recorded coordinates and
repeat-length floor. The check script is not distributed.

The two implementations returned identical sets on **52 benchmark genomes, 246 att
searches, 180 distinct flank windows and 1,033 repeats**. Neither reached its internal
cap. On those windows, a search with Vmatch's maximal-repeat semantics finds no repeat
that BacFlux's search misses; the Vmatch executable itself was not run. Agreement
between two implementations of one definition does not validate a biological
boundary.

### Licensing note

ICEfinder2 is distributed under CC BY-NC-SA 4.0. Its source was read to establish which
algorithm it uses and with which parameters — maximal exact repeats between the two
flanks, with a 15 bp minimum — and no code was copied or adapted. Maximal exact repeat
search is standard string processing (REPuter, Vmatch); `find_maximal_repeats` is an
independent standard-library Python implementation written from that definition.

---

## Part 2 — Why small plasmids go missing, and what to do

### Assembly method and small-plasmid recovery { #it-is-not-the-assembler }

Assembler choice can affect small-plasmid recovery. The cited comparison reported
the following recovery of plasmids **<10 kb**
([Microbial Genomics 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10272865/)):

| assembler | <10 kb recovery |
|---|---|
| Unicycler (hybrid) | 100% |
| Canu | 100% |
| **Flye `--nano-raw`** | **79%** |
| Flye `--meta` | 73% |
| **Flye `--nano-hq`** | **67%** |
| Miniasm/Minipolish | 64% |
| Raven | 39% |

Wick & Holt independently report NECAT "failed to assemble many plasmids",
NextDenovo "performed poorly on plasmid assembly", and Raven "wasn't good with
small plasmids"
([[Genome Biol]](https://pmc.ncbi.nlm.nih.gov/articles/PMC6966772/)).
These results describe the tested versions and datasets; they do not establish
that switching assemblers would improve or worsen a particular BacFlux run.

> **Actionable side-finding:** `--nano-hq` is the *worst* of Flye's three modes
> for small plasmids (67% vs 79% for `--nano-raw`), and BacFlux auto-selects the
> mode. Worth revisiting if small replicons matter more than consensus accuracy.
> Note also that Flye's `--plasmids` flag **no longer exists** — added in 2.4,
> removed in 2.9.

### The two causes that do explain it

**1. ONT ligation library prep.** Plasmids <20 kb are under-represented in
ligation read sets by a mean factor of **~4**, and up to **>100-fold** for the
smallest (one 2.4 kb plasmid). Unfragmented circular plasmids never receive a
blunt-end adapter. **Rapid (transposase) prep shows no such bias.**
A workflow cannot reconstruct a replicon reliably when the sequencing data do not
contain enough supporting sequence.
([Microbial Genomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC8549360/))

**2. Our own read QC.** BacFlux hard-coded `--keep_percent 90 --length_weight 10`;
filtlong's own default for the second is 1. Filtlong scores reads as

```
(Length^lw × MeanQ^mqw)^(1/(lw+mqw)) × WindowQ
```

and `--keep_percent` deletes the bottom of that ranking. Reads from small plasmids are
often short and can be disadvantaged when length dominates the score. Concatemers can
produce reads longer than one plasmid copy, so plasmid length is not an absolute upper
bound on read length.

Wick's Feb 2026 read-QC benchmark (5 genomes, 11 assemblers) found the same
failure and is worth reading carefully, because it does **not** say what we first
took it to say. His `Filtlong-defaults` arm ran filtlong at its **own defaults**
(`--target_bases` = 100 × genome size, so `length_weight` 1 — the value we now
ship), and that arm *"removed all reads below 10 kbp, essentially erasing these
plasmids"* in the 3 of 5 genomes that had sub-10 kb plasmids, with more
structural errors as a result. So it is not independent confirmation of our
`length_weight 10` finding: it is a warning that filtlong's default length
weighting is already enough to erase small plasmids once you cull hard to a
target. His recommendation is to go **below** the default —
`--length_weight 0 --window_q_weight 0`. Note this bears directly on `nanopore`
mode, which culls with `--target_bases` exactly as his arm did.
<https://rrwick.github.io/2026/02/05/read_qc_testing.html>

### Measured on a public clinical isolate

*Klebsiella pneumoniae* **TUM24772** — BioProject PRJNA1168299, closed genome
GCA_043950115.1, chromosome CP171785.1. A clinical test isolate, and **not** the
BAA-2146 mobilome positive control. Its 5,596 bp Col2 plasmid
(pMTY24772_Col2, CP171790.1) was missing from our assembly while sitting complete
in the Illumina data.

Reads mapping to that plasmid, for all four combinations of the two keys. Counts
are `minimap2 -x map-ont` against `CP171790.1`, keeping reads with an alignment
block of at least 1 kb; on that measure the **raw** read set holds **603**:

| | `length_weight 10` | `length_weight 1` |
|---|---|---|
| **`keep_percent 90`** — the old value | **93** | 413 |
| **`keep_percent 95`** | 500 | **502** |

**Both keys matter, and `keep_percent` is the one to reach for.** At the 95 that
now ships, changing `length_weight` is worth two reads; at the old 90 it is worth
320. The plasmid needed 90 and 10 acting together to disappear, and either change
on its own brings most of it back.

### So why is `length_weight` still 1?

Because the two keys protect different size classes, which only shows in the
read-length distribution. Reads per band, same sample — the raw ONT set holds
2,954 reads of 1–3 kb and 3,418 of 3–6 kb:

| | 1–3 kb | 3–6 kb |
|---|---|---|
| `keep_percent 90`, `lw 10` | **0** | 318 |
| `keep_percent 95`, `lw 10` | **0** | 2,666 |
| `keep_percent 90`, `lw 1` | 123 | 2,331 |
| `keep_percent 95`, `lw 1` | 605 | 2,935 |

`length_weight 10` empties the 1–3 kb band at **both** `keep_percent` values.
Raising `keep_percent` rescues the 3–6 kb band, which is why it rescued this
5.6 kb plasmid, but only lowering `length_weight` rescues 1–3 kb. So
the two settings affected different read-length bands in this sample. These
measurements do not define universal plasmid-size ranges protected by either setting.

Total bases barely moved across the arms (208 → 227 Mb), so none of this was ever
about depth — only about *which* reads survived. At the old 90/10, 3.1% of kept
reads were ≤6 kb, against 47.0% of the raw set.

**Fix applied:** `parameters.long_read_qc.keep_percent` now defaults to **95**
(was 90) and `length_weight` to **1** (was 10), both configurable. Note that
`length_weight` reaches filtlong in **hybrid mode only** — the `nanopore` rule
never passes `--length_weight`, so in two of the four modes the key does nothing.

### Tandem plasmid copies { #the-doubled-circle-artifact }

Our 23,928 bp contig is the 11,970 bp Col plasmid **duplicated within one
contig** (2 × 11,970 = 23,940). This is documented: start/end overlap in small
plasmids produces ~200% contiguity, and the Plassembler paper names Flye
specifically as having "multiplicated many small plasmids" on real data.

**dnaapler does not fix this.** dnaapler reorients a contig; rotation alone does not
remove duplicated sequence or resolve a whole-plasmid concatemer. `autocycler trim` does
handle it, but cannot distinguish artifact from genuine duplication and expects to run
inside the Autocycler pipeline.

### Options considered during development { #the-options-and-why-only-one-of-them-ships }

| option | what it would require | verdict |
|---|---|---|
| **`keep_percent` 90 → 95** (and `length_weight` 10 → 1) | nothing — a config change | **This is what BacFlux ships.** Recovers 502 of the plasmid's 603 reads — but *not* the assembled plasmid |
| **SPAdes rescue** | no new tool or database | Published approach: small plasmids "usually appear as circular contigs" in a short-read graph ([Wick/Judd/Holt](https://pmc.ncbi.nlm.nih.gov/articles/PMC9980784/)). We already write `assembly_graph_with_scaffolds.gfa` |
| **Plassembler** | a new tool, and the PLSDB database | **MIT, bioconda 1.8.3**, purpose-built. Pools reads that do *not* map to the Flye contigs and hybrid-assembles them with Unicycler — exactly where a Flye-absent plasmid lands. Can **reuse our existing Flye assembly** (`--flye_directory`), and writes empty outputs when it finds nothing, so the DAG never breaks. Cost: PLSDB is mandatory (`-d`, no skip flag) and states **no licence at all** — handle exactly as TnCentral/ICEberg: opt-in, and never redistributed (see [Licensing](about/licensing.md#components-you-switch-on-yourself)) |
| **Hybracter** | replacing the whole assembly stage | MIT, bioconda, wraps Plassembler — but replaces the *entire* assembly stage and nests Snakemake inside Snakemake, with resource coordination needed between the nested workflows |
| **Autocycler** | replacing the whole assembly stage | Not evaluated as a rescue method in this experiment. Aimed at chromosome consensus accuracy; its own paper reports the smallest plasmid tested (2.5 kb) was occasionally missed and needed manual curation. Its clustering discards any cluster seen in too few input assemblies — precisely the small-plasmid case. Its own docs say: *if small plasmids matter, add Plassembler* |
| **Trycycler** | manual curation per genome | Requires manual curation per genome; not part of BacFlux |
| **Union of multiple long-read assemblers** | running several assemblers per sample | Not tested in this experiment; all assemblers remain limited by the available read set |

### Scope of the assembler comparison { #why-newer-assemblers-are-not-in-the-table }

Switching assembler is the obvious alternative to the options above, and it cannot
be judged on evidence. The newer long-read assemblers — Myloasm, metaMDBG, LJA,
hifiasm — have been benchmarked for accuracy through 2025 and 2026, but no head-to-head
measurement of **small-plasmid recovery** for these assemblers was found (literature
checked in July 2026). That is the
number this section turns on, so there is nothing to compare against and none of
them can be recommended or ruled out here.

---

## Part 3 — Decontamination and optional short-read guidance { #part-3-the-compounding-trap-in-our-own-pipeline }

When `parameters.hybrid.short_read_guidance: true`, decontamination can also affect
which long reads reach Flye. Illumina reads selected by mapping to the retained draft
contigs become the Filtlong reference. Long reads from a discarded contig can then
lose support and be trimmed, split or removed.

This route is conditional: guidance is disabled by default. With guidance off, the
selected Illumina reads still affect Polypolish, but they do not guide Filtlong.
The final hybrid decontamination screen is advisory and does not remove sequence.

A plasmid contig is removed when the genus with the highest summed BLAST bit score for
it differs from the host's, which is not unusual because plasmids cross genus
boundaries. In a `contigs`-mode run of *K. pneumoniae* **ATCC BAA-2146** with
`mode: auto`, the 2,014 bp plasmid pMYS (`NZ_CP006660.1`) was assigned to
*Escherichia* (summed bit score 58,709) rather than *Klebsiella* (16,253), and was
removed, because `auto` keeps only the sample's majority genus. Each of up to 50
database sequences hit adds its bit scores to its genus, so the sums reflect how well
each genus is represented in `nt` as well as sequence similarity. In `contigs` mode
the contig is removed directly; the guidance route above applies only to `hybrid`.

**It did not cause the TUM24772 loss** — that plasmid survived decontamination
(verified: present in `contigs_filt.fasta`), and was lost at the read-filtering
step instead, to `keep_percent 90` and `length_weight 10` acting together. Both
mechanisms must be ruled out separately.

**How to investigate:** inspect `contig_taxonomy_decisions.tsv` and the underlying
taxonomic hits for any discarded plasmid-sized contig. If it belongs to the isolate,
re-run with `mode: include` naming both genera, or `mode: off`, for that sample
(`sample_overrides` allows a per-sample mode). `discard_no_hit: false` rescues only
contigs with no BLAST hit, under `exclude` or `off`, so it would not have kept pMYS.

Inspect a dry run to confirm which outputs will be rebuilt.
