# Validation

The mobilome module's ICE/IME caller has been measured three ways: against a set of
curated elements with published coordinates, head-to-head against the EBI Mobilome
Annotation Pipeline on the same genomes and with the same scoring code, and on
deliberately fragmented copies of those genomes to find out what a draft assembly
costs. This page is that work — what was tested, what it showed, and [what it does
not show](#what-the-benchmark-does-not-show).

Three terms used throughout. An **ICE** (integrative and conjugative element)
integrates into the chromosome and carries its own conjugation machinery, so it can
move itself into another cell. An **IME** (integrative mobilisable element)
integrates and carries a relaxase — the enzyme that nicks the DNA to start transfer —
but no apparatus of its own, so it needs a helper element to move. An ***att* site**
is the short direct repeat left at each end of an element when it integrates; finding
the pair is how an element's true edges are established.

!!! warning "Everything the module reports is a prediction from sequence"

    An element called **predicted self-transmissible** carries conjugation machinery
    that looks complete. It has not been shown to transfer. The confirmatory
    experiment is a filter or broth mating assay, and no number on this page
    substitutes for one. See [the mobility ladder](mobility-ladder.md).

## The benchmark set

Forty closed reference genomes downloaded from GenBank and used as-is, in three
groups. Two accessions carry both a curated ICE and a curated IME, so 40 genomes
yield 30 curated elements.

| Set | Genomes | What it is |
|---|--:|---|
| ICE pilot | 18 | Deposits carrying an ICE with curated coordinates in ICEberg, across 14 genera |
| IME pilot | 12 | The same for IMEs, 10 genera, sizes deliberately straddling the module's size floor (522 bp to 110 kb) |
| Negative control | 12 | Closed genomes with **no** curated ICEberg entry, 32,601,100 bp in total |

The negative set is three reduced endosymbionts, two streamlined marine genomes, a
deep-branching thermophile, five full-size reference strains (*E. coli* K-12 MG1655,
*B. subtilis* 168, *P. aeruginosa* PAO1, *S. aureus* N315, *S.* Typhimurium LT2) and
an archaeon, *Halobacterium* sp. NRC-1. The archaeon is a construction check: the
machinery models are bacterial, so any call there would be spurious by definition.

**How a call is scored.** ICEberg's curated interval is taken as truth. Among a
caller's calls on that genome the scorer keeps the one with the largest overlap, and
reports the start offset, the end offset, and the **recovered fraction** — overlapping
base pairs divided by the curated length. Recovered fraction is a recall measure and
nothing else; the scorer contains no precision or false-positive code at all, which
is why the negative control and the call-burden count below have to be read
alongside it.

!!! warning "The truth standard is a curation, not the truth"

    ICEberg's coordinates are neither perfect nor exhaustive, and this set contains
    one example of each failure: a curated entry whose deposited record cannot
    contain the element it stands for, and a negative-control genome that carries a
    real, well-known element ICEberg has no entry for. The first penalises every
    caller; the second rewards a correct call as if it were an error. Both are
    described where they bite, below.

## ICE pilot — 15 of 18

That is 9 of the 11 with real chromosomal context, and 6 of the 7 where the deposited
record *is* the element and there is no flanking sequence for an *att* site to sit in.

**The true ceiling is 17 of 18**, because one of the three misses is the broken
curated record. Each miss has a stated cause:

| Miss | Size | Why |
|---|--:|---|
| `SXT(HN1)` | 19 kb deposited | No relaxase and no conjugation apparatus in the record. Not findable; should be footnoted out of any future version of this benchmark |
| ICE*Psy*10 | 161 kb | No relaxase, no coupling protein and no VirB4 anywhere in the element. MacSyFinder's own rejected-candidates file records *"quorum of mandatory genes required (3) is not reached: 0"* |
| ICE*B2* (*Mycoplasma bovis*) | 37 kb | Machinery too divergent for the profile models to cover — a model-coverage limit that applies to any tool built on the same profiles |

**Boundaries.** On the nine detected chromosomal ICEs — the only ones where boundary
error is a meaningful question — the median absolute offset is **5,141 bp at the start
and 14,912 bp at the end**, and the median recovered fraction is **0.476**. On the six
standalone deposits the median recovered fraction is **0.669**.

Quote those two for BacFlux alone. The head-to-head below reports a different pair,
and the two are [not interchangeable](#the-094-figure-is-not-the-0669-figure).

## IME pilot — 5 of 12

Median recovered fraction 0.69. **This is the weaker half of the module and is
described that way everywhere it appears.**

The misses concentrate at the small end — 522 bp, 1,248 bp and 5,123 bp — where the
caller's 2 kb floor for IME-shaped machinery and its requirement for two anchor
classes both bite; and at SGI1 (42 kb) and MGI*Vvu*Tai1 (19 kb), which produced no
overlapping call at all.

There is a second, structural limit on this half of the module, and it is visible in
the code rather than in the score. The *att* search is always given the 8 kb ICE floor
and refuses any repeat pair implying an element smaller than that, while an IME is
admitted at 2 kb. Four of the five detected IMEs have machinery spans of 3,840,
4,959, 5,942 and 6,252 bp, all under the *att* floor, so **their boundaries can never
be resolved however clean the repeat pair is** — they keep `boundary_method = none`
and report the machinery span. The deliberate reading is that an *att* pair implying a
3 kb element is mostly noise; the honest reading is that IME extents are largely
unresolved. Both floors are constants in the script rather than config keys; see
[Tuning](tuning.md).

## Negative control — 2 calls in 32.6 Mb

Ten of the twelve genomes returned nothing at all, including *E. coli* K-12 MG1655,
*P. aeruginosa* PAO1 and *S. aureus* N315 — three large, well-annotated genomes dense
with recombinases, prophages and insertion sequences. The archaeon is among them,
which is the construction check passing.

Two calls were made, and neither survives as an error.

**_B. subtilis_ 168 → ICE*Bs1*, a real element ICEberg has not catalogued.** Called at
529,362–549,932, so 20,571 bp against an element of about 20.5 kb; all four anchor
classes present; a MOBT relaxase; and the boundary resolved by a tRNA-anchored *att*
site at a tRNA-Leu, which is where ICE*Bs1* is known to integrate. Bakta independently
annotated the integrase as "ICE*Bs1* integrase". This is the best-evidenced single
call in the whole validation.

**_S._ Typhimurium LT2 → an open question, not an answer.** An integrase sits beside a
protein Bakta labels only "DNA-binding protein", which the `T4SS_MOBM` relaxase profile
hits at i-evalue 1.8 × 10⁻⁸⁴ over 98.6% of the profile — a near-full-length relaxase
match. A strong relaxase beside an integrase is what an IME is. Whether this locus is
a *bona fide* IME or two co-located genes is a domain question, not a software one.

!!! warning "This is an upper bound of 2 candidate calls, not a false-positive rate"

    "Absent from ICEberg" is not "contains no element". The honest reading is **2
    candidates in 12 genomes, of which 0 are confirmed errors**.

    The set originally justified including *B. subtilis* 168 on the grounds that
    ICE*Bs1* is absent from 168 itself. **That is false — ICE*Bs1* was discovered in
    strain 168** (Auchtung *et al.* 2005). Had the calls been tallied instead of
    inspected, a textbook-correct detection would have been recorded as a false
    positive and a working detector would then have been "fixed" to stop making it.
    Every call in a negative control has to be examined, never merely counted.

## What the second model set adds

The optional ICEscan layer ([Turning it on](enabling.md)) supplies the integrase and
IME/AICE models that the default machinery model set has no equivalent for. It was
measured as a controlled two-arm experiment over the 28 pilot genomes: both arms
identical, the only difference being whether the ICE caller was given ICEscan's
table.

| | Default models only | With ICEscan |
|---|--:|--:|
| ICE pilot detected | 15/18 | **15/18 — unchanged** |
| IME pilot detected | 3/12 | **5/12** |
| Total calls | 51 | 63 |
| `ice` | 36 | **36 — unchanged** |
| `ime` | 10 | 21 |
| `aice` | 0 | 2 |
| Calls carrying a predicted self-transmissible claim | 36 | **36 — unchanged** |
| Confidence `high` | 16 | **16 — unchanged** |

**No element gained a self-transmissible claim.** Every net-new call in the ICEscan
arm is an IME or an AICE — tier 5 or no tier at all. That is the single most important
safety property of the layer: it can add elements at the mobilisable end of the ladder,
it cannot promote anything to the top rung.

Both new IME detections are honest recoveries rather than a bigger interval swallowing
the answer: Tn*4451* at 0.94 of its curated length (start offset +75 bp) and
`IME_CdiR20291_ND` at 0.69 (start offset +37 bp) — the two best-bounded IME recoveries
the caller has produced.

It costs something, and the cost is reported. One curated element, `CMGE(TZ080501)`,
loses 12,324 bp of recovered length because an ICEscan integrase sits 0 bp from the
machinery cluster and wins a "then the closest" tie-break against a more distant
integrase that happened to be nearer the true left edge. Several elements not covered
by either pilot shrink by 10–20 kb for the same reason. Whether "then the closest" is
the right tie-break when neither candidate yields an *att* pair is an open design
question, not a settled one.

## Head-to-head against the EBI Mobilome Annotation Pipeline

The same 40 genomes were run through the EBI Mobilome Annotation Pipeline v5.0.0, its
output converted into the same table layout, and **the same unmodified scoring script**
run over both. Any bias in that scorer applies equally to both callers. Two differences
are contained in the comparison rather than controlled for: the two pipelines call
genes with different tools (Prodigal there, Bakta here), which affects both the protein
set the shared models search and the integrase annotation; and the other pipeline was
run with its functional-annotation branches switched off, which is fair for a boundary
comparison but means nothing here evaluates its full output.

### These are not two independent callers

Both find conjugation machinery by running MacSyFinder with the ICEscan models — the
EBI pipeline inside its ICEfinder2-lite subworkflow, BacFlux in its `conjscan` /
`icescan` rules, with the optional layer switched on for this comparison. The
MacSyFinder versions differ slightly (2.1.4 there, 2.1.6 here) but
the profiles and the system definitions that decide *"is there a conjugation system
here at all"* are the same models.

| Layer | Shared or independent? |
|---|---|
| Machinery detection | **Shared** — same models, same engine |
| Clustering: which machinery genes belong to one element | **Independent** — different code, different rules |
| Boundary refinement: where the element starts and stops | **Independent** — vmatch plus an ICEfinder2-derived script there, an exact maximal-repeat search here |
| Class assignment: ICE / IME / AICE | **Partly independent** — same model vocabulary, different decision rules |

So the detection numbers are a consistency check — did two implementations of the same
idea behave the same way, and where did they diverge — and **must never be reported as
independent confirmation**. The boundary and class numbers are the real comparison.

### Detection

Over the 30 curated elements:

| | Count |
|---|--:|
| Found by both | 17 |
| Found by BacFlux only | **3** |
| Found by the EBI pipeline only | **0** |
| Found by neither | 10 |
| **Totals** | **BacFlux 20/30, EBI 17/30** |

The three are ICE*Vfl*Ind1, Tn*4451* and MTn*Pi*10; in each, the other pipeline's
prescan produced no candidate at all. Because the machinery models are shared, the
difference cannot be explained by a better model — it is attributable to what sits
around the shared engine: candidate seeding, a class-aware size floor, and the
decision to report low-quorum machinery at low confidence rather than discard it.

One honesty note on that count. ICE*Vfl*Ind1 is detected but classified
`cime_or_island` — integrase, coupling protein and apparatus, **no relaxase** — not
`ice`. It counts as detected because its call overlaps the curated interval, which it
does. It is not a correct class call and should not be presented as one.

### Class

Of the 17 elements both callers found, the classes agree on 14 and differ on 3. All
three are ICEberg-curated **ICEs that the other pipeline called IMEs**, and BacFlux
agrees with the curation on all three: SPI-7, `ICE_EfalV583_tRNALys` and
`ICE_Step12228_tRNAser`.

This is not cosmetic. ICE versus IME **is** the self-transmissible-versus-needs-a-helper
distinction — rungs 6 and 5 of the ladder, and the axis the intrinsic-versus-acquired
framing rests on. Calling a self-transmissible element "mobilisable with a helper"
understates transfer risk.

Two things keep this from being a claim that one caller is right about class in
general. There is a fourth element, TR2, where **both** callers say IME against an
ICEberg ICE curation — a shared miss, not counted among the three. And with three
disagreements out of seventeen shared elements, this is a consistent direction, not a
rate.

### Call burden

| | Calls | Total bp called |
|---|--:|--:|
| BacFlux | 65 | 2,630,059 |
| EBI pipeline | 64 | 2,881,859 |

One call apart, and BacFlux claims 8.7% fewer base pairs while recovering three more
curated elements. Neither caller buys recall with volume. This table has to be read
next to the recall figures, because the scorer keeps a caller's *largest overlapping*
call: a tool that emitted one enormous call per genome would score well on recovered
fraction and badly on nothing the scorer measures.

### Negative set, and the one call both callers made

Two calls each across the 32.6 Mb, on different genomes — BacFlux on *S.* Typhimurium
LT2, the EBI pipeline on *S. aureus* N315 — plus the shared call on *B. subtilis* 168.

On that shared call the two pipelines return **base-pair-identical coordinates, the
same 60 bp flanking repeat, the same MOBT relaxase and the same mating-pair type**, by
two different boundary algorithms. Because boundary refinement is the genuinely
independent layer, this agreement carries information, and what it evidences is that
ICE*Bs1* is a real element missing from the catalogue.

### Boundaries

Measured on the elements both callers found — nine chromosomal, five standalone.
Before the boundary rework, the other pipeline's boundaries were better than ours;
that is what prompted the fixes, and it is reported plainly.

| Metric | BacFlux, before | BacFlux, now | EBI pipeline |
|---|--:|--:|--:|
| Chromosomal, median \|start offset\| | 14,294 bp | **5,141 bp** | 7,427 bp |
| Chromosomal, median \|end offset\| | 34,239 bp | **14,912 bp** | 23,745 bp |
| Standalone, median recovered fraction | 0.59 | **0.94** | 0.94 |

The clearest single case: on ICE*Ec2* the two callers now land on exactly the same
interval, 27–92,263, from independent algorithms.

The fixes did not come from a better repeat finder. Of the four elements where the
other pipeline's boundary beat ours, **none** was a case where it found a repeat we
could not: two were repeats we found and then discarded through our own defects — a
ranking that put raw repeat length above tRNA anchoring, and a flank window too narrow
to reach the real *att* site — and two were candidate-extent differences with no
repeats on either side. The conservative policy was not loosened to achieve this: a
de novo repeat pair is still reported and never applied.

Both callers remain far from the curated coordinates on large chromosomal elements.
On TR2 both recover 0.055 of a 154 kb element; on IME*Ml*NZP2037-1 both recover 0.074
of 110 kb. Large chromosomal delimitation is an open problem for both, and the output
says so: those calls carry `boundary_method = none`, meaning the machinery span was
reported because no *att* pair was resolvable.

#### The 0.94 figure is not the 0.669 figure

Easy to quote wrongly, so it is spelled out. **0.94** is a like-for-like head-to-head
median over the **five** standalone ICEs *both* callers found, and the other pipeline's
0.94 is over exactly the same five. **0.669** is the pilot-wide median over **all six**
standalone ICEs BacFlux detected — the sixth is ICE*Vfl*Ind1 at 0.29, which the other
pipeline missed entirely. Detecting one extra, hard element lowers your own median.
Quote 0.669 for BacFlux alone; quote 0.94 only alongside the other caller's 0.94 on
the same five.

!!! warning "If you repeat this comparison: the `integron` term means two things"

    The EBI pipeline writes the Sequence Ontology term `integron` (SO:0000365) both
    for an ICEfinder2-lite **IME** and for a genuine **IntegronFinder integron** —
    across the 40 merged GFF files, 27 rows of the first kind and 6 of the second.
    Filtering on the type column alone silently imports unrelated integrons into the
    ICE/IME score, inflating both the call count and the base pairs called, and
    nothing errors. Filter on the `source` column first, then read the class from the
    `mobile_element_type` attribute.

## What a draft assembly costs

Every measurement above was made on closed, single-contig genomes, while most BacFlux
runs are short-read assemblies. So the same 40 genomes were cut into contigs at three
contiguities, re-annotated and re-run through the whole chain with identical
parameters: 120 assemblies, 191 calls. Breaks were placed at insertion sequences
first, then rRNA operons, then at random, because those are what actually break a real
assembly.

| | Closed | N50 150 kb | N50 50 kb | N50 20 kb |
|---|--:|--:|--:|--:|
| ICE detected (of 18) | 15 | **15** | **15** | 12 |
| IME detected (of 12) | 5 | 4 | 3 | 3 |
| Median called ÷ **true** length | 0.553 | 0.446 | 0.338 | **0.219** |
| Median called ÷ **largest surviving piece** | 0.553 | 0.574 | 0.771 | 0.721 |
| Confidence `high` | 16 (24%) | 12 (19%) | 7 (10%) | 2 (3%) |
| Confidence `low` | 7 (10%) | 22 (34%) | 42 (63%) | 39 (65%) |
| `boundary_method = none` | 30 (44%) | 49 (77%) | 55 (82%) | 54 (90%) |
| Calls on the 12 negative controls | 5 | 4 | 3 | 2 |

Three results matter more than the rest.

**Detection is flat down to 50 kb N50.** ICE recall is identical to the closed genomes
at both 150 kb and 50 kb, and only breaks at 20 kb. For calibration, a survey of 28
real BacFlux Illumina assemblies from this lab gives a median N50 of 307 kb, with 26
of the 28 at or above 150 kb — so the sweep is deliberately pessimistic.

**The length collapses, but the caller is not what collapses it.** Against the true
element, median recovered length falls from 0.55 to 0.22. Against the largest
surviving contig piece — the ceiling any per-contig caller can reach — it is flat or
better. The caller keeps recovering the same share of what is still visible; the
missing bases are the assembly's.

**Nothing is invented.** Across all 191 draft calls: no fabricated *att* sequence
(every *att* pair reported on a draft occurs at least twice in the corresponding
closed genome), no class promoted at high confidence, no new false positive on a
negative control. Calls on the negative controls went *down* under fragmentation, from
5 to 2. Contig ends removed evidence; they did not manufacture it.

What the sweep could not test: its cuts are cleaner than a real assembler's, it
contains no misassemblies, and no read-level effects — coverage variation, and the
collapse of multi-copy insertion sequences — are represented at all. Which columns to
read on a draft, and what not to conclude from them, is [its own
page](draft-assemblies.md).

## The whole deliverable, end to end

The pilots test the ICE caller. The module's actual output — one row per AMR gene with
a tier and a confidence — is exercised on *Klebsiella pneumoniae* ATCC BAA-2146
(`GCF_000364385.3`), a pan-resistant genome whose plasmids are described in the
published literature. That run reaches tiers 1, 2, 3, 5 and 6; it recovers the
chromosomal ICE at 0.946 of its curated length with a tRNA-anchored boundary and still
ends 3,138 bp inside the curated element, because a called interval is a floor. It does
not reach tier 4 — see [below](#tiers-4-and-5-and-tier-4-has-never-actually-been-assigned).
The walkthrough is the [worked example](worked-example.md); the columns are on
[Reading the output](output.md).

## Unit tests

The helper scripts are importable Python and are tested outside Snakemake:

```bash
pytest workflow/scripts
```

457 tests pass and 2 skip. 380 of them are the mobilome package's own. They test
the code, not the biology — that a repeat planted at a known offset is found, that a
confidence cap cannot be bypassed, that a malformed input is refused rather than
silently parsed. Nothing in them says a call is correct.

## What the benchmark does not show

Four limitations are invisible in the numbers above, and all four change how the
output should be read.

### A short boundary becomes a wrong tier, not just a wrong coordinate

The tables above report whether an element was found and how far its edges landed from
the curated ones. But the module's actual deliverable is a **tier per AMR gene**, and
a boundary that stops short silently demotes every gene beyond it.

Measured on a separate set of eight clinical genomes carrying ICEberg-curated
elements: of the **53 AMR genes that sit inside a curated ICE interval, only 15 (28%)
reached tier 6**, while **30 (57%) came out at tier 1, "chromosomal, intrinsic
candidate"** — the opposite conclusion. `bla`IMP-8 on `CP021851.1` is the clearest
case: it lies well inside curated Tn*6397*, but 24,007 bp past the right-hand edge the
caller drew, so it is reported as an intrinsic candidate.

**A tier 1 call in a genome that has any ICE or IME call in it is therefore weaker
evidence than a tier 1 call in a genome that has none.** Check the distance from the
gene to the nearest element call before reading "intrinsic candidate" at face value.

### Most calls report a machinery span, not the element's edges

`boundary_method = none` on 30 of the 68 calls over the closed 40-genome set (44%) —
the interval is the outermost machinery genes, not the element's edges. Of the calls
that do get an *att* pair, only a **tRNA-anchored** pair is ever applied to the
coordinates: over the 28 pilot genomes that is 5 of 36 ICE-class calls (14%). A de novo
pair is reported as a lead and never used.

That policy is deliberate and measured. On 300 randomly placed 15 kb non-element spans
of a clinical *K. pneumoniae* chromosome, with the real insertion-sequence mask
applied, **22%** returned a confident de novo "boundary"; adding a guard that counts
how often the repeat occurs elsewhere brought it to **16%**, and not applying de novo
pairs at all brings what reaches the AMR table to about **1%**. At a one-in-six error
rate a de novo pair is a reasonable lead for a human and an unacceptable basis for
silently redefining an element — widening an interval turns every gene inside it into
predicted cargo.

Read `boundary_method` before quoting any coordinate.

### About a third of calls correspond to nothing curated

Checked against every one of ICEberg's 1,677 curated entries on the 28 pilot
accessions, not just the 30 pilot elements: 17 of 51 calls (33%) without ICEscan, 21
of 63 (33%) with it, overlap no curated element; 9 self-transmissible claims and 4
high-confidence claims rest on nothing measurable in either arm.

This is an **upper bound on the error rate, not an error rate**. Every genome in both
pilots was chosen *because* it contains a curated element, so a call overlapping
nothing curated may be a genuine second element. Nothing in this benchmark can tell
the difference, and its scorer contains no code that could.

One class is unvalidated outright. **AICE** — the actinomycete class, which transfers
double-stranded DNA through a translocase instead of a relaxase — is called twice on
one *Streptomyces* genome, neither call overlaps a curated element, and neither pilot
set contains a curated AICE. Treat AICE calls as hypotheses.

### Tiers 4 and 5, and tier 4 has never actually been assigned

Tier 6 and tier 1 carry the benchmark work above. The middle of the ladder does not.

- **Tier 4** — inside a *named* transposon or integron — depends entirely on the
  opt-in TnCentral naming layer, and **no run retained on disk has ever assigned it**.
  The naming layer has produced exactly one curated hit on real data, `bla`KPC-2 inside
  Tn*7247* in a clinical *K. pneumoniae* isolate, and that gene scored **tier 6**,
  because the transposon sat on a conjugative plasmid and plasmid evidence outranks a
  transposon name. Tier 4 is reached only when a curated name is the *strongest*
  evidence available, which in practice means a chromosomal element recovered at ≥80%
  of its reference length, and no genome tested so far has produced that combination.
- **Tier 5** is exercised, but only by one of its two routes. Every tier 5 row on disk
  is a gene on a mobilisable plasmid. **No AMR gene in any run has been assigned an
  IME context**, so the "inside an IME" half of tier 5 is untested end to end — which
  matters, because the IME pilot is the weaker half of the module to begin with.

Neither gap is a reason to distrust tiers 1 and 6. Both are a reason to treat a tier 4
or tier 5 call as an unvalidated code path rather than a measured one, and to say so
if it appears in a dossier.

## On vmatch

Reviewers ask why the *att*-site search does not use vmatch, the string-matching tool
the reference ICE finders use for this step. The short answer is that BacFlux's own
search computes the same thing.

BacFlux's search finds **exact maximal repeats between the two flanking windows in
standard-library Python** — vmatch's own semantics. To test that rather than assert it,
a second implementation was built on vmatch's actual data structure (prefix-doubling
suffix array, Kasai LCP array, cross-flank maximal-match enumeration) and run against
the shipped one on every flank window the caller genuinely searched — after
insertion-sequence masking, at the real element coordinates, at the real computed
repeat-length floor:

> **52 genomes from the benchmark tree — a wider set than the 40 scored above — 246
> *att* searches, 180 distinct flank windows, 1,033 repeats, exact set equality on
> every window.** No window hit either implementation's internal cap.

Two searches over the same DNA return the same set of repeats, so vmatch would find
nothing that is not already found, and the search needs no external program at all.
The same measurement disposes of a tempting follow-up: shipping both searches and
treating their agreement as a confidence signal would be worthless, because identical
semantics are *expected* to agree exactly.

The licence position is secondary but real: bioconda's vmatch recipe declares
`license: Unknown / OTHER` and vmatch.de is unreachable, so its terms cannot be
established. With no measured benefit on one side and terms that cannot be read on the
other, it is not adopted.

!!! note "An earlier figure that should not be quoted"

    A previous version of this measurement was recorded as "35 of 35 windows, 122
    repeats". It does not reproduce under any scoping that could be reconstructed,
    because the script that produced it was never kept. The equality result itself
    reproduces, over roughly five times as many windows, and the numbers above are the
    ones to use. The script that produces them is retained with the benchmark
    artefacts and re-runs in about ten minutes.

## No code from the EBI pipeline is used here

Three of that pipeline's scripts are derived from ICEfinder2 and carry CC BY-NC-SA 4.0
headers. They were read to understand an approach, which is all that was done; what
was adopted from them is conventions — Sequence Ontology terms, an element ID format,
the discard-with-reason pattern and two numeric thresholds. Licence types for every
tool, model set and database the module can use are listed on
[Licensing](../about/licensing.md).

## Where the numbers come from

The benchmark tree — per-element tables, the scoring script, the fragmentation sweep
and the conversion audit — is kept with the benchmark run itself, outside this
repository. The method write-ups that these figures are drawn from are in the
repository and readable on GitHub: `docs/methods_ebi_comparison.md` (the head-to-head),
`docs/methods_icescan_union.md` (the two-arm model-set experiment),
`docs/methods_att_and_small_plasmids.md` (the *att* search and the vmatch measurement),
`docs/mobilome_draft_assemblies.md` (the fragmentation sweep) and
`docs/mobilome_worked_example.md` (the end-to-end run).

Call totals differ slightly between measurements — 51 and 63 in the two-arm
experiment over the 28 pilot genomes, 65 in the call-burden table over all 40, 68 in
the fragmentation sweep's closed baseline — because they were taken at different
commits over different subsets. The negative-control count moves the same way: 2 calls
in the head-to-head, 5 in the sweep's closed-genome row. Each figure is labelled with
what it was measured over; do not carry one into another's table.

## References

- Wang, M., et al. (2024). ICEberg 3.0: functional categorization and analysis of the
  integrative and conjugative elements in bacteria. *Nucleic Acids Research* 52(D1),
  D732–D737. <https://doi.org/10.1093/nar/gkad935> — the curated coordinates used as
  the truth standard, and the ICEscan models by way of ICEfinder2.
- Néron, B., et al. (2023). MacSyFinder v2: improved modelling and search engine to
  identify molecular systems in genomes. *Peer Community Journal* 3, e28.
  <https://doi.org/10.24072/pcjournal.250> — the engine both callers run.
- Coluzzi, C., Garcillán-Barcia, M. P., de la Cruz, F., & Rocha, E. P. C. (2022).
  Evolution of plasmid mobility: origin and fate of conjugative and nonconjugative
  plasmids. *Molecular Biology and Evolution* 39(6), msac115.
  <https://doi.org/10.1093/molbev/msac115> — the CONJscan model set.
- Auchtung, J. M., Lee, C. A., Monson, R. E., Lehman, A. P., & Grossman, A. D. (2005).
  Regulation of a *Bacillus subtilis* mobile genetic element by intercellular
  signaling and the global DNA damage response. *PNAS* — the ICE*Bs1* discovery paper,
  and the negative control's lesson.
- EBI Mobilome Annotation Pipeline v5.0.0, EMBL-EBI Microbiome Informatics team.
  <https://github.com/ebi-metagenomics/mobilome-annotation-pipeline>

Every tool, model set and database the module can use is listed with its citation on
[Citation and references](../about/citation.md).
