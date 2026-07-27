# Worked example: reading a mobilome report

*Destined for the MkDocs site, not the README — it is too long for a README and
works better as a page a user can read once and then keep as a reference.*

This page walks through a real BacFlux mobilome run on a genome whose answer is
already known from the literature, so every call can be checked against
something external. It doubles as the module's positive control: the isolates in
our own validation set are environmental and carry no acquired resistance, so
they only ever exercised the bottom of the mobility ladder.

## The test genome

**_Klebsiella pneumoniae_ ATCC BAA-2146** (`GCF_000364385.3`) — a complete
assembly: one chromosome plus four plasmids, carrying `blaNDM-1`, `blaCTX-M-15`,
`blaOXA-1`, `blaTEM-1`, `blaCMY-6`, `rmtC` and a long list of aminoglycoside,
sulfonamide, quinolone, macrolide and metal-resistance genes. It is pan-resistant
and its plasmids are described in the published literature, which is exactly what
a positive control needs.

Run as `mode: contigs` with `mobilome.run: true` and `phage.caller: genomad`.
Config: `BacFlux_v2_validation/kpnih1_positive_control/config.yaml`.

---

## What each output file is for

| File | What it answers |
|---|---|
| `{sample}_amr_mobility.tsv` | **The deliverable.** One row per AMR/stress gene: where it sits, what mobile element context it has, its mobility tier and confidence. |
| `{sample}_amr_mobility_audit.tsv` | Why. Every tier that was *not* raised, every confidence cap, with the evidence behind it. |
| `{sample}_replicon_calls.tsv` | Is each contig chromosome or plasmid, and if plasmid, how mobile? |
| `{sample}_is_elements.tsv` / `_is_summary.tsv` | Every insertion sequence found, plus the QC summary (how many, how many at contig ends). |
| `{sample}_ice_candidates.tsv` | Integrative elements (ICE/IME) and unbounded conjugative regions. |
| `{sample}_amrfinderplus.tsv` / `_mutations.tsv` | The raw AMRFinderPlus report, including point mutations. |
| `*_discarded.tsv` | Everything thrown away, with a reason column. |

---

## The mobility ladder, seen on real genes

Every tier the module currently implements fired on this genome, in a
biologically sensible place:

### Tier 1 — intrinsic candidate

```
gyrA_S83I        tier=1  intrinsic_candidate  conf=high
parC_S80I        tier=1  intrinsic_candidate  conf=high
ompK35_Q92Ter    tier=1  intrinsic_candidate  conf=high
ramR_Y59CfsTer13 tier=1  intrinsic_candidate  conf=high
```

These are **point mutations**, not acquired genes — chromosomal changes in the
bacterium's own `gyrA`, `parC`, `ompK35` and `ramR`. They are the cleanest
possible tier 1: nothing was acquired, so there is nothing to transfer.

They only appear at all because the isolate resolved to `Klebsiella_pneumoniae`,
one of AMRFinderPlus's 31 curated organisms, so the run passed
`--organism Klebsiella_pneumoniae`. Without that flag these rows would be
**silently absent** — see the GTDB↔NCBI organism-mapping page.

### Tier 2 — expression modulation, not mobilisation

```
blaCTX-M-15 (chromosomal copy)  tier=2  expression_modulation_not_mobilisation  conf=medium
```

An IS sits upstream, on the gene's own strand. That arrangement can supply an
outward-reading hybrid promoter and raise expression — but it does **not** make
the gene mobile. The tier deliberately stops here, and the confidence is capped
at medium because the mechanism is inferred from coordinates and strand alone;
no transcript was measured.

### Tier 3 — composite transposon

```
aadA2       tier=3  composite_mobilisable_within_cell  conf=high
qacEdelta1  tier=3  composite_mobilisable_within_cell  conf=high
sul1        tier=3  composite_mobilisable_within_cell  conf=high
```

Three genes flanked by two copies of the same IS family. `aadA2` +
`qacEdelta1` + `sul1` is the classic **class 1 integron 3′ conserved segment** —
finding them together, IS-flanked, on the chromosome is textbook.

### Tier 5 — on a plasmid

Two different flavours appeared, and the distinction matters:

```
# pHg — Platon: mobilisable (mobilization=1; inc=2)
blaOXA-1, blaTEM-1, blaCTX-M-15, aac(3)-IIe, sul2 ...
    tier=5  mobilisable_needs_helper  conf=high

# pCuAs — Platon: non-mobilisable (inc=2, no mobilisation genes)
tet(A), mph(A), the pco/sil/ars metal-resistance clusters
    tier=5  on_plasmid_typed_non_mobilisable  conf=high
```

Both are tier 5 because the gene *is* on a plasmid — acquired, which is the
regulatory question. But the second group is explicitly labelled
`on_plasmid_typed_non_mobilisable` rather than the generic
"mobilisable_needs_helper", because calling a plasmid we just typed
**non**-mobilisable "mobilisable" would contradict our own evidence in the one
column most people read.

### Tier 6 — predicted self-transmissible

```
blaNDM-1, blaCMY-6, rmtC, aac(6')-Ib3 ...   (pNDM-US)
    tier=6  predicted_self_transmissible  conf=MEDIUM
```

Note the confidence: **medium, not high.** The audit says why:

> capped at medium: tier 6 here rests on Platon's conjugation gene count
> (`conjugation=8;mobilization=1;oriT=1;inc=1`), which counts HMM hits rather
> than checking that a complete mating apparatus is present; the mating-pair
> apparatus was not verified on this contig

Eight conjugation-gene hits is a strong signal, and the tier is not lowered. But
"predicted self-transmissible" is the strongest claim this module can make about
a carbapenemase, so it is not asserted at high confidence without a verified
mating-pair apparatus. Always read the tier **and** the confidence together.

### Tier 4 — not currently reachable

Tier 4 (*inside a named unit transposon or integron cassette*) is implemented in
`colocalise.py` but nothing upstream produces a `unit_transposon`/`integron`
typed element yet: that needs the TnCentral naming cascade, which is not built.
Genes that are biologically in a named transposon therefore surface at tier 3
(pattern-based composite) or tier 5 (plasmid) instead. Not a bug, but worth
knowing before reading a report.

---

## Reading the "why not" trail

The audit file is where the module explains itself. Real example, for `blaNDM-1`:

```
tier_not_raised   flanking_is_pair_span_too_long
  1 candidate flanking pair(s) rejected. Examples:
  ...insertion_sequence-108774:110090 (IS110) / ...insertion_sequence-137590:140825 (new):
  span 32052 bp > 20000 bp limit
```

The composite-transposon detector *did* find two flanking IS copies around
`blaNDM-1` and *rejected* them, because 32 kb is beyond the configured
`max_composite_span_bp`. That is a deliberate, conservative miss rather than a
silent one — the measured span is right there in the audit, so a reader who
thinks 32 kb is defensible for their organism can raise the threshold and re-run.

---

## Two real gotchas this run exposed

### 1. `auto` decontamination can drop a genuine plasmid

The input had **five** sequences; the mobilome only ever saw **four**. The
missing one is `NZ_CP006660.1` (plasmid pMYS, 2,014 bp), removed by the
decontamination step:

```
contig          assigned_genus  action  reason
NZ_CP006660.1   Escherichia     remove  not_auto_genus:Klebsiella
```

The interesting part is *why* BlobTools called it *Escherichia*. Its single best
BLAST hit is **_Klebsiella pneumoniae_, 100% identity over the full 2,014 bp**
— the correct answer. But BacFlux runs `blobtools view --taxrule bestsum`, which
**sums bitscores per taxon** across all 95 hits:

| genus | summed bitscore |
|---|---|
| Escherichia | 58,709 |
| Klebsiella | 16,253 |
| Citrobacter | 7,696 |

A small broad-host-range plasmid hits *many* enterobacterial genomes, and
*E. coli* is vastly over-represented in `nt`. So the summed score follows
**database composition**, not biology, and outvotes the one perfect
species-level hit.

This is a real limitation of genus-consistency decontamination applied to mobile
elements: **a plasmid can look like contamination precisely because it is
mobile**. The related known case is a contaminant of the *same* genus but a
different species, which a genus-level rule cannot see at all.

`decontamination.mode: auto` remains the right default — it is correct for the
overwhelmingly common case, and every decision is written to
`contig_taxonomy_decisions.tsv` with its reason, so nothing is silent. But when
analysing plasmid-rich clinical isolates, check that file before trusting a
plasmid *absence*, and consider `mode: off` (or an explicit include list) if
small mobile replicons matter to the question being asked.

### 2. `genomad_fdr` is `NA`, and that is expected

Every row of `{sample}_plasmid_concordance.tsv` shows `genomad_fdr = NA`. This is
not a parsing failure. From geNomad's own `summary.py`:

```python
if not selected_classifier.startswith("calibrated"):
    max_fdr = None
```

geNomad only computes an FDR when its **score-calibration module** has run, and
that module is **off by default** (`--enable-score-calibration` opts in). Without
calibrated scores there is no calibrated error model, so geNomad writes `NA`
rather than a number it cannot justify — and BacFlux passes that through
unchanged instead of inventing one.

Use `plasmid_score` for confidence instead; it is populated either way. On this
genome all three surviving plasmids scored 0.9925–0.9943 and were independently
called plasmid by Platon as well:

```
contig         platon_call  platon_rds  genomad_call  genomad_score  agreement  confidence
NZ_CP006661.1  plasmid      31.8        plasmid       0.9925         both       high
NZ_CP006662.2  plasmid      20.1        plasmid       0.9937         both       high
NZ_CP006663.1  plasmid      27.1        plasmid       0.9943         both       high
```

Two methodologically independent callers agreeing at this strength is the best
outcome that table can produce.

---

## Summary of what this run demonstrates

- Point mutations are reported (tier 1) — the category a homology screen cannot see.
- IS-adjacent genes are called expression modulation, **not** mobilisation (tier 2).
- A class 1 integron's 3′ conserved segment is recovered as a composite (tier 3).
- Plasmid mobility is typed, and non-mobilisable plasmids are labelled honestly (tier 5).
- A carbapenemase on a conjugative plasmid reaches tier 6 — at *medium* confidence,
  because the mating apparatus was not verified.
- Every rejection carries its evidence in an audit file.
