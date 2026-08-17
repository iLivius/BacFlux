# Worked example

One genome, end to end, with the real rows. The isolate was chosen because its answer is
already known from the literature, so every call can be checked against something
external.

**_Klebsiella pneumoniae_ ATCC BAA-2146** (`GCF_000364385.3`) — a complete assembly, one
chromosome and four plasmids, carrying `blaNDM-1`, `blaCTX-M-15`, `blaOXA-1`, `blaTEM-1`,
`blaCMY-6`, `rmtC` and a long list of aminoglycoside, sulfonamide, quinolone, macrolide
and metal-resistance genes. Run as `mode: contigs` with `mobilome.run: true` and
`phage.caller: genomad`.

The isolates in BacFlux's own validation set are environmental and carry no acquired
resistance, so they only ever exercise the bottom of the ladder. This genome exercises
almost all of it.

## Tier 1 — intrinsic candidate

```text
gyrA_S83I        tier=1  intrinsic_candidate  conf=high
parC_S80I        tier=1  intrinsic_candidate  conf=high
ompK35_Q92Ter    tier=1  intrinsic_candidate  conf=high
ramR_Y59CfsTer13 tier=1  intrinsic_candidate  conf=high
```

**Point mutations**, not acquired genes: chromosomal changes in the bacterium's own
*gyrA*, *parC*, *ompK35* and *ramR*. Nothing was acquired, so there is nothing to
transfer — the cleanest possible tier 1.

They appear at all only because the isolate resolved to `Klebsiella_pneumoniae`, one of
AMRFinderPlus's curated organisms, so the run passed
`--organism Klebsiella_pneumoniae`. Without that flag these rows would be **silently
absent** ([Turning it on](enabling.md)).

## Tier 2 — expression modulation, not mobilisation

```text
blaCTX-M-15 (chromosomal copy)  tier=2  expression_modulation_not_mobilisation  conf=medium
```

An IS sits upstream, on the gene's own strand. That can supply an outward-reading hybrid
promoter and raise expression — it does not make the gene mobile. The tier stops here,
and the confidence is capped at medium because the mechanism is inferred from coordinates
and strand alone: no transcript was measured, and the audit line says exactly that.

## Tier 3 — composite transposon

```text
aadA2       tier=3  composite_mobilisable_within_cell  conf=high
qacEdelta1  tier=3  composite_mobilisable_within_cell  conf=high
sul1        tier=3  composite_mobilisable_within_cell  conf=high
```

Three genes flanked by two copies of the same IS family. `aadA2` + `qacEdelta1` + `sul1`
is the classic **class 1 integron 3′ conserved segment** — finding them together,
IS-flanked, on the chromosome is textbook.

## Tier 5 — on a plasmid, two flavours

```text
# pHg — Platon: mobilisable (mobilization=1; inc=2)
blaOXA-1, blaTEM-1, blaCTX-M-15, aac(3)-IIe, sul2 …
    tier=5  mobilisable_needs_helper  conf=high

# pCuAs — Platon: non-mobilisable (inc=2, no mobilisation genes)
tet(A), mph(A), the pco/sil/ars metal-resistance clusters
    tier=5  on_plasmid_typed_non_mobilisable  conf=high
```

Both are tier 5 because the gene *is* on a plasmid — acquired, which is the regulatory
question. The second group carries a different label because tier 5's generic wording,
`mobilisable_needs_helper`, would contradict our own evidence about a plasmid we had just
typed **non**-mobilisable, in the one column most people read.

## Tier 6 — predicted self-transmissible

```text
blaNDM-1, blaCMY-6, rmtC, aac(6')-Ib3 …   (pNDM-US)
    tier=6  predicted_self_transmissible  conf=high
    relaxase=MOBH  mpf=F  intact=yes
    machinery_source=same_replicon:NZ_CP006661.1|conjugative_region-57877:90473
```

Two independent methods agree here, and that is what the confidence records. Platon typed
the replicon conjugative from its gene counts; CONJscan independently typed a **complete**
F-type mating apparatus with a MOBH relaxase on the same contig, so the audit carries
`plasmid_conjugation_machinery_verified` and tier 6 is not capped.

Without that corroboration the same row comes out at **medium**, with the reason spelled
out: tier 6 would rest on a count of HMM hits rather than on a verified machine. The
effect is discriminating rather than blanket — pHg (relaxase MOBC only) and pCuAs (no
machinery) are untouched and still `NA` in the machinery columns.

`machinery_source` says `same_replicon`, not `containing_element`: the mating apparatus is
elsewhere on the plasmid, not wrapped around this gene. That is still the right evidence —
a mating apparatus makes the whole replicon transferable, and with it every gene on it —
and the column says so rather than implying containment.

## Tier 4 — what the naming layer is for, and why this run does not reach it

`mge_name` is `NA` on **all 66 rows**: this run did not have the TnCentral layer switched
on. Nothing else produces a `unit_transposon`, so a gene a curated transposon would have
explained falls back to tier 2 — "an IS is adjacent, so expression may change" — which is
true and can sell the situation short.

`blaCTX-M-15` is the standing example. IS*Ecp1* does not merely sit beside it supplying a
promoter: IS*Ecp1* **mobilises** it, capturing the gene and moving it as a unit, and
Tn*Ecp1.1* is that unit. Where a run can establish that, reporting "expression
modulation, not mobilisation" is precisely backwards — which is why a curated hit
**overrides** the pattern-based call rather than merely agreeing with it.

But note what happened on the clinical runs that *did* have the layer on: every
Tn*Ecp1.1* candidate was **discarded**, with the reason recorded — *"only 12% / 49% of the
3417 bp reference element is present (needs 80%). A fragment of a transposon is not that
transposon."* The refusal is the correct behaviour, and it is also why tier 4 is so hard
to reach on a fragmented assembly. Look in `{sample}_named_elements_discarded.tsv` before
concluding a genome has no named transposon in it.

**No run kept on disk contains a tier 4 row.** The naming layer has matched a curated
element on real data exactly once — `bla`KPC-2 inside Tn*7247* — and that gene scored
tier 6, because the transposon sat on a conjugative plasmid and the replicon evidence is
tested first ([The mobility ladder](mobility-ladder.md)).

## The chromosomal ICE, scored against its curation

ICEberg curates an element at exactly this locus, so the call can be scored:

| | |
|---|---|
| curated element | `ICEKpnATCCBAA-2146-1`, CP006659.2:4,603,840–4,661,887 (58,048 bp) |
| our call | 4,603,807–4,658,749 (**54,943 bp**) |
| start / end offset | **−33 bp** / **−3,138 bp** |
| fraction of the curated element recovered | **0.946** |
| `boundary_method` | **`tRNA`**, a 43 bp repeat at tRNA-Phe(gaa) |
| machinery | integrase + relaxase + T4CP + T4SS, all intact |
| confidence | **high** |

Even at 0.946 the interval is a **floor**: it ends 3,138 bp inside the curated element,
and any AMR gene in that last 3 kb is scored as though it were outside the ICE. That is
the general failure mode, and it is much larger on other genomes — see
[Validation](validation.md), where 57% of AMR genes inside curated ICE intervals come out
at tier 1 precisely because the called edge stopped short.

Two cautions come with that row. `CP006659.2` is **this isolate's own chromosome
accession** — ICEberg catalogued this ICE from this very genome, so a 100% identity match
is a self-match, not independent confirmation. And ICEs of one species are near-identical
across strains, so one real element matches dozens of curated entries (30 others fit about
as well here). A name from that layer is a **group label, not a unique identification**,
which is what the `-like` suffix is for.

## Reading the "why not" trail

The audit file is where the module explains itself. For `blaNDM-1`:

```text
tier_not_raised   flanking_is_pair_span_too_long
  1 candidate flanking pair(s) rejected. Examples:
  ...insertion_sequence-108774:110090 (IS110) / ...insertion_sequence-137590:140825 (new):
  span 32052 bp > 20000 bp limit
```

The composite detector found two flanking IS copies and rejected them: 32 kb is beyond
`mobilome.max_composite_span_bp`. A conservative miss rather than a silent one — the
measured span is in the file.

## Two things this run exposed

**A genuine plasmid can be dropped before the module ever sees it.** The input had five
sequences; the mobilome saw four. `NZ_CP006660.1` (plasmid pMYS, 2,014 bp) was removed by
the contamination screen as *Escherichia*, even though its single best BLAST hit was
*K. pneumoniae* at 100% identity over the full length. BlobTools sums bitscores per taxon
across all 95 hits, and *E. coli* is hugely over-represented in `nt`, so *Escherichia*
summed to 58,709 against *Klebsiella*'s 16,253. The assignment followed database
composition, not biology — and a broad-host-range plasmid can be discarded *precisely
because* it is mobile. Check `contig_taxonomy_decisions.tsv` before trusting a plasmid's
absence. See [Decontamination](../analysis/decontamination.md).

**`genomad_fdr` is `NA`, deliberately.** geNomad computes an FDR only when its
score-calibration module has run, and that module is off by default. Switching it on is
worse rather than better for an isolate: calibration needs enough sequences to estimate
the sample's composition, and below 1,000 *contigs* it silently substitutes a metagenome
preset whose assumed composition is 11% virus — which does not describe a finished
bacterial genome. Run on this genome, it collapsed three genuinely different plasmid
scores (0.9925, 0.9937, 0.9943) into one identical 0.9991 with the same FDR for all
three: precise to four decimals, and carrying no discriminating power at all. Use
`plasmid_score`, which is populated either way and does discriminate. See
[Plasmids](../analysis/plasmids.md).

## What this run demonstrates

- Point mutations are reported (tier 1) — the category a homology screen cannot see.
- IS-adjacent genes are called expression modulation, **not** mobilisation (tier 2).
- A class 1 integron's 3′ conserved segment is recovered as a composite (tier 3).
- Plasmid mobility is typed, and a non-mobilisable plasmid is labelled honestly (tier 5).
- `bla`NDM-1, on a conjugative plasmid, reaches tier 6 at high confidence — because
  CONJscan typed the mating apparatus rather than merely agreeing with a gene count.
- The chromosomal ICE is recovered at 0.946 of its curated length with a tRNA-anchored
  boundary — and still ends 3.1 kb short, because a called interval is a floor.
- **Tier 4 is not reached, and that is the honest result**, not an omission from the
  walkthrough.
- Every rejection carries its evidence in an audit file.

## What none of it establishes

Every mobility statement here is a **prediction made from sequence**. "Predicted
self-transmissible" means an element carries conjugation machinery that looks complete
and intact; it does not mean the element has been shown to move. The confirmatory
experiment is a filter or broth mating assay, and no output of this module substitutes
for one.

Nor is any of it "EFSA-compliant", a claim BacFlux does not make. What the module
produces is supporting evidence for the intrinsic-versus-acquired judgement. The
judgement remains the analyst's.

The longer write-up, with the full file inventory and the corrections made to it along
the way, is
[`mobilome_worked_example.md`](https://github.com/iLivius/BacFlux/blob/main/docs/mobilome_worked_example.md);
example output files are in
[`docs/validation/`](https://github.com/iLivius/BacFlux/tree/main/docs/validation).
