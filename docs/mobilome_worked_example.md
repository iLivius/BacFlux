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

### Tier 4 — `blaCTX-M-15`, inside Tn*Ecp1.1*

```
amr_gene            blaCTX-M-15
mge_context         unit_transposon
mge_name            TnEcp1.1
mobility_tier       4          named_element_mobilisable
confidence          high
identity             99.9      over 87% of the 3,417 bp reference
```

Tier 4 needs the **TnCentral naming layer** (`mobilome.tncentral.url`). Without
it nothing produces a `unit_transposon` element and this gene falls back to
**tier 2** — "an IS is adjacent, so expression may change" — which is true but
sells the situation short.

The difference is worth dwelling on, because it is the whole argument for the
naming layer. Tier 2/3 calls are **inferences** from IS positions: two copies of
one family, the right distance apart, a gene between them. That reasoning has a
famous blind spot — IS*26* forms translocatable units with its copies in *direct*
orientation, breaking the same-orientation rule the pattern depends on. A
TnCentral hit is not an inference at all: it matches an element somebody
characterised, named and deposited, whose architecture is already known.

And here the curated answer is the biologically right one. `blaCTX-M-15` is the
most prevalent ESBL in the world, and IS*Ecp1* does not merely sit beside it
providing a hybrid promoter — IS*Ecp1* **mobilises** it, capturing the gene and
moving it as a unit. Tn*Ecp1.1* is that unit. Reporting "expression modulation,
not mobilisation" would have been precisely backwards.

This is why spec §7 says a curated hit should **override** the pattern-based
call rather than merely agree with it.

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

### 2. `genomad_fdr` is `NA`, and BacFlux deliberately leaves it that way

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

**Why BacFlux does not simply switch calibration on.** It was tested on this
genome, and for a bacterial isolate it makes the output worse, not better.
Calibration converts scores into probabilities using an estimate of the sample's
composition, and it needs enough sequences to estimate that composition. From
`score_calibration.py`:

```python
if n_sequences < 1_000 and composition == "auto" and not force_auto:
    ...  # "The 'metagenome' preset will be used instead."
    composition = "metagenome"
```

Note this counts **sequences (contigs), not reads**. A finished bacterial isolate
has a handful of contigs — this genome has four — so for BacFlux's use case that
branch fires *every time*, silently substituting a preset whose assumed
composition is `[0.84 chromosome, 0.05 plasmid, 0.11 virus]`. An isolate assembly
contains essentially **no free viral sequences at all**, so the prior does not
describe the sample.

Running it anyway on this genome:

| plasmid | uncalibrated (shipped default) | calibrated |
|---|---|---|
| pCuAs | score **0.9943**, fdr `NA` | score 0.9991, fdr 0.0009 |
| pHg | score **0.9937**, fdr `NA` | score 0.9991, fdr 0.0009 |
| pNDM-US | score **0.9925**, fdr `NA` | score 0.9991, fdr 0.0009 |

Calibration collapsed three genuinely different scores into one identical value
and returned the same FDR for all three. The number looks precise to four
decimals and carries no discriminating power whatsoever — while destroying the
real variation that was there. This is the same trap the spec already flags for
PLSDB (§12.5): read the Mash *distance*, not the p-value, "which collapses to ~0
for any real match".

**So: leave calibration off, and use `plasmid_score`.** It is populated either
way and it actually discriminates. On this genome all three surviving plasmids
scored 0.9925–0.9943 and were independently called plasmid by Platon as well:

```
contig         platon_call  platon_rds  genomad_call  genomad_score  agreement  confidence
NZ_CP006661.1  plasmid      31.8        plasmid       0.9925         both       high
NZ_CP006662.2  plasmid      20.1        plasmid       0.9937         both       high
NZ_CP006663.1  plasmid      27.1        plasmid       0.9943         both       high
```

Two methodologically independent callers agreeing at this strength is the best
outcome that table can produce.

---

---

## An ICE cannot sit on a plasmid — and the classifier knows it

`{sample}_ice_candidates.tsv` reports elements on **every** replicon, chromosome
and plasmid alike: CONJscan is run over the whole proteome, so nothing restricts
it to the chromosome. On this genome it found three:

```
NZ_CP006659.2|ime-1979164:1994490                  chromosome  -> ime
NZ_CP006659.2|ice-4604073:4644558                  chromosome  -> ice   ICEKpnATCCBAA-2146-1-like
NZ_CP006661.1|conjugative_region-57877:90473       PLASMID     -> conjugative_region
```

That name comes from the **ICEberg naming layer** (`mobilome.iceberg.urls`), and
the audit line behind it is a small lesson in reading this module:

> named `ICEKpnATCCBAA-2146-1-like` (100.00% identity, covering 100% of our
> interval but only **70%** of the 58,048 bp curated element, accession
> CP006659.2). **30 other
> curated element(s) fit about as well**, so treat the exact name as one of a
> near-identical group. **NOTE** the curated element is 58,048 bp while our
> interval is 40,486 bp: our boundaries are a floor, not the element's true ends.

Three things to take from that. First, `CP006659.2` *is* KPNIH1's own chromosome
accession — ICEberg catalogued this ICE **from this genome**, so 100% identity is
a self-match and not independent confirmation. Second, ICEs of one species are
near-identical across strains, so a single real element matches dozens of
entries; the name is a group label, not a unique identification. Third, and most
useful: the curated record is **17.5 kb longer** than our call. With
`boundary_method = none` our interval is the machinery span — a floor — and here
is exactly how much we are missing. That shortfall is also why the name carries
a **`-like`** suffix: we have 70% of the curated element, not all of it, and the
suffix is the module refusing to claim the bare name for a partial match.

> **A worked example of how this table gets things wrong, and how it was caught.**
> The middle row used to read `conjugative_region-4610740:4644558` — the same
> machinery, but no integrase found, so it fell into the last row of the table
> below instead of the first. The integrase was there all along: Bakta annotates
> it `DNA integration/recombination/inversion protein` at 4,604,073–4,605,020,
> 5.7 kb away and comfortably inside the 15 kb clustering window, and the product
> regex simply did not match that wording. The element was being under-called on
> the very genome used to validate the module. Worth remembering when reading any
> `conjugative_region` row: the class turns on **one** annotation, and an
> annotation is a string match.

The plasmid one is typed `conjugative_region`, never `ice`, and that is a
definitional point rather than a threshold: **ICE** stands for *Integrative and
Conjugative Element* (spec §2.4 — "integrates into the **chromosome** and encodes
conjugation machinery"). Something that is already its own replicon has nothing
to integrate into; a conjugative element on a plasmid is simply a **conjugative
plasmid**. The classifier enforces this through the integrase requirement:

| integrase | relaxase | MPF | class |
|:-:|:-:|:-:|---|
| ✓ | ✓ | ✓ | **ICE** — predicted self-transmissible |
| ✓ | ✓ | ✗ | **IME** — mobilisable, needs a helper |
| ✓ | ✗ | ✗ | CIME / genomic island — passive |
| ✗ | ✓ | ✓ | **conjugative region** — *report it, do not call it an ICE* |

pNDM-US has a relaxase (MOBH), a coupling protein and a full typed F-type mating
apparatus, but **no integrase**, so it lands in the last row — correctly.

### Platon counts genes; CONJscan types the machinery

That same plasmid row says `machinery_intact = TRUE`, `mpf_typed_system = TRUE`,
`confidence = high`. geNomad independently agrees, listing `MOBH`, `t4cp1`,
`virb4` and a full `F_tra*` set in its own `conjugation_genes` column.

The two tools answer genuinely different questions about the same replicon:

- **Platon** *counts* hits to conjugation-related genes. Cheap, always available,
  but a count is not proof of a working machine.
- **CONJscan** asks whether a **complete, typed** mating-pair system is present,
  and `conjscan_to_ice.py` writes that verdict out — for every replicon, plasmids
  included, because CONJscan runs over the whole proteome.

So Platon sets the **tier** (it is the replicon-level mobility call) and CONJscan
**corroborates or contradicts** it, which moves the **confidence** and gets named
in the row. Three outcomes:

| Platon | CONJscan | result |
|---|---|---|
| conjugative | complete typed system | tier 6, **high** — audit `plasmid_conjugation_machinery_verified` |
| conjugative | nothing / relaxase only / degraded | tier 6, **medium** — audit `plasmid_conjugation_from_hit_counts_only` |
| *not* conjugative | complete typed system | tier **unchanged**, confidence dropped — audit `conjscan_typed_system_but_platon_did_not` |

That third row is deliberate: a disagreement is **not** silently promoted to
"predicted self-transmissible". The tier follows Platon's replicon-level call,
the conflict is stated plainly, and a human decides.

On this genome the effect is visible on the one call that matters most:

```
blaNDM-1   before:  tier=6  predicted_self_transmissible  conf=medium
           after:   tier=6  predicted_self_transmissible  conf=high
                    relaxase=MOBH  mpf=F  intact=yes
                    machinery_source=same_replicon:NZ_CP006661.1|conjugative_region-57877:90473
```

and it is discriminating rather than blanket — pHg (relaxase `MOBC` only) and
pCuAs (no machinery) are untouched, still `NA` in the machinery columns.

### The machinery columns

`relaxase_type`, `mpf_type`, `machinery_intact` and `boundary_method` / `attL` /
`attR` are in the deliverable (spec §9), so the evidence behind a transferability
claim is readable without opening the ICE table. `machinery_source` says where it
came from, and the distinction matters:

- `containing_element:<id>` — the gene sits **inside** that ICE/IME.
- `same_replicon:<id>` — the machinery is **elsewhere on the same plasmid**. Still
  the right evidence (a mating apparatus makes the whole replicon transferable,
  and with it every gene on it), but it is not wrapped around this gene, and the
  row says so rather than implying containment.

### `boundary_method`, and why it is not folded into the confidence

The att-site search (spec §8 Phase 3) **does** run, in every mode — the constraint
is assembly contiguity, not the sequencer. `boundary_method` says how the
element's ends were established, and it takes three values:

| value | meaning | does it move the element? |
|---|---|---|
| `tRNA` | attL/attR found at an annotated tRNA 3′ end — the site integrases actually target | **yes** |
| `denovo` | a bracketing direct repeat was found, with no tRNA behind it | **no** — reported only |
| `none` | no att pair; the interval is the machinery span | no |

**A de novo repeat never moves an element, and the reason is worth knowing.**
Measured on this genome's chromosome, 300 randomly placed 15 kb non-ICE spans
produced a confident de novo "boundary" **22%** of the time. Real chromosomes are
full of rRNA operons, REP elements and paralogues, so an exact 18–25 bp direct
repeat between two 30 kb windows is ordinary rather than remarkable. Requiring the
repeat to occur exactly twice — an integration scar is created once; an rRNA
operon exists seven times — cut that to 16%, and refusing to *act* on a de novo
repeat at all cuts the rate that can reach the AMR table to **1%**.

That restraint matters because widening an element makes every gene inside it
*cargo*. A fabricated 50 kb boundary would turn a `gyrA` point mutation — the
textbook intrinsic, non-transferable determinant — into "predicted
self-transmissible". So an unresolved boundary always leaves the interval at the
machinery span: a floor, never an invention.

**Confidence answers a different question.** `confidence` is about *is this an
ICE?* — anchor classes, intact machinery, one contig. `boundary_method` is about
*where does it stop?* On a fragmented short-read assembly the second usually has
no answer, so by default the two are read side by side rather than multiplied
together. Set `mobilome.require_trna_boundary_for_high: true` to apply the strict
spec §8 Phase 6 reading, which demands a tRNA-anchored boundary before any element
may be called `high` — sensible on **closed long-read assemblies**, where an
unresolved boundary is a warning sign rather than the norm.

---

## Summary of what this run demonstrates

- Point mutations are reported (tier 1) — the category a homology screen cannot see.
- IS-adjacent genes are called expression modulation, **not** mobilisation (tier 2).
- A class 1 integron's 3′ conserved segment is recovered as a composite (tier 3).
- Plasmid mobility is typed, and non-mobilisable plasmids are labelled honestly (tier 5).
- A curated TnCentral hit puts blaCTX-M-15 in TnEcp1.1 (tier 4), overriding the
  weaker IS-adjacency inference — and ICEberg names the chromosomal ICE.
- A carbapenemase on a conjugative plasmid reaches tier 6 — at *medium* confidence,
  because the mating apparatus was not verified.
- Every rejection carries its evidence in an audit file.
