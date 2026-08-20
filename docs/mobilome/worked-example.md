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

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **AMR** | antimicrobial resistance |
    | **IS** | insertion sequence — the smallest mobile element, carrying only the genes it needs to move itself |
    | **ICE** | integrative and conjugative element — carries its own conjugation machinery, so it can move itself into another cell |
    | **T4SS** | type IV secretion system — the mating apparatus that moves DNA between cells |
    | **T4CP** | type IV coupling protein — links the relaxase to the secretion system |
    | **HMM** | hidden Markov model — a statistical profile of a gene or protein family |
    | **FDR** | false discovery rate |
    | **EFSA** | European Food Safety Authority |

## Tier 1 — intrinsic candidate

```text
gyrA_S83I        tier=1  intrinsic_candidate  conf=high
parC_S80I        tier=1  intrinsic_candidate  conf=high
ompK35_Q92Ter    tier=1  intrinsic_candidate  conf=high
ramR_Y59CfsTer13 tier=1  intrinsic_candidate  conf=high
```

**Point mutations**, not acquired genes: chromosomal changes in the bacterium's own
*gyrA*, *parC*, *ompK35* and *ramR*. Nothing was acquired, so there is nothing to
transfer.

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
and strand alone: no transcript was measured, and the audit line says so.

## Tier 3 — composite transposon

```text
aadA2       tier=3  composite_mobilisable_within_cell  conf=high
qacEdelta1  tier=3  composite_mobilisable_within_cell  conf=high
sul1        tier=3  composite_mobilisable_within_cell  conf=high
```

Three genes flanked by two copies of the same IS family. `aadA2` + `qacEdelta1` + `sul1`
is the classic **class 1 integron 3′ conserved segment**.

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
typed **non**-mobilisable.

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

## Tier 4 — the naming layer, switched on

The run above was made without the TnCentral layer, so `mge_name` was `NA` on all 66
rows. Re-running the same genome with the layer on — nothing else changed — changes one
gene's answer completely.

`blaCTX-M-15` shows why. IS*Ecp1* does not merely sit beside it supplying a promoter:
IS*Ecp1* **mobilises** it, capturing the gene and moving it as a unit, and Tn*Ecp1.1* is
that unit. Reporting "expression modulation, not mobilisation" for that gene is wrong,
which is why a curated hit **overrides** the pattern-based call rather than merely
agreeing with it. The same gene, both ways:

| column | layer off | layer on |
|---|---|---|
| `mobility_tier` | **2** | **4** |
| `mge_context` | `is_adjacent` | `unit_transposon` |
| `mge_name` | `NA` | **Tn*Ecp1.1*** |
| `distance_bp` | 48 | 0 — the gene is *inside* it |
| `is_family` | IS1380 | `NA` |
| `confidence` | medium | **high** |

The match is 99.9% identity over 87% of the 3,417 bp reference element (`MF062700`). The
pattern-based read was not wrong — IS*Ecp1* is an IS1380-family element and it is 48 bp
away — but "an IS is adjacent, so expression may change" is a weaker and different claim
than "the gene sits inside a named transposon that moves it".

**Seven curated elements are named in this genome, and only one produces a tier 4 row.**
That is the precedence rule, not a shortfall:

```text
NZ_CP006659.2  chromosome        TnEcp1.1  87%   → tier 4   ← the only one
NZ_CP006661.1  plasmid, conj.    Tn7241    81%   → tier 6, replicon evidence wins
NZ_CP006661.1  plasmid, conj.    In781_p   82%   → tier 6
NZ_CP006661.1  plasmid, conj.    Tn3000   100%   → tier 6
NZ_CP006662.2  plasmid, mobil.   Tn3000    98%   → tier 5
NZ_CP006662.2  plasmid, mobil.   Tn6320    88%   → tier 5
NZ_CP006662.2  plasmid, mobil.   Tn1696.1 100%   → tier 5
```

Six of the seven sit on plasmids, and a plasmid is tested before a curated element, so
those genes score 5 or 6 and keep the name in `named_element` as supporting detail. Only
the chromosomal one is left for tier 4 to claim. A genome can be full of named
transposons and still show a single tier 4 row.

The opposite case is just as possible, and one genome later produced it: on
*Enterobacter hormaechei* a single chromosomal transposon, Tn*SMR478*, carried fifteen
resistance and stress genes, so one element produced **fifteen** tier 4 rows
([hybrid mode on five closed genomes](../methods_reference_genome_run.md)). The row
count follows the genes, not the elements — worth remembering before reading one as a
measure of the other.

**The 80% rule is doing most of the filtering.** 1,554 candidate hits were discarded to
produce those seven, and the reasons are recorded:

```text
1044  reference_coverage_below_threshold      ← the dominant one
 318  identity_below_naming_threshold
 103  interval_mostly_unaligned
  84  tncentral_hit_is_a_plain_is
   4  nested_or_overlapping_tncentral_hit
```

This is also why **tier 4 is much harder to reach on a draft**. On fragmented clinical
assemblies every Tn*Ecp1.1* candidate was discarded, with the reason recorded — *"only
12% / 49% of the 3417 bp reference element is present (needs 80%). A fragment of a
transposon is not that transposon."* Here the same element clears the bar at 87% because
this genome is closed. The refusal is correct behaviour in both cases; what differs is
the assembly. Look in `{sample}_named_elements_discarded.tsv` before concluding a genome
has no named transposon in it.

!!! note "Reproducing this"

    The naming layer was pointed at a TnCentral database fetched 2026-07-28 (533
    sequences, `fasta_sha256` starting `6a264191`). The endpoint is unversioned, so that
    digest — recorded in the database's own `PROVENANCE.txt` — is the only way to say
    which release these names came from.

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
accession** — ICEberg catalogued this ICE from this genome, so a 100% identity match
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
three, carrying no discriminating power. Use
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
- **Tier 4 is reached once, with the naming layer on** — `bla`CTX-M-15 inside Tn*Ecp1.1*,
  which the same run scores tier 2 without it. Six other curated elements are named and
  none of them produces a tier 4 row, because they sit on plasmids and the replicon is
  tested first.
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
