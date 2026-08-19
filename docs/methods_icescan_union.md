# Method reference: running CONJScan and ICEscan together

BacFlux can search for conjugation machinery with two different model sets. This page
explains why it runs **both and merges the results** rather than picking one, what the
second set genuinely adds, and which parts of it are deliberately distrusted.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **CONJScan** | the default model set, distributed with MacSyFinder. Conservative and well validated |
    | **ICEscan** | the optional second set, a fork of CONJScan taken from ICEfinder2. Off by default; adds classes CONJScan cannot see |
    | **MacSyFinder** | the engine both run on: it looks for a quorum of required genes, not single hits |
    | **ICE** | integrative and conjugative element — moves itself between cells |
    | **IME** | integrative mobilisable element — needs a helper element to move |
    | **AICE** | actinomycete ICE, which replicates rather than conjugates in the usual way |
    | **T4SS** | type IV secretion system, the mating apparatus |
    | ***att* site** | the short repeat marking an integrated element's ends |
    | **MGE** | mobile genetic element, the general term |
    | **IS** | insertion sequence, the smallest kind of mobile element |
    | **HMM** | hidden Markov model, a statistical profile used to recognise a protein family |

Measured on the 28-genome Phase 7 benchmark
(`<validation-root>/phase7_benchmark`),
comparing two arms that differ **only** by whether the ICE caller was given
`--icescan-tsv`.

> ### ⚠ Read this before quoting any coordinate from this page
>
> This is a **controlled two-arm experiment**, and it answers exactly one
> question: *what does adding ICEscan change?* Both arms were run on
> **2026-07-30**, before commit `4a93d89` reworked the *att*-site search.
>
> **The arm comparison is still valid** — both arms used the same att search, so
> the difference between them is unaffected, and that difference is what this
> document is for.
>
> **The absolute per-element coordinates and ratios in §8.1 and §8.2 are
> superseded.** `4a93d89` moved many of them substantially: ICE*Ec2*, for
> instance, is listed below at 55,055 bp (ratio 0.59) and is now called at
> 92,237 bp (0.98); SPI-7 is listed at 83,403 bp and is now 133,582 bp against a
> curated 133,500. Current per-element numbers live in `results_final.tsv` and
> `ime_results_final.tsv` in the benchmark tree, and the current headline figures
> are in `docs/mobilome/validation.md` and in `methods_ebi_comparison.md`.
>
> Treat them as a dated snapshot, not as the module's current output.
>
> **§8.3, §9 and the table in §10 were measured at commit `1651e6f`**
> (after the att-search rework, the confidence-cap fix, the loner fix that
> stopped one MacSyFinder system inventing a second machinery cluster, and the
> dead-code removal), by the same two-arm method: both arms
> run over the same 28 genomes, the only variable being `--icescan-tsv`.
> Three things differ from the July-30 figures:
>
> - The census counts moved (§8.3): `ime` is 21 in the union arm, not 22.
> - **Cost 3 in §9 — the nested double-report on CP011419.1 — no longer
>   happens.** It was fixed by a later commit (the loner fix) that the July-30
>   snapshot predates.
> - The fraction of calls resolved by a tRNA-anchored att pair is now much lower
>   in ABSOLUTE terms than the July-30 figure (12–13% against the original
>   34–35%). This is **not** a regression introduced by ICEscan — restricted to
>   `ice`-class calls only, the rate is identical in both arms today (5/36, 14%,
>   both control and union), so the union-vs-control comparison this document
>   exists to make is unaffected. What changed is some combination of the
>   att-search rework and the loner fix, applying equally to both arms. The
>   original per-genome run artefacts from July 30 no longer exist on disk, so
>   the exact mechanism cannot be reconstructed further.

**Short version.** ICEscan is not a second opinion — it is a fork of the tool we
already run. We add it *alongside* CONJScan rather than in place of it, take only
two things from it (integrase hits and the IME/AICE element classes), and
deliberately ignore four of its integrase profiles because they are not element
integrases at all. The measured effect is: integrative and conjugative elements
(ICEs) unchanged, integrative *mobilisable* elements (IMEs) meaningfully better,
boundaries barely moved.

---

## 1. What ICEscan is, and the correction that matters

ICEscan is the machinery-model set distributed with **ICEfinder2** and used by
the **EBI mobilome-annotation-pipeline**. It was adopted here after reading
their work — the approach is theirs, and the credit is theirs (§10).

It is natural to assume that adding ICEscan to a workflow that already runs
CONJScan gives you two independent opinions. **It does not.** ICEscan is a
*fork of CONJScan*, by the same authors at Institut Pasteur, one minor version
behind the release BacFlux installs.

The evidence is in the model set's own metadata — ICEscan's `metadata.yml`
does not describe itself as ICEscan at all:

```yaml
maintainer:
  name: Charles Coluzzi Eduardo Rocha
  email: ccoluzzi@pasteur.fr erocha@pasteur.fr
short_desc: CONJScan - Models for detection of conjugative and mobilisable elements.
vers: 2.0.1
doc: https://github.com/macsy-models/CONJScan
license: CC BY-NC-SA 4.0
copyright: 2018-2022, Institut Pasteur, CNRS
```

Its `README.md` is likewise the CONJScan README verbatim, down to the
`macsydata install CONJScan` instruction. MacSyFinder's own run log calls the
package `ICEscan-2.0.1`; BacFlux installs `CONJScan 2.1.0`.

**Why this matters practically.** Two forks of one model set share their
profiles, their scoring and most of their blind spots. Where they agree, that
agreement is largely structural and is *not* independent corroboration. Where
they disagree, the disagreement is informative — and it is only the
disagreements we use.

---

## 2. Exactly what ICEscan adds and removes

Measured by diffing the two installed model directories.

### 2.1 System models

| | CONJScan 2.1.0 | ICEscan 2.0.1 |
|---|---|---|
| `Chromosome/T4SS_type{B,C,F,FA,FATA,G,I,T}` | 8 | 8 |
| `Chromosome/dCONJ_type*` (decayed machinery) | 8 | **0 — removed** |
| `Chromosome/MOB` (relaxase-only, mobilisable) | 1 | **0 — removed** |
| `Plasmids/*` (the entire plasmid model set) | 17 | **0 — removed** |
| `Chromosome/IME` | **absent** | **1 — added** |
| `Chromosome/AICE` | **absent** | **1 — added** |
| **total definition files** | **34** | **10** |

### 2.2 Profile HMMs — 21 added, 1 removed

ICEscan carries 145 profiles against CONJScan's 125.

**Added (21).** Eight integrase/recombinase families —
`Phage_integrase`, `Recombinase`, `UPF0236`, `PB001819`, `rve`, `TIGR02224`,
`TIGR02225`, `TIGR02249`; the actinomycete AICE machinery —
`FtsK_SpoIIIE`, `RepSAv2`, `DUF3631`, `Prim-Pol`; Gram-positive and
IME-specific relaxase families — `Relaxase_firmi_{MOBL,Rep_2,Viral_Rep_A,Viral_Rep_B1,Viral_Rep_B2}`,
`Relaxase_PHA_IME_{A1,B}`, `Relaxase_profile_MOBT`; plus `T4SS_MOBL`.

**Removed (1).** `T4SS_MOBM`.

The eight integrase profiles are the substantive addition: CONJScan 2.1.0
carries **no integrase profile at all**, because conjugation machinery is its
subject and integration is not.

### 2.3 The mating-apparatus quorum is stricter in ICEscan

"Quorum" here means how many machinery genes MacSyFinder must find, and how
close together, before it will call a type IV secretion system (T4SS — the
apparatus that physically pushes DNA into the recipient cell).

| | CONJScan 2.1.0 chromosome T4SS | ICEscan 2.0.1 T4SS |
|---|---|---|
| `min_mandatory_genes_required` | 3 | **4** |
| `min_genes_required` | 4–6 (varies by type) | **8** (all types) |
| `inter_gene_max_space` | 60 genes | **20 genes** |

ICEscan is stricter on every axis. On a fragmented or a decayed element that is
the difference between a call and silence — which is the core reason a swap
loses elements.

---

## 3. Why we union rather than swap

A swap would trade CONJScan's `MOB`, `dCONJ` and plasmid models, and its looser
T4SS quorum, for ICEscan's `IME` and `AICE` models. Measured per genome, that
trade is bad:

| genome | curated element | CONJScan hits | ICEscan hits |
|---|---|---|---|
| `AF261825` | **SGI1** | 2 | **0 — genome emptied** |
| `GU725392` | **ICEEc2** | 24 | **3** |
| `AL513382` | SPI-7 | 24 | 3 |
| `AE009948` | — | 34 | 27 |
| `NC_004668_1` | ICE_EfalV583_tRNALys | 47 | 44 |

SGI1's genome goes to **zero** machinery hits under ICEscan alone: SGI1 is
mobilised by a helper plasmid and its signal in CONJScan is a `MOB` relaxase —
a model ICEscan deleted. ICEEc2 drops from 24 supporting hits to 3.

> **Not re-verifiable.** Earlier notes record a swap arm in which corroboration
> across the two model sets fell from **88% to 68%**. The per-element artefacts
> of that swap experiment are no longer on disk (only the `head_noice`,
> `new_noice` and `union` arms survive), so those two percentages could not be
> reproduced here and should not be quoted in a paper. The per-genome collapse
> in the table above *is* directly re-measurable and makes the same point.

**The union rule.** Both model sets run at `--coverage-profile 0.5` over the
same Bakta proteins; the two `best_solution.tsv` tables are merged before the
ICE caller sees them. From ICEscan we take **only**:

1. **integrase anchors** — from four trusted profiles (§4), and
2. **the IME and AICE element classes**.

We never take ICEscan's T4SS quorum, and we never use its system spans as
element boundaries.

---

## 4. Four integrase profiles we deliberately do not trust

ICEscan's `IME.xml` lists eight profiles as interchangeable integrases. Four of
them are not element integrases. Each is excluded, and each exclusion is written
to the audit TSV with its reason.

| profile | what the protein actually is | why it must not anchor an element |
|---|---|---|
| `TIGR02249` | **IntI1**, the class-1 **integron** integrase | An integron is cargo *inside* elements, not the element. Measured casualty on `CP042858.1`: an att-bounded **32,103 bp** tier-6 ICE became a **103,303 bp unbounded** one, at unchanged `high` confidence, anchored on IntI1. Fires 4× on this set, every time as IntI1. |
| `TIGR02224` | **XerC**, chromosomal *dif*-site recombinase | Housekeeping. Resolves chromosome dimers at *dif*; has nothing to do with mobile elements. A XerC was measured being attached **43,639 bp** from a 945 bp relaxase cluster, producing an element 6× its annotated size. |
| `TIGR02225` | **XerD**, the other half of the same *dif*-site pair | As XerC. |
| `rve` | the **rve/DDE catalytic domain**, shared by IS transposases | Not an integrase family — a domain. Of its used hits, the overwhelming majority were proteins Bakta independently annotates as transposases (IS3 family ×16, IS6 IS15DIV ×14, IS6 IS1216E ×8). Anchoring on one turns any insertion sequence into a spurious element. |

**We diverge from upstream on `rve` on purpose.** ICEscan's `IME.xml` lists
`rve` as an exchangeable of `Phage_integrase`, so a reader comparing BacFlux to
ICEscan will find a real disagreement here. It is deliberate: `rve` is the DDE
domain shared with IS transposases, which BacFlux detects separately, and
admitting it as an integrase would double-count IS elements as IMEs. The audit
text says so in as many words.

**Trusted (4):** `Phage_integrase`, `Recombinase`, `UPF0236`, `PB001819`.

### 4.1 A subtlety in how the exclusion is applied

MacSyFinder's `best_solution.tsv` has two gene columns. `gene_name` is the
profile that actually matched; `hit_gene_ref` is the top-level model gene it
stands in for. For a `Recombinase` hit these read:

```
gene_name = Recombinase        hit_gene_ref = Phage_integrase
```

Every excluded profile is an exchangeable of `Phage_integrase`, so **filtering
on `hit_gene_ref` would let all four back in.** The caller gates on `gene_name`.

### 4.2 What was actually discarded, and how we know

Verifying the exclusions by grepping the element table for profile names is
**vacuous** — that table never records profile names, so the grep returns zero
whether or not the fixes exist. The real evidence is the audit TSV, which
carries 15 `row_skipped | untrusted_integrase_profile` rows totalling:

| profile | hits discarded |
|---|---|
| `rve` | **103** |
| `TIGR02249` (IntI1) | 5 |
| `TIGR02224` (XerC) | 1 |
| `TIGR02225` (XerD) | 1 |

Raw availability confirms the damage was real rather than theoretical: the
un-filtered ICEscan output for this set contains 71 `rve` hits and 4
`TIGR02249` hits in `best_solution.tsv` alone, none of which reaches an element.

**Regression check, `CP042858.1`** — byte-identical in both arms:

```
4559721-4591823  len=32103  ice  bnd=denovo  conf=high
4911034-4988140  len=77107  ice  bnd=denovo  conf=high
```

The att-bounded 32,103 bp element survives at high confidence, and no 103,303 bp
element exists anywhere in the set.

---

## 5. What ICEscan does *not* solve

**It emits no coordinates.** MacSyFinder reports gene *ordinals*, never base
pairs. A real row from the Tn*4451* detection:

```
replicon  gene_name     hit_pos  model_fqn                locus_num
U15027    Recombinase   2        ICEscan/Chromosome/IME   1
U15027    T4SS_MOBV     6        ICEscan/Chromosome/IME   1
```

`hit_pos` 2 and 6 mean "the 2nd and 6th protein in the file". **Every base pair
BacFlux prints comes from its own join against the Bakta GFF3.** The clustering,
the coordinate work and `attach_nearby_integrases` cannot be delegated to
ICEscan and are not.

**Its IME systems sprawl and are not loci.** In `IME.xml` both mandatory genes
(`T4SS_MOBV` and `Phage_integrase`) carry `loner="1"`, which exempts them from
MacSyFinder's clustering entirely. Across the 28 benchmark genomes, in the
`best_solution.tsv` files the caller actually consumes:

| | count |
|---|---|
| ICEscan IME "systems" | 54 |
| …spanning ≤ 11 genes (usable as a locus) | 22 |
| …spanning > 1,000 genes | 12 |
| widest single "IME" | **8,784 genes** |

That widest one is on `NC_013929` (*Streptomyces scabiei*) and covers essentially
the **entire 10,148,695 bp chromosome as one "IME"**. BacFlux's own clustering is
what turns these into usable intervals.

> Earlier notes give these as 71 / 27 / 17. Neither `best_solution.tsv`
> (54 / 22 / 12) nor `all_systems.tsv` (97 / 33 / 24) reproduces those figures,
> so the counts above are the ones to quote. The headline example — 8,784 genes,
> *S. scabiei* — is identical in every file and is confirmed.

**It delimits nothing.** ICEscan performs no att-site search. Every boundary
BacFlux reports comes from its own att search, IS masking and tRNA scoring.
Bounded fraction moves from 17/50 (34%) to 22/63 (35%) — i.e. essentially not at all.

**The Bakta product-text regex stays a first-class source.** Twelve of the 30
curated pilot elements have an integrase found *only* by
`INTEGRASE_PRODUCT_PATTERN` matching Bakta's product text, with no usable
ICEscan hit. It is not a fallback.

---

## 6. The coverage-threshold confound, stated plainly

MacSyFinder's `--coverage-profile` sets how much of an HMM profile a protein
must cover to count as a hit. Lowering it finds more, and more marginal, hits.
This is a genuine confound and it is easy to credit the wrong cause.

**We adopt 0.5 for both model sets.** All 28 on-disk runs were verified:

```
28 × coverage_profile = 0.5   (ICEscan-2.0.1)
28 × coverage_profile = 0.5   (CONJScan-2.1.0)
```

**The IME gain is 100% model set, not threshold.** Tn*4451* is detected through
`Recombinase` + `T4SS_MOBV` as an `IME` system. `Recombinase.hmm` is one of the
21 profiles ICEscan adds and CONJScan does not have, and CONJScan has no `IME`
model at all — on the same genome CONJScan finds a single `MOB` relaxase hit and
no integrase. **No coverage threshold can make CONJScan find a profile it does
not ship.** The same holds for `IME_CdiR20291_ND`.

**The ICE improvements seen at 0.3 are 100% threshold.** They reproduce with
plain CONJScan at 0.3 and no ICEscan at all. That arm exists on disk
(`icescan_cov03/`) and is deliberately **not** used. Whether to lower the
threshold is a separate question, separately measurable, and is not decided here.

> ⚠ **`run_icescan.sh` does not reproduce the scored runs.** It hard-codes
> `--coverage-profile 0.3` at lines 13 (comment) and 34 (command). The runs that
> were scored are at 0.5. The script caches on `best_solution.tsv`, so the next
> `rm -rf work/*/icescan` would silently regenerate the whole benchmark at 0.3.
> **Fix the script before anyone re-runs it.**

---

## 7. AICE — a third class, and an honest statement of its validation

An **AICE** (actinomycete integrative and conjugative element) is found in
*Streptomyces* and relatives. It integrates into the chromosome like an ICE, and
it **does conjugate into another cell** — but it has **no relaxase and no
mating-pair apparatus**. A **TraB / FtsK-SpoIIIE translocase** carries it across
as **double-stranded** DNA, where classical conjugation transfers a single strand
nicked by a relaxase. Possoz *et al.* 2001 (PMID 11679075) concluded this for
pSAM2 from differential SalI methylation — hedged in their own words as
"probably transferred to the recipient as double-stranded DNA", but called "the
first experimental evidence for the transfer of double-stranded DNA during
bacterial conjugation". Vogelmann *et al.* 2011 (PMID 21505418) states it
without the hedge and measured the TraB pore at ~3.1 nm, wide enough for one.
(The translocase is TraB generically; pSAM2's copy is TraSA.)

> **An AICE does make single-stranded DNA, but not the kind that travels.** It
> copies itself by rolling-circle replication *inside* one cell, and that
> intermediate is single-stranded (te Poele *et al.* 2008, PMID 18523858).

The mobility ladder (spec §2.5) has no tier for this, and the reason is
bookkeeping rather than biology. Tiers 5 and 6 are **defined by machinery**: tier
5 is "a relaxase but no mating apparatus, so a helper must supply one", tier 6 is
"a relaxase plus its own". An AICE carries neither protein, so neither tier is
measuring anything about it — and with this branch's validation at nil (below),
reading its real transfer ability as tier 6 would be the worse of the two errors.
It is therefore reported as its own thing:

| field | value |
|---|---|
| `element_type` | `aice` — never folded into `ice` |
| `mge_class` | its own class |
| `mobility_tier` | `NA`, never blank |
| `mobility_tier_reason` | explicit: *"no tier: tiers 5 and 6 are DEFINED by relaxase + type IV secretion conjugation machinery, which an AICE does not carry…"* |
| `mobility` | *"predicted transferable into another cell as double-stranded DNA by a TraB translocase (FtsK/SpoIIIE family); not on the conjugation mobility ladder, whose tiers describe only relaxase + type IV secretion transfer"* |
| `missing_components` | *"none expected (an AICE has no relaxase and no mating-pair apparatus)"* — so the row does not read as a degraded ICE |

`colocalise.py` carries `"aice"` in both `ELEMENT_TYPE_SYNONYMS` and
`CONTEXT_ONLY_ELEMENT_TYPES`; without those it parsed as `None` and logged
*"unrecognised element_type"*.

### ⚠ Validation status at coverage 0.5: **zero**

Two AICEs are called, both on `NC_013929` (*S. scabiei*): 3,336,516–3,342,982
(6,467 bp, de novo boundary) and 5,262,523–5,273,184 (10,662 bp, unbounded).

**Neither overlaps any curated element,** and neither pilot set contains a
curated AICE. The one good AICE recovery on record — AICEScab56241, start +102 bp,
end +44 bp, recovered fraction 0.99, the best boundary this caller has ever
produced — exists **only at coverage 0.3**, which is not the arm we adopt.

So the AICE branch ships, and its validation on this benchmark is nil. Treat
AICE calls as hypotheses.

---

## 8. Full before/after results

Both arms differ only by `--icescan-tsv`. "Ratio" is called length ÷ published
length; "honest" means 0.4 ≤ ratio ≤ 2.0.

### 8.1 ICE pilot — 15/18 → 15/18, honest 9 → 9

| element | true bp | control bp | ratio | union bp | ratio | boundary | conf |
|---|---:|---:|---:|---:|---:|---|---|
| SXT(HN1) | 19,086 | — | — | — | — | | **missed** |
| R391 | 88,532 | 82,887 | 0.94 | 82,887 | 0.94 | none | high |
| ICEKpnQD23-1 | 190,855 | 77,107 | 0.40 | 77,107 | 0.40 | denovo | high |
| ICEEc2 | 93,895 | 55,055 | 0.59 | 55,055 | 0.59 | denovo | high |
| ICE_EfmISMMSVRE1_Tn916 | 16,894 | 45,485 | **2.69** | 45,485 | **2.69** | none | high |
| ICE(Tn4371)6061 | 43,841 | 43,585 | 0.99 | 43,585 | 0.99 | none | low |
| **CMGE(TZ080501)** | 125,779 | 62,749 | 0.50 | **50,425** | **0.40** | none | medium |
| ICEKpn16_GR_13-1 | 151,942 | 148,839 | 0.98 | 148,839 | 0.98 | **tRNA** | high |
| ICEMsp.M1D | 199,376 | 39,437 | 0.20 | 39,437 | 0.20 | none | high |
| ICEPsy10 | 161,009 | — | — | — | — | | **missed** |
| ICEPvuChnBC22 | 148,751 | 51,523 | 0.35 | 51,523 | 0.35 | none | high |
| ICEVflInd1 | 114,195 | 32,992 | 0.29 | 32,992 | 0.29 | none | low |
| SPI-7 | 133,500 | 83,403 | 0.62 | 83,403 | 0.62 | denovo | high |
| ICE_FprA2-165_rpsI | 82,411 | 39,187 | 0.48 | 39,187 | 0.48 | none | medium |
| ICE_EfalV583_tRNALys | 138,318 | 26,040 | 0.19 | 26,040 | 0.19 | denovo | medium |
| ICEB2 | 37,408 | — | — | — | — | | **missed** |
| ICE_Step12228_tRNAser | 38,415 | 24,222 | 0.63 | 24,222 | 0.63 | denovo | medium |
| TR2 | 154,097 | 8,458 | 0.05 | 8,458 | 0.05 | none | medium |

Exactly one scored element changes, and it changes **for the worse** — see §9.

### 8.2 IME pilot — 3/12 → 5/12, honest 1 → 3

| element | true bp | control bp | ratio | union bp | ratio | conf |
|---|---:|---:|---:|---:|---:|---|
| IME_Sag2603_rpmG | 522 | — | — | — | — | missed |
| MGIAmaMed1 | 1,248 | — | — | — | — | missed |
| IME_Ssal57I_tRNAlys | 5,123 | — | — | — | — | missed |
| **IME_CdiR20291_ND** | 5,551 | — | — | **3,840** | **0.69** | medium |
| **Tn4451** | 6,338 | — | — | **5,942** | **0.94** | low |
| MTnPi10 | 7,591 | 6,252 | 0.82 | 6,252 | 0.82 | medium |
| IME_EfalV583_rpsI | 9,405 | — | — | — | — | missed |
| IME_SsuNSUI002_NS | 11,112 | 181,279 | **16.31** | 179,889 | **16.19** | medium |
| MGIVvuTai1 | 19,039 | — | — | — | — | missed |
| IME_FprA2-165_tRNAlys_2 | 23,194 | — | — | — | — | missed |
| SGI1 | 42,451 | — | — | — | — | missed |
| IMEMlNZP2037-1 | 110,480 | 8,164 | 0.07 | 8,164 | 0.07 | medium |

Both new detections are **honest** rather than oversized intervals that merely
contain the curated element: Tn*4451* at ratio 0.94 (start offset +75 bp) and
`IME_CdiR20291_ND` at 0.69 (start offset +37 bp).
These are the two best-bounded IME recoveries the caller has produced.

### 8.3 Whole-set census (all 63 calls, not just scored ones)

**Measured at `1651e6f`** — both arms run over the same 28 genomes.

| | control | union |
|---|---:|---:|
| total elements | 51 | 63 |
| `ice` | 36 | **36 — unchanged** |
| `ime` | 10 | 21 |
| `aice` | 0 | 2 |
| `conjugative_region` | 3 | 2 |
| `genomic_island` | 2 | 2 |
| **predicted self-transmissible** | **36** | **36 — unchanged** |
| confidence `high` | **16** | **16 — unchanged** |
| bounded by a tRNA-anchored att pair | 6 (12%) | 8 (13%) |
| … of `ice`-class calls only | 5/36 (14%) | 5/36 (14%) — identical |

**No element gained a self-transmissible claim.** The net-new calls in the
union arm are all IME or AICE — tier 5 or no tier.

"Bounded by an att pair" now means `boundary_method=tRNA` specifically, because
that is the only method the caller ever acts on — a `denovo` pair is reported
but never applied to the coordinates (see the tuning guide's note on this).
Restricting to `ice`-class calls, where a real element and a fragment are not
conflated, the rate is identical between the two arms: adding ICEscan neither
helps nor hurts how often an ICE gets its true edges. See the box at the top
of this document for why the absolute rate looks lower than the July-30
figure.

---

## 9. Costs — what got worse, re-checked at `1651e6f`

All four were re-run against current code: one no longer happens, the other
three still do, exactly as first measured.

**1. `CMGE(TZ080501)` loses 12,324 bp. Still true.** On `KX077897`:
`48,091–110,839` → `60,415–110,839`; recovered fraction 0.4988 → **0.4009**.
With ICEscan on, an ICEscan `Recombinase` sits **0 bp** from the cluster and
wins the "then the closest" tie-break over a Bakta product-text `Site-specific
recombinase` **11,032 bp** away — and the *more distant* one was nearer the
curated left edge. This lands at 0.4009: roughly 1.1 kb from dropping out of
the "honest" band.

**2. One of two boundary losses still holds; the other has since resolved
itself.**

| | control | union |
|---|---|---|
| `AE009948` | 929,751–986,650 · 56,900 bp · **denovo** · medium | 923,639–941,160 · 17,522 bp · **none** · medium |
| `CP048437_1` | 150,670–187,352 · 36,683 bp · **denovo** · medium | 158,404–187,352 · 28,949 bp · **denovo** · medium |

`AE009948` is unchanged from the original measurement: the union arm still
loses this de novo boundary. `CP048437_1` is not — re-measured, **both
arms now report `denovo`**, so this element no longer loses anything. Nothing
here was deliberately fixed for this locus; it moved as a side effect of the
att-search rework and the loner fix. `high` + `boundary_method=none` is still
10 in both arms.

**3. The nested double-report on `CP011419.1` — RESOLVED, no longer happens.**
In the July-30 run, the union arm emitted *both*:

```
98,234–278,122   179,889 bp  ime   <- the blob, ratio 16.19
246,807–251,765    4,959 bp  ime   <- the honest call, ratio 0.45, start +63 bp
```

Re-run, the union arm gives three clean, non-overlapping calls on this
genome — the honest 4,959 bp element, a separate compact 2,407 bp `ime`
elsewhere on the same contig, and an unrelated `ice` — and **no 179,889 bp
blob at all**. The cause, diagnosed at the time, was a MacSyFinder *loner* gene
(one admitted without the normal co-localisation test) being treated as
grounds to merge two distant clusters into one. That merge rule was corrected
by a later commit specifically because of this locus; its own commit message
names `CP011419.1` as the case that exposed the defect.

**4. Unscored ICE intervals shrink by 10–20 kb. Still true.** Several loci not
covered by either pilot move as ICEscan supplies a closer integrase — e.g.
`NC_004668_1` 89,312 → 69,537 bp (2,177,344–2,266,655 → 2,197,119–2,266,655)
at `high` confidence, `CP048437_1` 69,163 → 59,955 bp
(1,318,921–1,388,083 → 1,328,129–1,388,083). These are unvalidated collateral
of the same tie-break rule.

**Open design question, not resolved here:** when neither candidate integrase
yields an att pair, is "then the closest" the right second key? It is what cost
`CMGE` its 12 kb and what cost the `AE009948` boundary above.

---

## 10. Limitations — read before trusting any output

### There is no false-positive measurement, and this benchmark cannot produce one

`score.py` computes **recall only**. It contains no precision, false-positive or
specificity code at all, and its own docstring says so. This is not an oversight
that can be patched: **every genome in both pilots was chosen *because* it
contains a curated element**, so a call that overlaps nothing curated may be a
genuine second element rather than an error. Nothing here can tell the difference.

Counting them anyway, as an upper bound on the error rate and not an error rate.
**Recomputed at `1651e6f`** against every one of ICEberg's curated entries on
these 28 accessions (1,677 rows, not just the pilots' 30), since that is the
correct denominator for "does this call correspond to something curated":

| | calls | overlap no curated element | unvalidated self-transmissible | unvalidated `high` |
|---|---:|---:|---:|---:|
| control | 51 | 17 (33%) | 9 | 4 |
| union | 63 | 21 (33%) | 9 | 4 |

This is a substantially better picture than the July-30 figures (64%/63%
unvalidated, 23 unvalidated self-transmissible, 9 unvalidated high in both
arms) — corroboration against the full curation roughly doubled. That
improvement tracks the same code changes noted throughout this document (the
att-search rework, the confidence-cap fix, the loner fix). **In both arms**,
about a third of calls still land on nothing curated, and 9 self-transmissible / 4 high-confidence
claims still rest on nothing measurable — the union does not worsen that ratio,
it adds unvalidated IME and AICE calls at medium/low confidence on top of it.

### The negative control — 2 calls in 32.6 Mb, neither invented

Spec §8 Phase 7 asks for "genomes with no reported ICE". Twelve closed genomes
were screened, every one verified absent from all 1,677 ICEberg entries:

| Genome | Why it is here | Calls |
|---|---|---|
| *Buchnera aphidicola* Sg / APS, *Wigglesworthia glossinidia* | reduced endosymbionts, ~0.65 Mb, no MGE traffic | **0** |
| *Prochlorococcus marinus*, *Synechococcus* sp. CC9605 | streamlined marine genomes | **0** |
| *Aquifex aeolicus* VF5 | deep-branching thermophile | **0** |
| *Staphylococcus aureus* N315 | clinical isolate, prophages + transposons | **0** |
| *Escherichia coli* K-12 MG1655 | the most-studied bacterial genome there is | **0** |
| *Pseudomonas aeruginosa* PAO1 | 6.3 Mb reference strain | **0** |
| *Halobacterium* sp. NRC-1 | an ARCHAEON — the models are bacterial, so any call would be spurious by construction | **0** |
| *Bacillus subtilis* 168 | reference strain | **1** `ice` |
| *Salmonella* Typhimurium LT2 | reference strain | **1** `ime` |

**Total: 2 calls across 32,601,100 bp. Neither survives as an error.**

*B. subtilis* 168 → **ICE*Bs1*, a real element**, and the best-evidenced call in
this entire validation: 529,362–549,932 (20,571 bp, ICE*Bs1* is ~20.5 kb), all
four anchor classes, MOBT relaxase (NicK), boundary resolved by a tRNA-anchored
att at **tRNA-Leu(gag)** — ICE*Bs1* integrates at *trnS-leu2* — and Bakta
annotated the integrase literally as "ICE*Bs1* integrase". ICEberg simply has no
entry for NC_000964.3.

⚠ **A warning about this benchmark's own design.** The set originally justified
including strain 168 with "ICE*Bs1* is absent from 168 itself". That is false — ICE*Bs1* was *discovered*
in strain 168 (Auchtung *et al.* 2005). Had the call been counted rather than
inspected, a textbook-correct detection would have been recorded as a false
positive, and a working detector might then have been "fixed". Every call in a
negative control must be examined, never merely counted.

*Salmonella* Typhimurium LT2 → 2,900,443–2,906,075 (5,633 bp), an integrase plus
a protein Bakta labels only "DNA-binding protein" but which CONJscan's `T4SS_MOBM`
profile hits at **i-eval 1.8e-84, coverage 0.986** — a near-full-length relaxase
match. A strong relaxase beside an integrase is what an IME *is*. Whether this
locus is a *bona fide* IME or two co-located genes is a question for a
specialist; it is listed among the open domain questions rather than scored
either way.

**What this does and does not establish.** It does show the caller is not
promiscuous: nine full genomes, including three large well-annotated ones dense
with recombinases, prophages and IS elements, produced nothing at all. It does
**not** give a false-positive rate — "absent from ICEberg" is not "contains no
element", since ICEberg's curation is partial. The honest reading is an upper
bound of **2 candidate false positives in 12 genomes, of which 0 are confirmed
errors**.

### Boundaries: measured before `4a93d89`, and again after

> ⚠ **Superseded numbers.** The two-arm comparison in this document was run on
> 2026-07-30, *before* commit `4a93d89` reworked the *att* search. The arm comparison itself is unaffected — it asks what ICEscan adds, and
> both arms were measured with the same att search — but the absolute boundary
> numbers it produced are no longer current. Both are given below, labelled.

**Before `4a93d89`.** About **one third** of detected elements got a direct
repeat (6 of 18 in the CONJScan-only arm, 6 of 20 in the union arm). On the nine
ICE pilot elements that are both detected and *chromosomal* (the only ones where
boundary error is a meaningful question — see §10, "standalone" deposits), the
median absolute offsets were **14,294 bp at the start and 34,239 bp at the end**,
**identical in both arms** — which was this section's actual finding: ICEscan
moved boundaries not at all.

**After `4a93d89`.** Re-measured on the union arm, which is what ships
(`results_final.tsv` / `ime_results_final.tsv`): **10 of the 20 detected curated
elements (50%)** now carry an *att* pair, and across all calls the rate is
**34 of 63 (54%)** on this 28-genome benchmark. The chromosomal median absolute
offsets are **5,141 bp at the start and 14,912 bp at the end** — roughly a third
and a half of the previous values.

Two cautions on those newer figures. First, mind the denominator: 50% is over
*detected curated elements*, 54% is over *all calls*, and the two are different
questions. Second, the CONJScan-only arm has **not** been re-run since
`4a93d89`, so the "identical in both arms" result above is established only for
the older att search; it is likely but not measured that it still holds.

Either way the shape of the problem is unchanged: half
of all elements are still reported unbounded, their start and end being the
outermost machinery genes rather than the true edges, and the residual error on
the ones that are bounded is still measured in kilobases. Read `boundary_method`
before quoting any coordinate.

### One benchmark entry is broken

`SXT(HN1)` / `AB450045` is deposited as a **19,086 bp** record containing **17
proteins, no relaxase and no T4SS**, standing in for an element that is really
~100 kb. Verified directly: both CONJScan and ICEscan return **zero** hits on it.
No machinery-based caller can find machinery that is not in the record. It is
counted as a miss throughout this document for consistency, but the ICE pilot's
**real ceiling is 17/18, not 18/18**.

### Other standing limitations

- ICEscan and CONJScan are **forks of one another**, so agreement between them
  is weaker evidence than agreement between independent tools would be (§1).
- The AICE branch is **entirely unvalidated at the adopted threshold** (§7).
- Two ICE pilot entries (`ICEPsy10`, `ICEB2`) are missed in both arms;
  seven of the eighteen ICE entries are "standalone" deposits where the element
  *is* the record, so there is no flanking sequence for an att site to sit in and
  the boundary question is trivial.
- Agreement with ICEberg's coordinates is not the same as truth — several
  ICEberg entries are themselves predictions, or simply "the whole deposited
  record".

---

## 11. Licensing

The ICEscan models are **CC BY-NC-SA 4.0** (Institut Pasteur / CNRS) — academic,
non-commercial use only. This is the same licence as the CONJScan models the
module already fetches.

**BacFlux ships neither.** Both are fetched at run time from a URL in
`config/config.yaml`, default `run: false`, exactly as the spec's hard rule
(§11) requires for every licence-encumbered database. Setting `run: true` means
*you* download them under your own agreement with the licensor. Commercial users
should leave this off or obtain permission from the model authors.

No ICEscan model, profile or code is present in this repository, and none of
EBI's CC BY-NC-SA-derived scripts has been copied or adapted (see CITATIONS.md).

---

## 12. Citations

We found this approach by reading other people's work. The credit is theirs.

**ICEfinder2 / ICEscan model set**
> Wang M, Goh Y-X, Tai C, Wang H, Deng Z, Ou H-Y. (2024)
> *ICEberg 3.0: functional categorization and analysis of the integrative and
> conjugative elements in bacteria.* Nucleic Acids Research 52(D1):D732–D737.
> <https://doi.org/10.1093/nar/gkad935>

**CONJScan / MacSyFinder** — the models both sets descend from
> Coluzzi C, Garcillán-Barcia MP, de la Cruz F, Rocha EPC. (2022)
> *Evolution of plasmid mobility: origin and fate of conjugative and
> non-conjugative plasmids.* Molecular Biology and Evolution 39(6):msac115.

> Cury J, Touchon M, Rocha EPC. (2017) *Integrative and conjugative elements and
> their hosts: composition, distribution and organization.*
> Nucleic Acids Research 45(15):8943–8956. <https://doi.org/10.1093/nar/gkx607>

> Abby SS, Cury J, Guglielmini J, Néron B, Touchon M, Rocha EPC. (2016)
> *Identification of protein secretion systems in bacterial genomes.*
> Scientific Reports 6:23080. <http://dx.doi.org/10.1038/srep23080>

> Néron B, Denise R, Coluzzi C, Touchon M, Rocha EPC, Abby SS. (2023)
> *MacSyFinder v2: Improved modelling and search engine to identify molecular
> systems in genomes.* Peer Community Journal 3:e28.

**EBI Mobilome Annotation Pipeline** — the source of the approach
> EBI-Metagenomics `mobilome-annotation-pipeline` (Apache-2.0), which carries its
> own ICEfinder2 attribution. <https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline>

Reading their pipeline is how we learned that ICEscan existed, that it ships the
IME and AICE models, and what coverage threshold it is run at. **No code was
copied** — three of their scripts are CC BY-NC-SA-derived and copying them would
contaminate BacFlux's MIT licence (spec §11). Adopting facts, parameters and
design decisions is permitted and is exactly what was done.

**Element ground truth**
> ICEberg 3.0 curated ICE/IME coordinates, as above.

---

## 13. Bottom line

Adopt the union. ICE detection (15/18), predicted self-transmissible counts (36)
and high-confidence counts (17) are untouched; IME detection goes 3/12 → 5/12
with two genuinely well-bounded new recoveries; all four integrase exclusions
demonstrably fire; the `CP042858.1` regression check holds byte-identically in
both arms.

`run_icescan.sh` needs its threshold fixed; the FIX-4 tie-break's second key
is worth revisiting; and this benchmark still cannot measure a false-positive
rate.
