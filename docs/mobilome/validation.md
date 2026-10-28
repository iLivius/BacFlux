# Validation

This page describes three evaluations of the ICE/IME caller: recovery of curated
elements, comparison with the EBI Mobilome Annotation Pipeline, and sensitivity to
simulated assembly fragmentation. The results describe the tested datasets and software
versions. Read them alongside the [limitations](#what-the-benchmark-does-not-show),
especially when interpreting AMR mobility tiers.

Terms used throughout. An **ICE** (integrative and conjugative element) integrates into
the chromosome and carries its own conjugation machinery, so it can move itself into
another cell. An **IME** (integrative mobilisable element) integrates and carries a
relaxase — the enzyme that nicks the DNA to start transfer — but no complete conjugation
apparatus of its own, so it needs a helper element. BacFlux calls an element `ime` only
when it finds an integrase and a relaxase without a typed mating-pair system; elements
mobilised without an encoded relaxase are outside what it detects
([Tuning](tuning.md#what-no-knob-fixes)). An ***att* site** is the short direct repeat
left at each end of an element when it integrates; a supported pair can help locate an
element's boundaries. **AMR** is antimicrobial resistance, an **MGE** is a mobile
genetic element of any class, and an **AICE** is an actinomycete integrative and
conjugative element.

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
| IME pilot | 12 | The same for IMEs, 10 genera, sizes deliberately straddling the 8,000 bp ICE floor (522 bp to 110 kb) |
| Negative control | 12 | Closed genomes with **no** curated ICEberg entry, 32,601,100 bp in total |

Every accession, with the reason it was chosen, is listed under
[the benchmark genomes](#the-benchmark-genomes) below. The negative set spans reduced
endosymbionts, streamlined marine genomes, a deep-branching thermophile, five full-size
reference strains and one archaeon — the archaeon being a construction check, since the
machinery models are bacterial and a call there would require separate evaluation
outside the models' intended scope.

**How a call is scored.** ICEberg's curated interval is used as the reference. Among a
caller's calls on that genome the scorer keeps the one with the largest overlap, and
reports the start offset, the end offset, and the **recovered fraction** — overlapping
base pairs divided by the curated length. Recovered fraction is a recall measure and
nothing else; the scorer contains no precision or false-positive code, which is why
the negative control and the call-burden count below have to be read alongside it.

!!! warning "The truth standard is a curation, not the truth"

    ICEberg's coordinates are neither perfect nor exhaustive, and this set contains
    one example of each failure: a curated entry whose deposited record cannot
    contain the element it stands for, and a negative-control genome that carries a
    real, well-known element ICEberg has no entry for. The first is counted as a miss
    for every caller; the second scores a correct call as a false positive. Both are
    described below.

## The benchmark genomes

Every genome is a public GenBank deposit, used as downloaded. The sets were chosen
deliberately rather than sampled — the IME set includes small elements and architectures
that challenge the caller. Its recovery rate should not be extrapolated to an unselected
population of IMEs.

### ICE pilot — 18 genomes

Six are classic positive controls; the rest widen the set to 14 genera.

| Accession | Organism | Curated element | Length | Why this one |
|---|---|---|--:|---|
| `AB450045` | *Vibrio cholerae* | SXT(HN1) | 19.1 kb | classic positive control: SXT |
| `AY090559` | *Providencia rettgeri* | R391 | 88.5 kb | classic positive control: R391 |
| `CP042858.1` | *Klebsiella pneumoniae* | ICEKpnQD23-1 | 190.9 kb | classic positive control: ICEKp |
| `GU725392` | *Escherichia coli* | ICEEc2 | 93.9 kb | classic positive control: ICEEc2 |
| `CP012430.1` | *Enterococcus faecium* | ICE_EfmISMMSVRE1_Tn916 | 16.9 kb | classic positive control: Tn916 |
| `KP299160` | *Pseudomonas aeruginosa* | ICE(Tn4371)6061 | 43.8 kb | classic positive control: Tn4371 |
| `KX077897` | *Streptococcus suis* | CMGE(TZ080501) | 125.8 kb | genus diversity: Streptococcus |
| `CP027036.1` | *Klebsiella pneumoniae* | ICEKpn16_GR_13-1 | 151.9 kb | genus diversity: Klebsiella |
| `CP034444.1` | *Mesorhizobium* sp. | ICEMsp.M1D | 199.4 kb | genus diversity: Mesorhizobium |
| `CP024712.1` | *Pseudomonas syringae* | ICEPsy10 | 161.0 kb | genus diversity: Pseudomonas |
| `MH160822.1` | *Proteus vulgaris* | ICEPvuChnBC22 | 148.8 kb | genus diversity: Proteus |
| `GQ463144` | *Vibrio fluvialis* | ICEVflInd1 | 114.2 kb | genus diversity: Vibrio |
| `AL513382` | *Salmonella enterica* | SPI-7 | 133.5 kb | genus diversity: Salmonella |
| `CP048437.1` | *Faecalibacterium duncaniae* | ICE_FprA2-165_rpsI | 82.4 kb | genus diversity: Faecalibacterium |
| `NC_004668.1` | *Enterococcus faecalis* | ICE_EfalV583_tRNALys | 138.3 kb | genus diversity: Enterococcus |
| `CP002188` | *Mycoplasma bovis* | ICEB2 | 37.4 kb | genus diversity: Mycoplasma |
| `NC_004461.1` | *Staphylococcus epidermidis* | ICE_Step12228_tRNAser | 38.4 kb | genus diversity: Staphylococcus |
| `NC_013929` | *Streptomyces scabiei* | TR2 | 154.1 kb | genus diversity: Streptomyces |

### IME pilot — 12 genomes

The selection straddles the 8,000 bp ICE size floor: six of the twelve elements are
smaller. IME-architecture clusters are admitted from a 2,000 bp machinery span, and
element length and machinery-span length are different quantities, so the six are not
all misses: with the ICEscan layer on, three of them are detected.

| Accession | Organism | Curated element | Length | Why this one |
|---|---|---|--:|---|
| `AE009948` | *Streptococcus agalactiae* | IME_Sag2603_rpmG | 522 bp | extreme sub-floor: 522 bp; smallest of the 395 IMEs in ICEberg 3.0 (`ground_truth.tsv`) |
| `CP001103` | *Alteromonas mediterranea* | MGIAmaMed1 | 1,248 bp | deep sub-floor: 1.2 kb MGI, Gram-negative counterpart to the 522 bp case; genus diversity: Alteromonas |
| `CP002888.1` | *Streptococcus salivarius* | IME_Ssal57I_tRNAlys | 5,123 bp | modal sub-floor architecture: S. salivarius tRNALys ~5.1 kb family, 11 near-identical copies in ICEberg |
| `NC_013316.1` | *Clostridioides difficile* | IME_CdiR20291_ND | 5,551 bp | sub-floor, no named integration site (ND) so no tRNA hint for the att search; genus diversity: Clostridioides |
| `U15027` | *Clostridium perfringens* | Tn4451 | 6,338 bp | named IME control: tnpX serine recombinase + oriT, no T4SS; curated IME used to test the [classification rules](output.md#an-ice-cannot-sit-on-a-plasmid); 1.7 kb below the 8,000 bp ICE floor (`--min-element-bp`) |
| `NC_017860` | *Prevotella intermedia* | MTnPi10 | 7,591 bp | floor edge case: 7591 bp, only 409 bp below min_element_bp; sits on a 580 kb second chromosome |
| `NC_004668.1` | *Enterococcus faecalis* | IME_EfalV583_rpsI | 9,405 bp | just above the floor; same chromosome as ICE pilot entry ICE_EfalV583_tRNALys, so it tests ICE vs IME discrimination in one genome |
| `CP011419.1` | *Streptococcus suis* | IME_SsuNSUI002_NS | 11.1 kb | above floor, no named integration site; genus diversity: Streptococcus suis |
| `BA000037` | *Vibrio vulnificus* | MGIVvuTai1 | 19.0 kb | 15-30 kb MGI stratum; genus diversity: Vibrio |
| `CP048437.1` | *Faecalibacterium duncaniae* | IME_FprA2-165_tRNAlys_2 | 23.2 kb | same chromosome as ICE pilot entry ICE_FprA2-165_rpsI, so it tests ICE vs IME discrimination; genus diversity: Faecalibacterium |
| `AF261825` | *Salmonella enterica* | SGI1 | 42.5 kb | named IME control: canonical AMR-carrying SGI1 (ACSSuT / In104), mobilised by an IncA/C helper plasmid; borderline deposit, one-sided flank |
| `CP016079.1` | *Mesorhizobium loti* | IMEMlNZP2037-1 | 110.5 kb | largest stratum: 110 kb element on a 7 Mb chromosome, tests 15 kb anchor clustering across a long element; genus diversity: Mesorhizobium |

### Negative control — 12 genomes, 32.6 Mb

Genomes with **no** curated ICEberg entry. Any call here is a candidate false positive.

| Accession | Organism | Why this one |
|---|---|---|
| `NC_004061.1` | *Buchnera aphidicola* Sg | reduced endosymbiont genome (0.64 Mb): no conjugation machinery expected at all |
| `NC_004344.2` | *Wigglesworthia glossinidia* brevipalpis | reduced endosymbiont (0.70 Mb): obligate intracellular, no MGE traffic |
| `NC_002528.1` | *Buchnera aphidicola* APS | second Buchnera lineage, independent confirmation |
| `NC_005042.1` | *Prochlorococcus marinus* SS120 | streamlined marine cyanobacterium, 1.75 Mb |
| `NC_007516.1` | *Synechococcus* sp. CC9605 | marine picocyanobacterium, no curated ICE |
| `NC_000918.1` | *Aquifex aeolicus* VF5 | deep-branching thermophile, 1.55 Mb |
| `NC_002745.2` | *Staphylococcus aureus* N315 | well-characterised clinical isolate, no ICEberg entry |
| `NC_000913.3` | *Escherichia coli* K-12 MG1655 | the reference lab strain; well-characterised genome with no curated ICEberg entry in this benchmark |
| `NC_002516.2` | *Pseudomonas aeruginosa* PAO1 | reference strain, extensively studied |
| `NC_000964.3` | *Bacillus subtilis* 168 | **Mis-assigned when the set was built.** ICE*Bs1* is present in 168 — it was discovered there (Auchtung 2005) — but ICEberg has no entry for it. The caller finds it, so this is a true positive scored as a false alarm. It is left in place rather than quietly removed |
| `NC_003197.2` | *Salmonella Typhimurium* LT2 | reference strain, no curated ICEberg entry |
| `NC_002607.1` | *Halobacterium* sp. NRC-1 | archaeon: CONJscan/ICEscan models are bacterial, so any call is suspect |

## ICE pilot — 15 of 18

That is 9 of the 11 with real chromosomal context, and 6 of the 7 where the deposited
record *is* the element and there is no flanking sequence for an *att* site to sit in.

One of the three misses cannot be found by any caller: the deposited `SXT(HN1)` record
carries no relaxase and no conjugation apparatus, so 17 of 18 is the ceiling on this
set. Each miss has a stated cause:

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

The median recovered fraction is 0.69. IME recovery was lower than ICE recovery in
these selected pilot sets.

The misses concentrate at the small end — 522 bp, 1,248 bp and 5,123 bp — and at SGI1
(42 kb) and MGI*Vvu*Tai1 (19 kb). Elements under about 2 kb rarely give any relaxase
signal. The 5.1 kb element is lost to integrase assignment instead: a nearer integrase
joins the machinery cluster, and the element's own integrase is left with nothing to
attach to (see [Tuning](tuning.md#an-element-you-know-is-there-was-not-reported)).

The boundary search uses the 8 kb ICE floor even for candidates admitted under the
2 kb IME machinery floor: it rejects every repeat pair implying an element shorter than
8 kb. Three of the five detected curated IMEs are shorter than that (5,551, 6,338 and
7,591 bp), so a pair at their true ends would be rejected. An unresolved call retains
its machinery span and `boundary_method=none`. Neither floor is a configuration key,
and a per-class *att* floor has not been measured. See [Tuning](tuning.md) for the two
separate thresholds.

## Negative control — 2 calls in 32.6 Mb, ICEscan layer off { #negative-control-2-calls-in-326-mb }

Ten of the twelve genomes returned nothing at all, including *E. coli* K-12 MG1655,
*P. aeruginosa* PAO1 and *S. aureus* N315 — three large, well-annotated genomes dense
with recombinases, prophages and insertion sequences. With the default models the
archaeon is among them, which is the construction check passing; with the ICEscan layer
on it is not (below).

With the default machinery models alone, two calls were made. One corresponds to a
known ICE; the other remains unresolved. With the ICEscan layer on, the same genomes
gave five: the same two plus three low-confidence `cime_or_island` calls on the archaeon
*Halobacterium* sp. NRC-1, built from integrase and type IV secretion profile hits with
no assembled system (`evidence_level=profile_hits_only`). Both arms, re-run with the
v2.1.0 code and models, give the same calls at the same coordinates.

**_B. subtilis_ 168 → ICE*Bs1*, a real element ICEberg has not catalogued.** Called at
529,362–549,932, so 20,571 bp against an element of about 20.5 kb; all four anchor
classes present; a MOBT relaxase; and the boundary resolved by a tRNA-anchored *att*
site at a tRNA-Leu, which is where ICE*Bs1* is known to integrate. Bakta independently
annotated the integrase as "ICE*Bs1* integrase". This is the best-evidenced single
call in the whole validation.

**_S._ Typhimurium LT2 → an open question, not an answer.** An integrase sits beside a
protein Bakta labels only "DNA-binding protein", which the `T4SS_MOBM` relaxase profile
hits at i-evalue 1.8 × 10⁻⁸⁴ over 98.6% of the profile — a near-full-length relaxase
match. A strong relaxase match beside an integrase supports an IME candidate, but does not
establish a functional element. Further examination is needed to distinguish an IME
from co-located genes.

!!! warning "These are candidate calls, not a false-positive rate"

    "Absent from ICEberg" is not "contains no element". With the default models the
    result is **two candidates in 12 genomes: one known element and one unresolved
    candidate**.
    This does not establish a zero false-positive rate.

    The set originally justified including *B. subtilis* 168 on the grounds that
    ICE*Bs1* is absent from 168 itself. **That is false — ICE*Bs1* was discovered in
    strain 168** (Auchtung *et al.* 2005). Had the calls been tallied instead of
    inspected, a textbook-correct detection would have been recorded as a false
    positive and a working detector would then have been "fixed" to stop making it.
    Every call in a negative control has to be examined, never merely counted.

## What the second model set adds

The optional ICEscan layer ([Enabling the module](enabling.md)) supplies the integrase and
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

This table is the original run. With the current code the default models find 2 of the
12: the third detection above was a 181,279 bp call laid over an 11,112 bp IME, counted
only because the scorer accepts any overlap. The ICEscan arm still finds 5, so the layer
now adds three IMEs: the two below, and `IME_SsuNSUI002_NS`, now recovered as a
4,959 bp call starting 63 bp from its curated start.

**No element gained a self-transmissible claim.** Every net-new call in the ICEscan
arm is an IME or an AICE — tier 5 or no tier at all — and both arms hold 36 `ice`
calls. ICEscan cannot supply a mating apparatus, because its mating-pair and
coupling-protein hits are never used as anchors. Its integrase or relaxase hits can
complete an ICE that the default models alone left as a conjugative region or island;
on these 28 genomes that did not happen.

The two new detections of the original run fit the curated element closely rather than
being oversized intervals that merely contain it: Tn*4451* at 0.94 of its curated length
(start offset +75 bp) and `IME_CdiR20291_ND` at 0.69 (start offset +37 bp) — the two
best-bounded IME recoveries the caller has produced.

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
run over both. Using one scorer makes the evaluation procedure consistent, but does not remove
its limitations or guarantee that they affect both callers equally. Two differences
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
difference is consistent with implementation choices around that engine, including
candidate seeding, size floors and low-quorum reporting. Gene prediction and other
version differences also remain part of the comparison.

ICE*Vfl*Ind1 is detected but classified `cime_or_island` — integrase, coupling protein
and apparatus, **no relaxase** — not `ice`. It counts as detected because its call
overlaps the curated interval, which it does. It is not a correct class call and should
not be presented as one.

### Class

Of the 17 elements both callers found, the classes agree on 14 and differ on 3. All
three are ICEberg-curated **ICEs that the other pipeline called IMEs**, and BacFlux
agrees with the curation on all three: SPI-7, `ICE_EfalV583_tRNALys` and
`ICE_Step12228_tRNAser`.

ICE and IME labels imply different predicted transfer requirements and can affect
tiers 6 and 5 when an AMR gene is assigned to the interval. They do not by themselves
establish whether the gene was acquired, whether it is expressed, or whether transfer
occurs under the conditions of interest.

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

BacFlux's 65 calls are 63 on the 28 pilot genomes with the ICEscan layer on and 2 on
the 12 negative genomes with CONJScan alone. With the layer on all 40 genomes, as in the
EBI pipeline, BacFlux makes 68 calls over 2,687,359 bp: 6.7% fewer base pairs than the
EBI pipeline while recovering three more curated elements. Neither caller achieves
recall by volume. This table has to be read
next to the recall figures, because the scorer keeps a caller's *largest overlapping*
call: a tool that emitted one enormous call per genome would score well on recovered
fraction and badly on nothing the scorer measures.

### Negative set, and the one call both callers made

With BacFlux's default models, two calls each across the 32.6 Mb, on different genomes
— BacFlux on *S.* Typhimurium LT2, the EBI pipeline on *S. aureus* N315 — plus the
shared call on *B. subtilis* 168. With the ICEscan layer on, BacFlux adds three
low-confidence calls on the archaeon.

On that shared call the two pipelines return **base-pair-identical coordinates, the
same 60 bp flanking repeat, the same MOBT relaxase and the same mating-pair type**, by
two different boundary algorithms. Because boundary refinement is the genuinely
independent layer, this agreement carries information, and what it evidences is that
ICE*Bs1* is a real element missing from the catalogue.

### Boundaries

Measured on the elements both callers found — nine chromosomal, five standalone.
Before the boundary rework, the other pipeline's boundaries were better than ours;
that is what prompted the fixes.

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

**ICE detection was unchanged down to 50 kb N50 in this experiment.** ICE recall is
identical to the closed genomes at both 150 kb and 50 kb, and only breaks at 20 kb.

**The length collapses, but the caller is not what collapses it.** Against the true
element, median recovered length falls from 0.55 to 0.22. Against the largest
surviving contig piece — the ceiling any per-contig caller can reach — it is flat or
better. The caller keeps recovering the same share of what is still visible; the
missing bases are the assembly's.

**No new sequence or high-confidence class promotion was observed in this sweep.**
Across the 191 draft calls, no reported *att* sequence was absent from the closed genome
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
published literature. With all four optional layers on, that run reaches tiers 1, 3, 4,
5 and 6. Tier 2 is missing only because the TnCentral layer names the transposon
Tn*Ecp1.1* around the chromosomal `bla`CTX-M-15. With the layer that gene scores tier 4;
without it, tier 2 — see [below](#tiers-4-and-5-the-thinnest-tested-part-of-the-ladder).
The run recovers the chromosomal ICE at 0.946 of its curated length with a
tRNA-anchored boundary and still ends 3,138 bp inside the curated element, because a
called interval is a floor.
The walkthrough is the [worked example](worked-example.md); the columns are on
[Reading the output](output.md).

## Unit tests

The helper scripts are importable Python and are tested outside Snakemake:

```bash
pytest workflow/scripts
```

Three tests are skipped unless pointed at real GTDB-Tk and GTDB taxonomy files
(`BACFLUX_TEST_GTDBTK_SUMMARY`, `_SAMPLE` and `_GENUS`, and
`BACFLUX_TEST_GTDB_TAXONOMY`); most tests cover the mobilome helpers. These tests check
software behaviour rather than biological validity — that a repeat planted at a known
offset is found, that a confidence cap cannot be bypassed, that a malformed input is
refused rather than silently parsed. Nothing in them says a call is correct.

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

A tier-1 label means that the available evidence did not assign a stronger mobile
context. It does not establish intrinsic resistance, whether or not another element
was detected elsewhere in the genome. Only a gene within 5 kb of an element call is
flagged (`mge_context = near_<class>`); `bla`IMP-8, 24 kb past the called edge, carries
no flag. In a genome with any ICE or IME call whose `boundary_method` is `none`, check
the distance from a tier-1 gene to the nearest call before reading "intrinsic
candidate" at face value.

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
of 63 (33%) with it, overlap no curated element. In each arm, 9 predicted
self-transmissible calls and 4 high-confidence calls overlap no curated element.

This is the **fraction without curated overlap, not a false-positive rate**. Every
genome in both pilots was chosen *because* it contains a curated element, so a call
overlapping nothing curated may be a genuine second element. Nothing in this benchmark
can tell the difference, and its scorer contains no code that could.

One class is unvalidated outright. **AICE** — the actinomycete class, which transfers
double-stranded DNA through a translocase instead of a relaxase — is called twice on
one *Streptomyces* genome, neither call overlaps a curated element, and neither pilot
set contains a curated AICE. Treat AICE calls as hypotheses.

### Tiers 4 and 5, the thinnest-tested part of the ladder

Tier 6 and tier 1 carry the benchmark work above. The middle of the ladder does not.

- **Tier 4** — inside a *named* transposon or integron — depends entirely on the opt-in
  TnCentral naming layer, and has now been measured on two genomes. On the
  *K. pneumoniae* ATCC BAA-2146 positive control, `bla`CTX-M-15 scores tier 4 inside
  Tn*Ecp1.1* at 99.9% identity over 87% of the 3,417 bp reference; the same run without
  the layer scores that gene tier 2, so the override behaves as designed
  ([Worked example](worked-example.md)). On *Enterobacter hormaechei* it produced **15
  rows at once**, all fifteen genes inside a single curated chromosomal transposon,
  Tn*SMR478* ([hybrid mode on five closed genomes](../methods_reference_genome_run.md)).

  Those two results say different things. One gene in one element is the path working;
  fifteen genes in one element is a reminder that a tier 4 count reports **elements as
  much as genes** — a single transposon carrying a resistance cluster produces as many
  rows as it has genes, and reading fifteen rows as fifteen independent findings would
  overstate the evidence considerably.

  Tier 4 stays uncommon for two structural reasons rather than one measurement. It is
  reached only when a curated name is the *strongest* evidence available — six other
  curated elements in the BAA-2146 genome sit on plasmids and score 5 or 6 instead — and
  the ≥80% reference-coverage rule is strict: on fragmented clinical assemblies the same
  Tn*Ecp1.1* was refused at 12% and 49% coverage. Expect tier 4 on closed genomes far
  more than on drafts.
- **Tier 5** is exercised, but only by one of its two routes. Every tier 5 row measured
  so far is a gene on a plasmid (mobilisable, or typed non-mobilisable like pCuAs in the
  worked example). **No AMR gene in any run has been assigned an IME context**, so the
  "inside an IME" half of tier 5 is untested end to end — which matters, because the
  IME pilot is the weaker half of the module to begin with.

Validation coverage differs between routes. Tier 4 has limited end-to-end examples;
the IME route to tier 5 has not been exercised with an assigned AMR gene in these runs.
Report those gaps alongside the boundary-related limitations of tiers 1 and 6.

## On vmatch

The *att*-site search does not use vmatch, the string-matching tool the reference ICE
finders use for this step. BacFlux's own search targets the same class of exact maximal
repeat.

BacFlux's search finds **exact maximal repeats between the two flanking windows in
standard-library Python** — vmatch's own semantics. To test that rather than assert it,
a second implementation was built on vmatch's actual data structure (prefix-doubling
suffix array, Kasai LCP array, cross-flank maximal-match enumeration) and run against
the shipped one on every flank window the caller genuinely searched — after
insertion-sequence masking, at the real element coordinates, at the real computed
repeat-length floor:

> **52 benchmark genomes — a wider set than the 40 scored above — 246 *att* searches,
> 180 distinct flank windows, 1,033 repeats, exact set equality on every window.** No
> window hit either implementation's internal cap.

The two implementations returned the same set of repeats on every window tested.
vmatch itself was not run: the comparison partner reimplements its method, so the
result shows that the shipped search returns the exact maximal repeats it is defined
to return. Agreement between two implementations of one definition says nothing about
whether a boundary is biologically real.

The Bioconda recipe for vmatch 2.3.1 declares `license: Unknown / OTHER`, and no
licence text could be retrieved from the upstream site when this was checked in July
2026, so its terms could not be established. BacFlux does not use vmatch.

!!! note "An earlier figure that should not be quoted"

    A previous version of this measurement was recorded as "35 of 35 windows, 122
    repeats". It does not reproduce under any scoping that could be reconstructed,
    because the script that produced it was never kept. The equality result itself
    reproduces, over roughly five times as many windows, and the numbers above are the
    ones to use. The script that produces them is not distributed.

## No code from the EBI pipeline is used here

Three of that pipeline's scripts are derived from ICEfinder2 and carry CC BY-NC-SA 4.0
headers. They were read to understand an approach, which is all that was done; what
was adopted from them is conventions — Sequence Ontology terms, an element ID format,
the discard-with-reason pattern and two numeric thresholds. Licence types for every
tool, model set and database the module can use are listed on
[Licensing](../about/licensing.md).

## Where the numbers come from

The benchmark design is in the repository, under
[`benchmark/mobilome/`](https://github.com/iLivius/BacFlux/tree/main/benchmark/mobilome):
the three genome sets, with the reason each genome was chosen, and the scoring script,
`score.py`. The genomes and the full run outputs are not in the repository. The genomes
are public accessions, so the closed-genome scoring can be repeated with `score.py`. The
scripts that cut the draft assemblies and converted the EBI output are not distributed.
The method write-ups these figures are drawn from are
[Head-to-head vs the EBI pipeline](../methods_ebi_comparison.md),
[CONJScan + ICEscan together](../methods_icescan_union.md) (the two-arm model-set
experiment), [att sites and small plasmids](../methods_att_and_small_plasmids.md) (the
*att* search and the vmatch measurement), [Draft assemblies](draft-assemblies.md) (the
fragmentation sweep) and the [Worked example](worked-example.md) (the end-to-end run).

Call totals depend on the genome set and on the ICEscan layer. On the 28 pilot
genomes: 51 calls without ICEscan and 63 with it. On all 40 genomes: 65 in the
head-to-head with the EBI pipeline, where the pilot genomes ran with ICEscan and the 12
negative genomes without it (2 calls), and 68 in the fragmentation sweep's closed
baseline, all with ICEscan (5 negative-control calls). Each figure is labelled with
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
