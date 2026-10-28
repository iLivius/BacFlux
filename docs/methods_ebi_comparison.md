# Comparing BacFlux's ICE/IME caller against the EBI Mobilome Annotation Pipeline

This page compares BacFlux's ICE/IME caller with the EBI Mobilome Annotation
Pipeline on the same genomes, using one scoring procedure. It assesses curated
element recovery, class assignment and boundaries within the tested configuration.

The two callers share ancestry, so this is **not** an independent check —
treating it as one would overstate what agreement between them proves.


!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **MAP** | the EBI **M**obilome **A**nnotation **P**ipeline, the caller compared against |
    | **ICE** | integrative and conjugative element — integrates into the chromosome and carries its own conjugation machinery, so it can move itself |
    | **IME** | integrative mobilisable element — integrates and can be moved, but has no machinery of its own and needs a helper |
    | **AICE** | actinomycete integrative and conjugative element, a third class using a different transfer mechanism |
    | **T4SS / T4CP** | type IV secretion system and its coupling protein — the apparatus that pushes DNA into the next cell |
    | ***att* site** | the short repeat marking an integrated element's ends |
    | **SO** | Sequence Ontology, the controlled vocabulary used for feature types in GFF output |
    | **GFF** | the tab-separated annotation format both pipelines write features to |
    | **AMR** | antimicrobial resistance |
    | **CDS** | coding sequence — one predicted protein-coding gene |
    | **IS** | insertion sequence, the smallest kind of mobile element |
    | **MPF** | mating-pair formation, the class of the mating apparatus a conjugative element carries |
    | **MAG** | metagenome-assembled genome, a genome reconstructed from a mixed community |

---

## 0. The two callers are not independent

BacFlux's ICE/IME caller and the EBI Mobilome Annotation Pipeline (MAP) both find
conjugation machinery by running **MacSyFinder over closely related model sets**.
MAP runs the ICEscan model set inside its `ICEFINDER2_LITE` subworkflow. BacFlux
runs CONJScan 2.1.0, of which ICEscan is a fork, and on the 28 pilot genomes of this
comparison the optional ICEscan layer as well, in the `conjscan` / `icescan` rules of
`workflow/rules/shared/80_mobilome.smk`. The MacSyFinder versions differ slightly —
2.1.4 in MAP's container, 2.1.6 in BacFlux's conda environment
(`workflow/envs/macsyfinder.yaml`) — and the profile HMMs and system definitions that
decide *"is there a conjugation system here at all"* overlap heavily without being
identical.

The consequence must not be glossed over:

> **Detection agreement is not independent confirmation.** When both callers
> find the same element, that is substantially the same engine agreeing with
> itself. It is *not* two independent lines of evidence, and it must never be
> reported as independent confirmation of detection.

What *is* genuinely independent, and therefore what this comparison can actually
measure:

| Layer | Shared or independent? |
|---|---|
| Machinery detection (does a conjugation system exist here) | **Largely shared** — overlapping model sets, same engine |
| Clustering (which machinery genes belong to one element) | **Independent** — different code, different rules |
| Boundary refinement (where does the element start and stop) | **Independent** — MAP runs `vmatch` plus an ICEfinder2-derived script; BacFlux runs its own exact maximal-repeat search |
| Class assignment (ICE vs IME vs AICE) | **Partly independent** — same model vocabulary, different decision rules |

So read the detection numbers as a *consistency check* — did two implementations
of the same idea behave the same way, and where did they diverge — and read the
boundary and class numbers as the real comparison.

One place where the shared engine does **not** apply, and where agreement
therefore does carry information, is §6.2: on *Bacillus subtilis* 168 the two
pipelines converge on byte-identical coordinates and the same 60 bp repeat by two
different boundary algorithms. That agreement is about **boundaries**, which are
independent, so it means something.

---

## 1. What each pipeline is for

These are not competing implementations of one task. They answer different
questions, and the overlap is one subworkflow.

**EBI Mobilome Annotation Pipeline (MAP)** — a broad **mobilome census**. Given an
assembly, it inventories every class of mobile element it can: insertion
sequences (ISEScan), integrons (IntegronFinder), plasmids and prophages
(geNomad, CheckV), ICEs and IMEs (ICEfinder2-lite), plus compositional outliers,
and merges them into one reconciled GFF. Its design centre is **metagenomes and
MAGs**, where the question is *what mobile elements are present in this
community* and the input is fragmented, taxonomically mixed and un-curated. It is
Nextflow, container-based, and carries optional functional-annotation branches
(AMRFinderPlus, DeepARG, RGI, antiSMASH, SanntiS, GECCO) for downstream biology.

**BacFlux's mobilome module** — an **AMR-mobility assessment on a single
isolate**. It starts from the AMR genes actually detected in one characterised
strain and asks, per gene, *is this gene inside a mobile element, and how
transferable is that element*, placing each on six mobility tiers from
"chromosomal, intrinsic candidate" to "**predicted** self-transmissible". The
regulatory framing is EFSA's intrinsic-versus-acquired distinction for strains
entering the food and feed chain. The module is **off by default** and is one
stage of a workflow whose main business is assembly, decontamination, taxonomy,
annotation and AMR calling.

**These are complementary scopes.** MAP goes wider across element classes and is
built for community-level data; BacFlux goes narrower but ties every call back to
a specific resistance gene in a specific isolate and reports a mobility tier with
a confidence level. Neither is a reduced version of the other, and this document
should not be read as ranking them. The difference is the unit of analysis, not
the tool list: MAP runs AMRFinderPlus too, on by default. What differs is that it
is built for metagenomes and MAGs, where a contig's host is itself uncertain,
while BacFlux assumes one isolate and can therefore attach a mobility tier to a
named gene on a named replicon.

A practical consequence for users: for a full
mobilome inventory of an isolate, or for anything metagenomic, MAP is the right
tool and BacFlux writes `{sample}_contigs.fna` and `{sample}.gbk` ready to feed it.
Using MAP's `mobilome.gff.gz` output does not incorporate MAP code into BacFlux;
output and source-data terms should be considered separately (§10).

---

## 2. How the comparison was run

### 2.1 The genome set

Forty closed reference genomes, in three groups, all downloaded from GenBank and
used as-is:

| Set | n | What it is |
|---|---|---|
| ICE pilot | 18 | Genomes/deposits carrying an ICE with **curated coordinates in ICEberg**, spanning 14 genera |
| IME pilot | 12 | Same, for IMEs, 10 genera, sizes deliberately straddling the module's size floor (522 bp to 110 kb) |
| Negative control | 12 | Closed genomes with **no** curated ICEberg entry: 32,601,100 bp total |

Two accessions (`CP048437.1`, `NC_004668.1`) carry both a curated ICE and a
curated IME, so 40 genomes yield **30 curated elements**. The negative set is
three reduced endosymbionts, two streamlined marine genomes, a deep-branching
thermophile, five full-size reference strains (*E. coli* K-12 MG1655, *B.
subtilis* 168, *P. aeruginosa* PAO1, *S. aureus* N315, *S.* Typhimurium LT2) and
an archaeon (*Halobacterium* sp. NRC-1). The archaeon is a **construction
check**: CONJscan and ICEscan are bacterial models, so any call there would be
outside the models' intended scope and would require separate evaluation.

Both callers were given **the same FASTA files**: the benchmark genomes, fetched by
accession from the set files in `benchmark/mobilome/`. No re-assembly, no
re-annotation on our side of the fence — MAP does its own gene calling (Prodigal),
BacFlux uses Bakta, which is one of the differences the comparison contains rather
than controls for.

### 2.2 Scoring

Both callers' output was reduced to one table per genome,
`work/{sample}/ice_elements.tsv`, with the columns `benchmark/mobilome/score.py`
reads. `score.py` was then run **unmodified** over both. It:

- uses the curated ICEberg interval as the reference;
- among a caller's calls on that genome, selects the one with the **largest
  overlap** with the curated interval;
- reports start offset, end offset, and `recovered_fraction` = overlapping bp ÷
  curated length.

Using one scorer makes the procedure consistent. Its limitations (§7) still apply
and may affect callers differently depending on their output.

### 2.3 Which MAP, and the two environment fixes needed to run it

MAP **v5.0.0** (`ebi-metagenomics/mobilome-annotation-pipeline`, repo checkout at
`3aa408d`), run with `-profile singularity`. Recorded tool versions inside the
containers: MacSyFinder 2.1.4, ISEScan 1.7.3, IntegronFinder 2.0.6, geNomad
1.11.1, HMMER 3.4, Aragorn 1.2.41.

Two environment problems can make MAP fail. The fixes for both are in
`benchmark/ebi_comparison/nf_env.sh` and `benchmark/ebi_comparison/run_map.sh`.
Without them, a fresh setup may hit either.

**(a) Apptainer instead of an old Singularity.** Singularity **2.x** predates the
SIF container format that all current biocontainers use, so it cannot run them at
all. If your host's Singularity is that old, install Apptainer (the renamed, modern
continuation of Singularity) into a conda environment of its own and put it ahead of
the system binary on `PATH`. Running it without root needs a kernel with
unprivileged user namespaces enabled. The conda package ships a `singularity` →
`apptainer` symlink, which lets Nextflow's stock `singularity` profile work
unmodified.

**(b) `NXF_VER=24.10.6`.** On the tested Nextflow 26.04.6 build, MAP v5.0.0 dies inside
its own ICEfinder2-lite subworkflow:

```
ERROR ~ Invalid method invocation `call` with arguments:
        [[id:pos_test], null, .../pos_test_uniprot_names.tsv] on _closure11 type
```

This is the `.join(..., remainder: true)` feeding `REFINE_BOUNDARIES` in
`subworkflows/local/icefinder2lite.nf` handing a three-element tuple to a closure
written for six. Newer Nextflow stopped tolerating that arity mismatch. MAP's
manifest only requires `!>=24.04.0` and its CI does not pin a version, so this is
a genuine forward-compatibility gap in MAP rather than a local misconfiguration.
Pinned to 24.10.6, the same commit runs `REFINE_BOUNDARIES` to completion.

**(c) Functional branches switched off.**
AMRFinderPlus, DeepARG, RGI-CARD, PathoFact2, SanntiS, GECCO and antiSMASH were
disabled: they answer a different question, each needs its own multi-GB database,
and they contribute nothing to a boundary comparison. geNomad and CheckV are
**not** skippable in MAP and feed its overlap integrator, so they ran. Note these
have to be set in a config file, not as `--skip_x` on the command line: MAP
validates parameters with nf-schema in strict mode and a bare CLI `--skip_x`
arrives as the string `"true"`, failing the boolean type check.

### 2.4 Completeness of the MAP run

**40/40 genomes reached `COMPLETED`** (`benchmark/ebi_comparison/map_status_40.tsv`). Median
wall-clock 30.3 min per genome (range 23.4–39.9). Twelve genomes produced no
`ices.tsv` because no ICE candidate survived MAP's prescan; that is a real
result — no call — not a failure, and it is scored as zero calls.

Three run-level caveats:

- `CP016079.1` — a **non-ICE** stage aborted (`GFF_MAPPING:ABORTED`). The
  ICEfinder2-lite table was written and is what we scored, so the ICE/IME result
  stands, but MAP's merged GFF for that genome is incomplete.
- `NC_000964.3` and `U15027` — some stages were reused from an earlier run via
  Nextflow `-resume`. `U15027` shows 0.1 min elapsed because it was served almost
  entirely from cache. The status table records which run directory is canonical
  for each accession, and the converter reads that table rather than guessing.

---

## 3. Converting MAP's output so the same scorer could read it

MAP writes a merged GFF that pools every predictor it runs. Turning that into the
two-column-plus-metadata table `score.py` expects required three judgement calls.
Each is described below. The conversion script, which is not in the repository,
does **no** filtering, no re-thresholding and no coordinate arithmetic — only
translation.

### 3.1 The SO-term `integron` trap

MAP labels rows with Sequence Ontology terms. It writes the term **`integron`
(SO:0000365) for two completely different things**: an ICEfinder2-lite **IME**,
and a genuine **IntegronFinder integron**. Filtering MAP's GFF on the `type`
column alone therefore imports integrons into MAP's ICE/IME score.

Counted directly across all 40 merged GFFs:

| Source | GFF `type` | Rows |
|---|---|---|
| `ICEfinder` | `conjugative_integron` | 37 |
| `ICEfinder` | `integron` | **27** ← these are IMEs |
| `IntegronFinder` | `integron` | **6** ← these are integrons |
| `ICEfinder` | `direct_repeat_element` | 22 |
| `ISEScan` | `insertion_sequence` | 926 |
| `IntegronFinder` | `attC_site` | 279 |
| `geNomad` | `prophage` / `plasmid` | 67 / 7 |
| `MAP` | `compositional_outlier` | 657 |

A naive `type == "integron"` filter would have pulled **6 unrelated
IntegronFinder integrons** into MAP's element set — inflating its call count and
its called base pairs, and scoring it against ICE coordinates it never claimed.
Nothing would have errored.

The fix: **filter on the GFF `source` column first** (`ICEfinder` only), then read
the class from the `mobile_element_type` **attribute**, never from the type
column. `direct_repeat_element` rows are excluded too — those are the att sites
bounding an element, not elements.

### 3.2 Class vocabulary

MAP's `mobile_element_type` maps one-to-one onto ours for the three classes it
can emit:

| MAP | BacFlux |
|---|---|
| `T4SS-type_ICE` | `ice` |
| `IME` | `ime` |
| `AICE` | `aice` |

The `_with_DRs` suffix is **not** a fourth class — it means ICEfinder2-lite
refined that element's boundary onto a flanking direct-repeat pair. It is
stripped from the class and recorded as `boundary_method=direct_repeat`, which is
where BacFlux keeps the same fact. We deliberately do **not** claim
`boundary_method=tRNA` for MAP: its `close_to_RNA` column is a separate
annotation and does not mean the boundary was derived from the tRNA.

BacFlux has two classes MAP has no equivalent for (`cime_or_island`,
`conjugative_region`); MAP has none we lack. Noted, not scored.

### 3.3 Conversion cross-check

The converter reads MAP's merged GFF, then cross-checks against ICEfinder2-lite's
own `{sample}_ices.tsv` — the same elements in a different layout. **All 40
genomes agree**: 64 elements in the GFF, 64 rows in the TSV, decomposing as 32
ICEs + 27 IMEs + 5 AICEs. All five AICE calls fall in *Streptomyces scabiei*,
which is consistent with the class's taxonomic scope but does not validate
the individual calls.

---

## 4. Detection

Thirty curated elements, both callers, same genomes.

| | Count |
|---|---|
| Found by **both** | 17 |
| Found by **BacFlux only** | **3** |
| Found by **MAP only** | **0** |
| Found by neither | 10 |
| **Totals** | **BacFlux 20/30, MAP 17/30** |

The three BacFlux-only detections are **ICEVflInd1** (*Vibrio fluvialis*),
**Tn4451** (*Clostridium perfringens*) and **MTnPi10** (*Prevotella intermedia*).
In each, MAP's prescan produced no ICE candidate and no table.

No curated element was found by MAP alone. BacFlux ran CONJScan 2.1.0 alongside
ICEscan, while MAP runs ICEscan alone, so the machinery models overlap but are not
identical. ICEVflInd1 rests on CONJScan profile hits alone, reported at low confidence
because no system was assembled. Tn4451 and MTnPi10 carry hits to ICEscan's IME and
CONJScan's MOB models, and MAP's prescan produced no candidate for either.

**Class-assignment caveat for ICEVflInd1.** BacFlux detects the locus but classifies it
`cime_or_island` — integrase + T4CP + T4SS, **no relaxase** — not `ice`. It is
counted as detected because the scorer asks whether a call overlaps the curated
interval, which it does (0.289 recovered). It is *not* a correct ICE class call,
and should not be presented as one.

### 4.1 Per-pilot detection, BacFlux alone

Scored by `benchmark/mobilome/score.py` on the ICE and IME pilot sets:

**ICE pilot: 15/18** — chromosomal 9/11, standalone 6/7.

One of the 18 benchmark records does not contain the expected machinery. `SXT(HN1)`
(AB450045) is a broken benchmark record: the deposit is **19,086 bp carrying 17 CDS**
(verified by counting the Bakta `.faa`), standing in for a ~100 kb element, and contains
**no relaxase and no T4SS**. No machinery-based caller can find a conjugation system
that is not in the sequence, and MAP does not find it either. It is counted as a miss
throughout, so the highest achievable score on this pilot is 17/18; it should be removed
from, or footnoted in, any future version of this benchmark.

The other two misses are real:

- **ICEPsy10** (*P. syringae*, 161 kb) — no relaxase, T4CP or VirB4 anywhere in
  the element. MacSyFinder's own `rejected_candidates.tsv` says
  *"quorum of mandatory genes required (3) is not reached: 0"*.
- **ICEB2** (*Mycoplasma bovis*) — machinery divergent enough that the models do
  not cover it. A genuine model-coverage limitation, and one that will apply to
  any tool built on the same profiles, MAP included.

**IME pilot: 5/12**, median recovered fraction 0.69. The misses concentrate at the
small end (522 bp, 1,248 bp, 5,123 bp), where the size floor and the two-anchor
requirement apply, and at SGI1 (42 kb) and MGIVvuTai1 (19 kb), which produced no
overlapping call at all. IME detection is the weaker half of this module and
should be described that way.

---

## 5. Class assignment

Of the 17 elements both callers found, the classes agree on **14** and differ on
**3**. Every disagreement is the same shape:

| Element | ICEberg curation | BacFlux | MAP |
|---|---|---|---|
| SPI-7 (*S. enterica*) | ICE | **ice** | ime |
| ICE_EfalV583_tRNALys (*E. faecalis*) | ICE | **ice** | ime |
| ICE_Step12228_tRNAser (*S. epidermidis*) | ICE | **ice** | ime |

All three are ICEberg-curated **ICEs** that MAP called **IMEs**; BacFlux agrees
with the curation on all three.

ICE versus IME is the difference between tier 6 (predicted self-transmissible) and
tier 5 (mobilisable, needs a helper). Calling a curated ICE an IME therefore lowers the
predicted transferability of every AMR gene inside it. Neither label establishes
whether a gene is intrinsic or acquired, and neither demonstrates transfer.

**Limits of this comparison.**

1. There is a fourth element, **TR2** (*Streptomyces scabiei*), where **both**
   callers say `ime` against an ICEberg **ICE** curation. It is not counted among
   the three disagreements because the two callers agree with each other. It is a
   shared miss, and it means the score above is "BacFlux agrees with ICEberg on 3
   elements where MAP does not", **not** "BacFlux is right about class and MAP is
   wrong in general".
2. Class was assessed only on elements both callers found, so this is a
   17-element comparison, not 30. With n=3 disagreements, the direction is
   consistent but the sample is small. It should be reported as a consistent
   direction, not a rate.

---

## 6. Call burden, negative controls, boundaries

### 6.1 Call burden

Across all 40 genomes:

| | Calls | Total bp called |
|---|---|---|
| BacFlux | 65 | 2,630,059 |
| MAP | 64 | 2,881,859 |

BacFlux's 65 calls are 63 on the 28 pilot genomes, run with the ICEscan layer on,
and 2 on the 12 negative genomes, run with CONJScan alone. With the ICEscan layer on
all 40 genomes, as in MAP, BacFlux makes 68 calls over 2,687,359 bp: four more calls
than MAP and **6.7% fewer** base pairs, while recovering three more curated elements.
**Neither caller achieves recall by volume** — the higher recall in §4 does not come
from making many more calls.

> **Provenance of these two rows.** BacFlux's row was recomputed over the same 40
> samples after a fix to the att search. Three samples changed: `AL513382`
> 95,548 → 145,727 bp; `GU725392` 2 calls / 57,721 bp → 1 call / 92,237 bp;
> `NC_013929` 200,047 → 272,744 bp. **MAP's row is unchanged** — its run was never
> re-done; its figures describe the same 40 genomes.

### 6.2 Negative controls

Twelve genomes with no curated ICEberg entry, 32,601,100 bp (listed in
`benchmark/mobilome/negative_set.tsv`):

| Genome | BacFlux | MAP |
|---|---|---|
| *B. subtilis* 168 | `ice` 529,362–549,932 (20.6 kb) | `ice` 529,362–549,932 (20.6 kb) |
| *S.* Typhimurium LT2 | `ime` 2,900,443–2,906,075 (5.6 kb) | — |
| *S. aureus* N315 | — | `ime` 36,435–42,455 (6.0 kb) |
| *Halobacterium* sp. NRC-1 | — (3 low-confidence `cime_or_island` with the ICEscan layer on) | — |
| Other eight (incl. *E. coli* K-12, *P. aeruginosa* PAO1) | — | — |

**Two calls each, with BacFlux's default models. The difference between the callers
is one call each**, and they are on different genomes. BacFlux's column comes from
CONJScan alone; with the ICEscan layer on, as in MAP, BacFlux adds three low-confidence
`cime_or_island` calls on the archaeon, built from integrase and type IV secretion
profile hits with no assembled system.

**The *B. subtilis* 168 call is the same element, found twice.** Both pipelines
return **base-pair-identical coordinates**, the same 60 bp flanking repeat, the
same MOBT relaxase and the same MPF type — reached by two *different* boundary
algorithms (MAP: vmatch plus an ICEfinder2-derived refiner; BacFlux: its own
exact maximal-repeat search). Because boundary refinement is the genuinely
independent layer (§0), this agreement is real evidence, and what it evidences is
that **ICEBs1 is a real element that ICEberg has not catalogued**: 20,571 bp
against a ~20.5 kb element, all four anchor classes present, and the boundary
resolved by a tRNA-anchored att site at tRNA-Leu(gag) — ICEBs1 integrates at
*trnS-leu2*. Bakta independently names the integrase "ICEBs1 integrase".

BacFlux's other call, on *S.* Typhimurium LT2, is an integrase beside a protein
Bakta annotates only as "DNA-binding protein" but which the T4SS_MOBM profile
hits at i-eval 1.8 × 10⁻⁸⁴ over 98.6% of the profile. Whether that locus is a
bona fide IME is a domain question this benchmark cannot settle.

**Two caveats, both mandatory when quoting these numbers.**

> **These are candidate calls, NOT a false-positive rate.**
> "Absent from ICEberg" is not "contains no mobile element". ICEberg is a curated
> catalogue, not an exhaustive one. The counts describe call burden on the selected control genomes.
> Interpreting individual calls still requires evidence beyond catalogue absence.

> **Every call in a negative control must be examined, never merely counted.**
> This set originally justified including *B. subtilis* 168 on the grounds that
> "ICEBs1 is absent from 168 itself". That is **false** — ICEBs1 was *discovered*
> in strain 168 (Auchtung et al. 2005). Had the calls been tallied instead of
> inspected, a textbook-correct detection would have been recorded as a false
> positive, and the pipeline would have been "fixed" to stop making it.

### 6.3 Boundaries

Measured on the elements both callers found (chromosomal n=9, standalone n=5):

| Metric | BacFlux | MAP |
|---|---|---|
| Chromosomal ICEs, median \|start offset\| | **5,141 bp** | 7,427 bp |
| Chromosomal ICEs, median \|end offset\| | **14,912 bp** | 23,745 bp |
| Standalone ICEs, median recovered fraction | 0.94 | 0.94 |

The callers match on standalone span recovery, and BacFlux is closer on both
chromosomal offsets. On **ICEEc2** (GU725392) they land on **exactly the same
interval**, 27–92,263, from independent algorithms.

!!! warning "This comparison changed the code it was measuring"

    An earlier version of BacFlux's att search scored worse than MAP on both
    chromosomal offsets, and this comparison is what prompted two fixes: a ranking
    that put raw repeat length above tRNA anchoring (on SPI-7 a 51 bp repeat in
    ordinary sequence was beating the real 24 bp pair at tRNA-Phe), and a flank
    window of 30 kb when SPI-7's true *attR* sits 5.8 kb outside it. The numbers
    above are measured after those fixes.

    Stated because tuning against a benchmark and then reporting that benchmark is
    a methodological weakness, and a reader is entitled to know it happened. Two
    things limit it: the fixes were to defects identifiable without the comparison,
    and the conservative policy was **not** loosened to achieve them — a de novo
    repeat pair is still *reported and not applied*, and both recovered boundaries
    are tRNA-anchored, so they apply under the pre-existing rule.

**Neither fix was a better repeat finder.** Of the four elements where MAP's boundary
beat ours, **zero** were cases where MAP found a repeat BacFlux could not. Two were
repeats BacFlux found and then discarded through the defects above; two were
candidate-extent differences with no repeats on either side. That is what §8 tests.

> ### Denominator of the 0.94 figure { #honesty-note-what-the-094-figure-covers }
>
> The standalone median of **0.94** is computed over the **five** standalone ICEs
> that **both** callers found — it is a like-for-like head-to-head figure and MAP's
> 0.94 is computed over exactly the same five.
>
> Across **all six** standalone ICEs BacFlux detected, the median recovered
> fraction is **0.669**, because the sixth is ICEVflInd1 at 0.289 — an element MAP
> did not find at all. Detecting one extra, hard element *lowers* your own median.
> When quoting a pilot-wide number for BacFlux alone, use **0.669**; when quoting
> a head-to-head number, use 0.94 and say it is on the five shared elements.
> The BacFlux-only ICE-pilot summary is therefore: chromosomal median recovered
> **0.476** (n=9), standalone median recovered **0.669** (n=6).

Both callers remain far from the curated coordinates on large chromosomal ICEs.
On TR2 both recover 0.055 of a 154 kb element; on IMEMlNZP2037-1 both recover
0.074 of 110 kb; on ICEMsp.M1D BacFlux recovers 0.198 and MAP 0.080. Large
chromosomal element delimitation is an open problem for both, and the module's
own output flags this — those calls carry `boundary_method=none`, meaning the
machinery span was reported because no att pair was resolvable.

---

## 7. Limitations

Read these as constraints on every number above.

1. **The callers share a detection engine.** §0. Detection agreement is based partly on
   shared models. Never cite it as independent confirmation.
2. **Small n.** 30 curated elements; 17 in the class comparison; 9 chromosomal and
   5 standalone in the boundary comparison. Medians over five to nine values move
   a long way on one element. Report directions, not rates, and never a
   significance claim.
3. **ICEberg coordinates are the truth standard, and they are neither perfect nor
   exhaustive.** One record in this set is demonstrably broken (SXT(HN1), §4.1),
   and one negative-control genome contains a real element ICEberg lacks (ICEBs1,
   §6.2). The first counts as a miss for both callers; the second scores a correct
   call as a false positive.
4. **The negative control counts candidate calls; it is not a false-positive
   rate.** §6.2.
5. **`score.py` selects the largest-overlap call.** Where a caller emits several
   calls on one genome, the metric rewards the biggest overlapping interval. A
   caller that emitted one enormous call per genome would score well on
   `recovered_fraction` and badly on nothing the scorer measures. This is why the
   call-burden table (§6.1) must be reported alongside recall, and it applies
   symmetrically to both callers.
6. **Closed reference genomes, not draft assemblies.** These are single-contig
   sequences. BacFlux's short-read entry point produces fragmented SPAdes assemblies
   where IS elements themselves cause the contig breaks, and fragmentation can reduce
   recovery or boundary resolution. Nothing in this comparison measures that. Any
   element spanning contigs is capped at `low` confidence by design.
7. **Different gene callers.** MAP uses Prodigal, BacFlux uses Bakta. That is a
   real difference contained in the comparison rather than controlled for, and it
   affects both the protein set the shared HMMs search and the integrase
   annotation regex.
8. **One run per genome, no replicates**, and three genomes carry run caveats
   (§2.4).
9. **MAP was run with functional branches disabled.** Fair for a boundary
   comparison, but it means nothing here evaluates MAP's full output.
10. **Bacterial models only.** ICEscan and CONJscan are bacterial. The archaeal
    genome is a construction check, not coverage. With the ICEscan layer on,
    BacFlux makes three low-confidence profile-only calls on it.

---

## 8. Why vmatch was evaluated and not adopted

MAP refines boundaries with Vmatch; BacFlux implements exact maximal-repeat search
in Python (`workflow/scripts/80_mobilome/att_search.py`). On the four elements where
MAP's boundaries were initially closer, MAP found no repeat that BacFlux missed. Two
were repeats BacFlux had found and then discarded through ranking and search-window
defects, since fixed; the other two were differences in candidate extent with no
repeat on either side (§6.3).

A separate implementation using a suffix array, LCP array and cross-flank MEM
enumeration was compared with BacFlux on **52 genomes, 246 att searches, 180 distinct
flank windows and 1,033 repeats**. The repeat sets were identical on every window.
The Vmatch executable itself was not run: on those windows, a search with Vmatch's
maximal-repeat semantics finds no repeat that BacFlux's search misses. Agreement
between two implementations of one definition does not validate an element boundary.

Vmatch is packaged in Bioconda (2.3.1), but its recipe declares `license: Unknown /
OTHER` and the upstream site could not be reached when this was checked in July 2026,
so its terms could not be established. With no measured gain and unverifiable terms,
BacFlux does not use it.

---

## 9. Summary

| Question | Answer |
|---|---|
| Are these independent callers? | **No.** Shared MacSyFinder + ICEscan detection engine. Boundaries and clustering are independent. |
| Detection, 30 curated elements | BacFlux 20/30, MAP 17/30; MAP found nothing BacFlux missed |
| Class, 17 shared elements | Agree on 14; all 3 differences are ICEberg ICEs that MAP called IMEs |
| Call burden, 40 genomes | BacFlux 65 calls / 2.63 Mb (68 / 2.69 Mb with ICEscan on all 40); MAP 64 calls / 2.88 Mb |
| Negative controls, 32.6 Mb | 2 candidate calls each (BacFlux 5 with ICEscan on); identical coordinates on the one shared call |
| Boundaries | Comparable — chromosomal median \|start offset\| 5,141 bp vs MAP 7,427; standalone recovery 0.94 vs 0.94. §6.3 records that this comparison prompted two fixes to our att search |
| Is vmatch required by BacFlux? | No. Two implementations of the intended repeat search agreed on the 180 tested windows; this was not a direct Vmatch executable benchmark. |

Nothing here supports a claim that either pipeline is better than the other. It
supports three narrower claims: BacFlux's caller recovers the curated elements MAP
recovers plus three more; it assigns the ICE/IME class in line with ICEberg's
curation in the three cases where the two callers differ; and its boundary
refinement is comparable to a vmatch-based refiner without
requiring vmatch.

Every mobility statement derived from these calls remains a **prediction**. An
element reported as **predicted self-transmissible** has conjugation machinery
that looks intact; it has not been shown to transfer. The confirmatory experiment
is a **filter or broth mating assay**, and no output of this module substitutes
for one.

---

## 10. Citations and credit

The tools and resources this comparison depends on, credited plainly. Reading
EBI's pipeline is how we found the ICEscan model set, which materially improved
BacFlux's IME and AICE handling; that debt should be acknowledged in any
publication.

- **EBI Mobilome Annotation Pipeline** — ebi-metagenomics/mobilome-annotation-pipeline,
  v5.0.0, **Apache-2.0**. EMBL-EBI Microbiome Informatics team.
  https://github.com/ebi-metagenomics/mobilome-annotation-pipeline
- **ICEfinder2** — Wang M., *et al.* ICEberg 3.0: functional categorization and
  analysis of the integrative and conjugative elements in bacteria. *Nucleic Acids
  Research* (2024). ICEfinder2 provides the ICEscan models and the boundary-
  refinement approach MAP's ICEfinder2-lite implements.
- **CONJScan / MacSyFinder** — Cury J., Abby S.S., Doppelt-Azeroual O.,
  Néron B., Rocha E.P.C. Identifying conjugative plasmids and integrative
  conjugative elements with CONJscan. *Methods in Molecular Biology* (2020);
  and Abby S.S., Néron B., Ménager H., Touchon M., Rocha E.P.C. MacSyFinder:
  a program to mine genomes for molecular systems with an emphasis on CRISPR-Cas
  systems. *PLoS ONE* (2014).
- **Conjugation system classification** — Coluzzi C., Garcillán-Barcia M.P.,
  de la Cruz F., Rocha E.P.C. Evolution of plasmid mobility: origin and
  fate of conjugative and non-conjugative plasmids. *Molecular Biology and
  Evolution* (2022).
- **ICEberg** (curated coordinates used as truth) — Wang M. *et al.*, *NAR* (2024).
- **ICEBs1 discovery** (the negative-control lesson, §6.2) — Auchtung J.M.,
  Lee C.A., Monson R.E., Lehman A.P., Grossman A.D. Regulation of a *Bacillus
  subtilis* mobile genetic element by intercellular signaling and the global DNA
  damage response. *PNAS* (2005).

### Licensing position

No code from the EBI pipeline was copied or adapted. Three of its scripts —
`bin/ice_boundary_refinement.py`, `bin/map_tools/icefinder_process.py` and
`bin/prescan_to_fasta.py` — derive from ICEfinder2 and carry CC BY-NC-SA 4.0 headers.
Their NonCommercial and ShareAlike conditions are incompatible with distribution under
BacFlux's MIT licence, and making such code optional would not change that.

CONJScan and ICEscan models are fetched separately when enabled and are not bundled
with BacFlux. Their attribution, NonCommercial and share-alike terms remain relevant
to use and distribution. See [Licensing](about/licensing.md).

Users can analyse the same assembly with MAP for a broader inventory. Check
terms relevant to its outputs and source databases separately from the licence
of BacFlux's implementation.
