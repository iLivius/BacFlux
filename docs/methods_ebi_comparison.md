# Comparing BacFlux's ICE/IME caller against the EBI Mobilome Annotation Pipeline

The obvious question about any new caller is whether an established one does better.
This page answers it directly: BacFlux's mobilome module and the EBI Mobilome
Annotation Pipeline, run on the same genomes and scored by the same code.

Read section 0 first. The two callers share ancestry, so this is **not** an
independent check — treating it as one would overstate what agreement between them
proves.

Numbers current as of BacFlux commit `4a93d89` (branch `release/v2.0.0`), and every
figure was recomputed from the run artefacts while writing, not carried over from a
draft.

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

## 0. Read this first: the two callers are not independent

**The most important caveat in this document, stated before any result.**

BacFlux's ICE/IME caller and the EBI Mobilome Annotation Pipeline (MAP) both find
conjugation machinery by running **MacSyFinder with the ICEscan model set**. MAP
does it inside its `ICEFINDER2_LITE` subworkflow; BacFlux does it in the
`conjscan` / `icescan` rules of `workflow/rules/shared/80_mobilome.smk`. The
versions differ slightly — MacSyFinder 2.1.4 in MAP's container, 2.1.6 in
BacFlux's conda environment (`workflow/envs/macsyfinder.yaml`) — but the profile
HMMs and the system definitions that decide *"is there a conjugation system here
at all"* are the same models.

The consequence must not be glossed over:

> **Agreement at the detection layer is partly tautological.** When both callers
> find the same element, that is substantially the same engine agreeing with
> itself. It is *not* two independent lines of evidence, and it must never be
> reported as independent confirmation of detection.

What *is* genuinely independent, and therefore what this comparison can actually
measure:

| Layer | Shared or independent? |
|---|---|
| Machinery detection (does a conjugation system exist here) | **Shared** — same models, same engine |
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
transferable is that element*, placing each on a six-rung ladder from
"chromosomal, intrinsic candidate" to "**predicted** self-transmissible". The
regulatory framing is EFSA's intrinsic-versus-acquired distinction for strains
entering the food and feed chain. The module is **off by default** and is one
stage of a workflow whose main business is assembly, decontamination, taxonomy,
annotation and AMR calling.

**These are complementary scopes.** MAP goes wider across element classes and is
built for community-level data; BacFlux goes narrower but ties every call back to
a specific resistance gene in a specific isolate and reports a mobility tier with
a confidence level. Neither is a reduced version of the other, and this document
should not be read as ranking them. MAP's own documentation notes that its
current release does not run gene-level AMR association, which is precisely the
step BacFlux's module exists to perform.

A practical consequence for users: for a full
mobilome inventory of an isolate, or for anything metagenomic, MAP is the right
tool and BacFlux writes `{sample}_contigs.fna` and `{sample}.gbk` ready to feed it.
Consuming MAP's `mobilome.gff.gz` carries no licensing consequence for BacFlux
(§10).

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
thermophile, four full-size reference strains (*E. coli* K-12 MG1655, *B.
subtilis* 168, *P. aeruginosa* PAO1, *S. aureus* N315, *S.* Typhimurium LT2) and
an archaeon (*Halobacterium* sp. NRC-1). The archaeon is a **construction
check**: CONJscan and ICEscan are bacterial models, so any call there would be
spurious by definition.

Both callers were given **the same FASTA files** from
`phase7_benchmark/genomes/`. No re-assembly, no re-annotation on our side of the
fence — MAP does its own gene calling (Prodigal), BacFlux uses Bakta, which is one
of the differences the comparison contains rather than controls for.

### 2.2 Scoring

Both callers' output was reduced to one table per genome,
`work/{sample}/ice_elements.tsv`, with the columns `phase7_benchmark/score.py`
reads. `score.py` was then run **unmodified** over both. It:

- takes the curated ICEberg interval as truth;
- among a caller's calls on that genome, selects the one with the **largest
  overlap** with the curated interval;
- reports start offset, end offset, and `recovered_fraction` = overlapping bp ÷
  curated length.

Using the same scoring code on both sides is the point. Any bias in `score.py` —
and there is at least one, see §7 — applies equally to both callers.

### 2.3 Which MAP, and the two environment fixes needed to run it

MAP **v5.0.0** (`ebi-metagenomics/mobilome-annotation-pipeline`, repo checkout at
`3aa408d`), run with `-profile singularity`. Recorded tool versions inside the
containers: MacSyFinder 2.1.4, ISEScan 1.7.3, IntegronFinder 2.0.6, geNomad
1.11.1, HMMER 3.4, Aragorn 1.2.41.

Two host-specific problems had to be solved first. Both are recorded in
`ebi_map/nf_env.sh` and `ebi_map/run_map.sh`; anyone repeating this will hit them.

**(a) Apptainer instead of the host's Singularity.** This machine carries
Singularity **2.6.1**, which predates the SIF container format that all current
biocontainers use — it cannot run them at all. Apptainer (the renamed, modern
continuation of Singularity) was installed into a private conda environment and
put ahead of the system binary on `PATH`. It runs rootless here because the kernel
has unprivileged user namespaces enabled, and the conda package ships a
`singularity` → `apptainer` symlink, which lets Nextflow's stock `singularity`
profile work unmodified.

**(b) `NXF_VER=24.10.6`.** On Nextflow 26.04.6 — the current release, and what a
fresh conda install gives you — MAP v5.0.0 dies inside its own ICEfinder2-lite
subworkflow:

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

**(c) Functional branches switched off** (`ebi_map/skip_functional.config`).
AMRFinderPlus, DeepARG, RGI-CARD, PathoFact2, SanntiS, GECCO and antiSMASH were
disabled: they answer a different question, each needs its own multi-GB database,
and they contribute nothing to a boundary comparison. geNomad and CheckV are
**not** skippable in MAP and feed its overlap integrator, so they ran. Note these
have to be set in a config file, not as `--skip_x` on the command line: MAP
validates parameters with nf-schema in strict mode and a bare CLI `--skip_x`
arrives as the string `"true"`, failing the boolean type check.

### 2.4 Completeness of the MAP run

**40/40 genomes reached `COMPLETED`** (`ebi_map/map_status_40.tsv`). Median
wall-clock 30.3 min per genome (range 23.4–39.9). Twelve genomes produced no
`ices.tsv` because no ICE candidate survived MAP's prescan; that is a real
result — no call — not a failure, and it is scored as zero calls.

Three run-level caveats, disclosed rather than buried:

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
All three are written into `ebi_map/scoring/convert_map_to_ice_elements.py` at the
point where they take effect. The converter does **no** filtering, no
re-thresholding and no coordinate arithmetic — only translation.

### 3.1 The SO-term `integron` trap

**This is the single most important thing to know before repeating this
comparison, and it fails silently.**

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
ICEs + 27 IMEs + 5 AICEs (`ebi_map/scoring/map_conversion_audit.tsv`). All five
AICE calls fall in *Streptomyces scabiei*, which is biologically correct — AICEs
are the actinomycete class.

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

**MAP found nothing BacFlux missed.** Given §0, this is the strongest form the
detection result can take: the difference cannot be explained by a better
machinery model, because the machinery models are shared. It is attributable to
what sits around the shared engine — candidate seeding, the class-aware size
floor, and BacFlux's decision to report low-quorum machinery at low confidence
rather than discard it.

**Honesty note on ICEVflInd1.** BacFlux detects the locus but classifies it
`cime_or_island` — integrase + T4CP + T4SS, **no relaxase** — not `ice`. It is
counted as detected because the scorer asks whether a call overlaps the curated
interval, which it does (0.289 recovered). It is *not* a correct ICE class call,
and should not be presented as one.

### 4.1 Per-pilot detection, BacFlux alone

From `phase7_benchmark/results_final.tsv` and `ime_results_final.tsv`:

**ICE pilot: 15/18** — chromosomal 9/11, standalone 6/7.

The **true ceiling is 17/18**, not 18/18. `SXT(HN1)` (AB450045) is a broken
benchmark record: the deposit is **19,086 bp carrying 17 CDS** (verified by
counting the Bakta `.faa`), standing in for a ~100 kb element, and contains **no
relaxase and no T4SS**. No machinery-based caller can find a conjugation system
that is not in the sequence, and MAP does not find it either. It should be
removed from, or footnoted in, any future version of this benchmark.

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

This is not a cosmetic difference. ICE versus IME **is** the
self-transmissible-versus-needs-a-helper distinction, which is rungs 6 and 5 of
the mobility ladder and the exact axis the EFSA intrinsic/acquired framing rests
on. Calling a self-transmissible element "mobilisable with a helper" understates
transfer risk.

**Two honesty notes.**

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

One call apart, and BacFlux claims **8.7% fewer** base pairs while recovering
three more curated elements. **Neither caller achieves recall by volume** — the
higher recall in §4 does not come from making many more calls.

> **Provenance of these two rows.** `ebi_map/scoring/call_burden.tsv` was written
> at 09:29 on 2026-07-31 from a snapshot of our calls taken at 09:26 — *before*
> `4a93d89` regenerated them at ~11:12. Its BacFlux row therefore reads 66 calls /
> 2,472,667 bp, which is stale. The figures above are recomputed from the live
> `phase7_benchmark/work/*/ice_elements.tsv` over the same 40 samples. Three
> samples changed: `AL513382` 95,548 → 145,727 bp; `GU725392` 2 calls / 57,721 bp
> → 1 call / 92,237 bp; `NC_013929` 200,047 → 272,744 bp. **MAP's row is
> unchanged** — its run was never re-done, so the comparison remains like for
> like.

### 6.2 Negative controls

Twelve genomes with no curated ICEberg entry, 32,601,100 bp
(`ebi_map/scoring/negative_side_by_side.tsv`):

| Genome | BacFlux | MAP |
|---|---|---|
| *B. subtilis* 168 | `ice` 529,362–549,932 (20.6 kb) | `ice` 529,362–549,932 (20.6 kb) |
| *S.* Typhimurium LT2 | `ime` 2,900,443–2,906,075 (5.6 kb) | — |
| *S. aureus* N315 | — | `ime` 36,435–42,455 (6.0 kb) |
| Other nine (incl. *E. coli* K-12, *P. aeruginosa* PAO1, *Halobacterium*) | — | — |

**Two calls each. The difference between the callers is one call each**, and they
are on different genomes.

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

> **This is an upper bound of 2 candidate calls, NOT a false-positive rate.**
> "Absent from ICEberg" is not "contains no mobile element". ICEberg is a curated
> catalogue, not an exhaustive one. The only quantity here that is genuinely
> comparable is the *difference* between the two columns, because the imperfect
> negative applies equally to both callers.

> **Every call in a negative control must be examined, never merely counted.**
> This set originally justified including *B. subtilis* 168 on the grounds that
> "ICEBs1 is absent from 168 itself". That is **false** — ICEBs1 was *discovered*
> in strain 168 (Auchtung et al. 2005). Had the calls been tallied instead of
> inspected, a textbook-correct detection would have been recorded as a false
> positive, and the pipeline would have been "fixed" to stop making it.

### 6.3 Boundaries — where MAP was better, and by how much

**Before BacFlux commit `4a93d89`, MAP's boundaries were better than ours. That
is what prompted the fixes.**

Measured on the elements both callers found
(`ebi_map/scoring/boundary_agreement.tsv`; chromosomal n=9, standalone n=5):

| Metric | BacFlux **before** `4a93d89` | BacFlux **at** `4a93d89` | MAP |
|---|---|---|---|
| Chromosomal ICEs, median \|start offset\| | 14,294 bp | **5,141 bp** | 7,427 bp |
| Chromosomal ICEs, median \|end offset\| | 34,239 bp | **14,912 bp** | 23,745 bp |
| Standalone ICEs, median recovered fraction | 0.59 | **0.94** | 0.94 |

After the fixes BacFlux matches MAP on standalone span recovery and is closer on
both chromosomal offsets. The single clearest case: on **ICEEc2** (GU725392) the
two callers now land on **exactly the same interval**, 27–92,263, from
independent algorithms.

**The fixes did not come from a better repeat finder.** Of the four elements where
MAP's boundary beat ours, **zero** were cases where MAP found a repeat BacFlux
could not. Two were repeats BacFlux found and then discarded through its own
defects; two were candidate-extent differences with no repeats on either side. The
two defects — a ranking that put raw repeat length above tRNA anchoring (on SPI-7
a 51 bp repeat in ordinary sequence was beating the real 24 bp pair at tRNA-Phe),
and a flank window of 30 kb when SPI-7's true attR sits 5.8 kb outside it — are
documented in the commit message of `4a93d89`.

The conservative policy was **not** loosened to achieve this: a de novo
repeat pair is still *reported and not applied*. Both recovered boundaries are
tRNA-anchored and so apply under the pre-existing rule.

> ### Honesty note: what the 0.94 figure covers
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

1. **The callers share a detection engine.** §0. Detection agreement is partly
   tautological. Never cite it as independent confirmation.
2. **Small n.** 30 curated elements; 17 in the class comparison; 9 chromosomal and
   5 standalone in the boundary comparison. Medians over five to nine values move
   a long way on one element. Report directions, not rates, and never a
   significance claim.
3. **ICEberg coordinates are the truth standard, and they are neither perfect nor
   exhaustive.** One record in this set is demonstrably broken (SXT(HN1), §4.1),
   and one negative-control genome contains a real element ICEberg lacks (ICEBs1,
   §6.2). The first counts as a miss for both callers; the second scores a correct
   call as a false positive.
4. **The negative control gives an upper bound on candidate calls, not a
   false-positive rate.** §6.2.
5. **`score.py` selects the largest-overlap call.** Where a caller emits several
   calls on one genome, the metric rewards the biggest overlapping interval. A
   caller that emitted one enormous call per genome would score well on
   `recovered_fraction` and badly on nothing the scorer measures. This is why the
   call-burden table (§6.1) must be reported alongside recall, and it applies
   symmetrically to both callers.
6. **Closed reference genomes, not draft assemblies.** These are single-contig
   sequences. BacFlux's short-read entry point produces fragmented SPAdes
   assemblies where IS elements themselves cause the contig breaks, and calls
   there will be materially worse. Nothing in this comparison measures that. Any
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
    genome is a construction check, not coverage.

---

## 8. Why vmatch was evaluated and not adopted

MAP refines boundaries with **vmatch**; BacFlux does not. When MAP's boundaries
were measurably better (§6.3), the obvious hypothesis was that vmatch was the
reason. It was tested rather than assumed, and the hypothesis was wrong twice
over.

**First, the two searches are the same search.** BacFlux's att search has not used
`blastn` since commit `9c05c3e`; it computes **exact maximal repeats in standard-
library Python** (`workflow/scripts/80_mobilome/att_search.py`), which is vmatch's
own semantics. To test that properly rather than argue it, a second search was
built on vmatch's actual data structure — prefix-doubling suffix array, Kasai LCP
array, cross-flank MEM enumeration — and run against ours on every real flank
window in the benchmark:

> **35/35 windows, 122 repeats, exact set equality.**

Two searches over the same DNA return the same set of repeats. vmatch would find
nothing BacFlux does not already find. The same measurement also disposes of a
related idea: **treating agreement between two repeat searches as evidence is empty**,
because when the two searches have identical semantics the agreement is
guaranteed and therefore carries no information.

**Second, the boundary gap was not a repeat-finding gap at all.** Of the four
elements where MAP's boundary beat ours, none involved a repeat MAP found and we
could not (§6.3).

**The licensing position is secondary but real.** An earlier version of the spec
rejected vmatch as "not on bioconda", which was **wrong** — `bioconda/vmatch
2.3.1` exists for linux-64 and osx-64. The rejection stands on different ground:
the recipe declares `license: Unknown / OTHER` and vmatch.de is unreachable, so
its terms **cannot be established**. For an MIT workflow that is weaker ground
than a known-restrictive licence, because there is nothing to comply with. With
no measured benefit on one side and unverifiable terms on the other, vmatch is not
adopted.

> **Reproducibility gap, stated rather than hidden.** The 35/35 set-equality
> measurement is recorded in the commit message of `4a93d89`, but the diagnostic
> script that produced it is not retained in the repository or the benchmark tree.
> Before this appears in a publication it should be re-run and the script and its
> output archived alongside the other validation artefacts.

---

## 9. Summary

| Question | Answer |
|---|---|
| Are these independent callers? | **No.** Shared MacSyFinder + ICEscan detection engine. Boundaries and clustering are independent. |
| Detection, 30 curated elements | BacFlux 20/30, MAP 17/30; MAP found nothing BacFlux missed |
| Class, 17 shared elements | Agree on 14; all 3 differences are ICEberg ICEs that MAP called IMEs |
| Call burden, 40 genomes | BacFlux 65 calls / 2.63 Mb; MAP 64 calls / 2.88 Mb |
| Negative controls, 32.6 Mb | 2 candidate calls each; identical coordinates on the one shared call |
| Boundaries before `4a93d89` | **MAP better** — chromosomal median \|start offset\| 7,427 bp vs our 14,294 |
| Boundaries at `4a93d89` | Comparable — 5,141 bp vs MAP 7,427; standalone recovery 0.94 vs 0.94 |
| Is vmatch needed? | No — exact set equality with our search on 35/35 flank windows |

Nothing here supports a claim that either pipeline is better than the other. It
supports three narrower claims: BacFlux's caller recovers the curated elements MAP
recovers plus three more; it assigns the ICE/IME class in line with ICEberg's
curation in the three cases where the two callers differ; and its boundary
refinement, after `4a93d89`, is comparable to a vmatch-based refiner without
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
- **ICEfinder2** — Wang J., *et al.* ICEberg 3.0: functional categorization and
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
- **ICEberg** (curated coordinates used as truth) — Wang J. *et al.*, *NAR* (2024).
- **ICEBs1 discovery** (the negative-control lesson, §6.2) — Auchtung J.M.,
  Lee C.A., Monson R.E., Lehman A.P., Grossman A.D. Regulation of a *Bacillus
  subtilis* mobile genetic element by intercellular signaling and the global DNA
  damage response. *PNAS* (2005).

### Licensing position

- **No code from EBI's pipeline was copied.** Three of its scripts
  (`bin/ice_boundary_refinement.py`, `bin/map_tools/icefinder_process.py`,
  `bin/prescan_to_fasta.py`) are ICEfinder2-derived and carry **CC BY-NC-SA 4.0**,
  which is non-commercial and share-alike and therefore incompatible with
  BacFlux's MIT licence. Share-alike attaches to distributed source, not to
  execution, so a default-off flag would not cure it. Reading that code to
  understand an approach is permitted and is how the ICEscan model set was found;
  the separation between reading and copying was kept real.
- **CONJScan and ICEscan model sets are both CC BY-NC-SA 4.0** (Institut
  Pasteur/CNRS; ICEscan is an ICEfinder2 fork of CONJScan 2.0.1 by the same
  authors). They are **fetched at runtime from a configured URL and never
  vendored**, the mobilome module is **off by default**, and
  `docs/about/licensing.md` carries the non-commercial notice — the same pattern
  BacFlux already uses for every
  licence-encumbered database.
- **Consuming MAP's `mobilome.gff.gz` output carries no licensing consequence**
  for BacFlux. Users who want a full mobilome census are pointed at MAP in the
  README, and BacFlux writes `{sample}_contigs.fna` and `{sample}.gbk` ready for it.
