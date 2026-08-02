# Tuning the mobilome module

*Destined for the MkDocs site. A companion to `mobilome_worked_example.md`: that
page teaches you to read the output, this one teaches you to change it.*

*Every number on this page was measured from artefacts on disk
(`/media/data/antonielli_dir/BacFlux_v2_validation/`) or read out of the source,
on 2026-07-31. Where a knob has never been measured, it says so rather than
offering a plausible-sounding recommendation.*

---

## How to use this page

You are here because you looked at a result and wanted something different. So
the page is organised around what you **saw**, not around the parameter list:

| I want to… | Go to |
|---|---|
| find out why something was dropped | [§1 Start at the audit file](#1-start-at-the-audit-file-not-at-the-config) |
| write a filter that survives the next release | [§1.1 Never hard-code a column number](#11-never-hard-code-a-column-number-in-a-result-table) |
| know where a given knob physically lives | [§2 Where the knobs live](#2-where-the-knobs-live) |
| look a parameter up | [§3 The inventory](#3-the-inventory) |
| tighten everything for a dossier | [§4.1 More conservative](#41-more-conservative) |
| loosen everything for a survey | [§4.2 Wider net](#42-wider-net) |
| match a symptom to a knob | [§5 Symptom → knob](#5-symptom--knob) |
| know whether a knob can fix this at all | [§6 What is not tunable](#6-what-is-not-tunable-and-why) |
| see it done | [§7 Worked examples](#7-worked-examples) |
| avoid the traps we fell into | [§8 Heads-up](#8-heads-up-the-traps-we-hit) |

---

## 1. Start at the audit file, not at the config

Every filtering decision in this module writes a row with a **reason**. That
audit trail is the whole tuning story: it tells you *why* something was dropped,
and therefore *which* knob to reach for. Changing a parameter before reading it
is guessing.

Per sample, under `08.mobilome/{sample}/`:

| File | Written by | What it explains |
|---|---|---|
| `{sample}_ice_discarded.tsv` | `conjscan_to_ice.py` | every ICE/IME candidate dropped, downgraded or flagged, and why |
| `{sample}_amr_mobility_audit.tsv` | `colocalise.py` | every AMR gene that did **not** get the context or tier you expected |
| `{sample}_is_discarded.tsv` | `isescan_to_table.py` | IS calls dropped or flagged |
| `{sample}_named_elements_discarded.tsv` | `name_transposons.py` | TnCentral hits refused a name (tier 4 layer) |
| `{sample}_ice_naming.tsv` | `name_ice_elements.py` | ICEberg hits refused a name |
| `{sample}_is_copy_number_audit.tsv` | `isosdb_copy_number.py` | database entries excluded from the copy-number estimate |

Four of the six share one shape — `sample, contig, start, end, action, reason,
detail`. Two do not:

- `{sample}_amr_mobility_audit.tsv` keys on the gene rather than an interval:
  `sample, contig, amr_gene, amr_start, amr_end, mobility_tier, mge_context,
  decision, reason, detail`.
- `{sample}_is_copy_number_audit.tsv` keys on a database entry, and has no
  coordinates at all: `sample, element, action, reason, detail`.

In every one of them the **last three** columns work the same.
`action`/`decision` says what happened, `reason` is a short token you filter on,
and `detail` carries the numbers in plain language — usually including the
threshold that was applied and the value that failed it. Only the column
*numbers* differ, which is why the recipes below say which column to cut.

```bash
# What happened to the ICE/IME candidates in this sample?
cut -f5,6 08.mobilome/S1/S1_ice_discarded.tsv | sort | uniq -c | sort -rn

# Why did this AMR gene not get a mobile-element context?
# (this file has a wider key block: decision is column 8, reason 9, detail 10)
awk -F'\t' '$8=="no_mge_context" || $8=="tier_not_raised"' \
    08.mobilome/S1/S1_amr_mobility_audit.tsv | cut -f3,9,10

# Everything that was actually thrown away, across all samples
awk -F'\t' 'FNR>1 && $5=="dropped"' 08.mobilome/*/*_ice_discarded.tsv
```

**The `action` words are not the same in every file**, which trips people up
when they grep across all six. In `{sample}_ice_discarded.tsv` and
`{sample}_is_discarded.tsv` a thing that is gone from the table says `dropped`;
in `{sample}_named_elements_discarded.tsv` and `{sample}_is_copy_number_audit.tsv`
the same event says `discarded`. The full vocabulary, read out of the source:

| File | `action` / `decision` values it can write |
|---|---|
| `{sample}_ice_discarded.tsv` | `dropped`, `kept_flagged`, `evidence_recorded`, `boundaries_resolved`, `not_applicable`, `row_skipped`, `input_missing` |
| `{sample}_is_discarded.tsv` | `dropped`, `kept_flagged`, `input_missing` |
| `{sample}_named_elements_discarded.tsv` | `discarded`, `not_applicable`, `summary` |
| `{sample}_ice_naming.tsv` | `named`, `kept_flagged`, `not_applicable` |
| `{sample}_is_copy_number_audit.tsv` | `discarded`, `not_applicable`, `summary` |
| `{sample}_amr_mobility_audit.tsv` | `input_missing`, `row_skipped`, `no_mge_context`, `tier_not_raised`, `evidence_recorded` |

`dropped`/`discarded` means the thing is gone. `kept_flagged` means it is in the
table but something about it was hedged — usually a confidence cap.
`evidence_recorded`, `boundaries_resolved`, `not_applicable` and `summary` are
neither: they are notes, so that a decision is never invisible. `row_skipped`
and `input_missing` are about the *inputs*, not the biology — a hit that could
not be joined to the annotation, or a file that was not there.

### 1.1 Never hard-code a column number in a RESULT table

The audit layouts above are stable, and the recipes on this page cut them by
number. The two **result** tables are a different matter.
`{sample}_amr_mobility.tsv` and `{sample}_ice_candidates.tsv` have both gained
columns as the module grew, so tables kept from different builds do not line up.
Checked on 2026-07-31: the positive-control table in `docs/validation/` has 44
columns and no `named_element`; a clinical run from 2026-07-28 has 46; the
current code writes 44 in a third order. The ICE tables are worse — the benchmark
ones have 51 columns where the code writes 49, and `evidence_level` is missing
from the clinical and positive-control copies altogether. Every recipe on this
page that touches a result table therefore looks its columns up **by name**:

```bash
awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) col[$i] = i; print; next}
            $col["confidence"] == "high"' any_result_table.tsv
```

One caveat with that idiom: if the column is **absent**, `col[...]` is 0 and `$0`
is the whole line, so the test quietly passes for every row instead of failing
loudly. On an old table that is how a filter silently stops filtering. Check the
header first — `head -1 file | tr '\t' '\n' | nl` — when the row count looks too
generous.

---

## 2. Where the knobs live

There are three places, and they differ in how much you are signing up for.

### 2.1 Layer A — `config/config.yaml`, `mobilome:` block

Supported, documented in the config file itself, survives a `git pull`. Anything
here you can change with confidence.

### 2.2 Layer B — the shell block of a rule in `workflow/rules/shared/80_mobilome.smk`

The Python scripts expose more command-line flags than the rules currently pass.
These are real, tested knobs — the rule just uses the script's default. To change
one you add the flag to the rule's `shell:` block.

The flags in this category, all on `rule conjscan_ice` (line 829) unless noted:

`--window-bp`, `--min-element-bp`, `--min-ime-element-bp`, `--max-element-bp`,
`--boundary-bp`, `--integrase-window-bp`, `--att-flank-window-bp`;
plus `--min-length-bp` on `rule isescan_table` and `--min-chromosome-bp` on
`rule mobilome_replicons`.

Concretely, to raise the IME size floor you edit the shell block:

```python
        shell:
            """
            python {params.script} \
              --sample {wildcards.sample} \
              --conjscan-tsv {input.conjscan_dir}/best_solution.tsv \
              ...
              --replicons {input.replicons} \
              --min-ime-element-bp 8000 \
              {params.strict_boundary} \
              --out-table {output.table} \
              --out-audit {output.audit} > {log} 2>&1
            """
```

Nothing else changes — but two warnings, because both are easy to assume the
other way round.

**The script does not sanity-check these values.** `argparse` enforces only that
they are integers. `--min-ime-element-bp 80000` is accepted, silently drops every
IME-architecture call, and the result is indistinguishable from a genome that has
none. (`coverage_profile` is guarded, but that guard lives in the Snakefile, not
in the script — see [§8.6](#86-a-profile-coverage-value-outside-0-1-is-silently-catastrophic).)

**A threshold reaches the audit file only when it actually fires.** The value
appears inside the `detail` of the row it caused — `cluster_shorter_than_min`
names the floor and the measured span, `integrase_without_conjugation_machinery`
names the integrase window. There is no unconditional "settings used" row: on a
sample where nothing was dropped, the audit and the run log say nothing at all
about which floors were in force. If you change a Layer B knob, write the value
down yourself.

### 2.3 Layer C — module constants in the Python

Constants at the top of `colocalise.py`, `att_search.py` and `conjscan_to_ice.py`
that no flag reaches. Changing one is a source edit. `pytest workflow/scripts`
(365 mobilome tests, 420 in the full suite) is the check that you did not break
an invariant, but the tests encode today's thresholds in places, so expect to
update a test alongside a deliberate change.

### 2.4 If you change a Layer B knob twice, promote it to Layer A

That is exactly how `coverage_profile` got there, and the pattern is four small
steps:

1. read the value into a `MOBILOME_*` constant. Most of them live in
   `00_common.smk` as `_mobilome_cfg.get("your_key", <default>)`;
   `MOBILOME_COVERAGE_PROFILE` instead sits at the top of `80_mobilome.smk`
   (line 130) as `(config.get("mobilome") or {}).get("coverage_profile", 0.5)`,
   because only the two MacSyFinder rules read it and the whole ICEscan layer is
   deliberately kept readable in one file. Either home is fine — pick the one
   whose rules use the knob;
2. validate it immediately and `sys.exit()` with a plain-language message if it
   is nonsense — a silently wrong threshold looks exactly like a genome with
   nothing in it;
3. add it to the rule's `params:` and reference it in `shell:`;
4. document the default **and the measurement behind it** in
   `config/config.yaml`, next to the key.

Step 4 is the one that matters. Every threshold in this module is a convention
rather than a biological constant, and the config file is where that gets said.

---

## 3. The inventory

Defaults are as shipped. "Measured" means there is an artefact on disk behind
the claim; "not measured" means nobody has swept it, and you should treat any
change as an experiment with the audit file as your instrument.

### 3.1 Master switches

| Knob | Layer | Default | What it means |
|---|---|---|---|
| `mobilome.run` | A | `false` | The whole module. Turning it on downloads the CONJscan models, which are CC BY-NC-SA (non-commercial). |
| `mobilome.icescan.run` | A | `false` | Second MacSyFinder model set, unioned with CONJscan. Adds the IME and AICE *models*, which CONJScan 2.1.0 does not carry, plus 21 profiles. Also CC BY-NC-SA. Needs `url` or `dir` as well: `run: true` with neither exits at parse time rather than running without the layer. |
| `mobilome.tncentral.url` / `.dir` | A | empty | The curated transposon/integron naming layer. Empty means tier 4 is unreachable. |
| `mobilome.iceberg.urls` / `.dir` | A | empty | ICE naming. Adds no elements and changes no tier — it turns "predicted self-transmissible element" into a name you can look up. |
| `mobilome.isosdb.fasta_url` | A | empty | Read-mapped IS copy number. Illumina/hybrid only. Changes no tier. |

**Measured, ICEscan on vs off** (18 curated ICEs, 12 curated IMEs, 12 negative
genomes, coverage 0.5; both arms re-run from the same MacSyFinder searches so the
comparison is like for like):

| | ICEs | IMEs | calls on the negative set |
|---|---|---|---|
| ICEscan **off** (default) | 15/18 | 2/12 | 2 |
| ICEscan **on** | 15/18 | 5/12 | 5 |

So ICEscan buys three IMEs and costs three extra calls on genomes with no curated
element. All three extras are in *Halobacterium* sp. NRC-1 — an **archaeon**,
where the bacterial models have no business firing — and all three come out
`cime_or_island` / `passive` / `low` / `evidence_level=profile_hits_only`, i.e.
carrying no mobility claim and no tier. ICE recall is identical either way.

> **A stale 3/12 is in circulation** — it appears in `config/config.yaml`
> ("3 detected before, 5 after") and in older notes. It came from a run made
> before the current nesting and size-floor rules, whose third "detection" was a
> **181,279 bp** call laid over the **11,112 bp** IME_SsuNSUI002_NS. `score.py`
> counts any overlap, so it scored as a hit; it is the same artefact
> [§8.1](#81-lowering-a-threshold-changes-what-detected-means) describes. Re-run
> with today's code the honest figure is **2/12**.

**ICEscan does not create the `ime` class**, and it is worth being exact about
this because the config comment reads that way. BacFlux's own rules already call
an element an IME when the anchors are an integrase plus a relaxase with no
mating apparatus, and CONJscan's MOB models supply plenty of those relaxases.
Measured over the 28 benchmark genomes for which both arms were run: **10 IME
rows without ICEscan, 22 with**. So turning ICEscan off roughly halves the IME
rows rather than removing them. The `aice` class *is* ICEscan-only — both AICE
rows disappear without it.

### 3.2 Machinery detection

| Knob | Layer | Default | Meaning |
|---|---|---|---|
| `mobilome.coverage_profile` | A | `0.5` | Fraction of an HMM profile a protein must align to before the hit is kept. Applied to **both** MacSyFinder searches, deliberately: their hit tables are merged, so different stringencies would make an element's class depend on which model set was looser. |

Not identity, not an E-value. A hit can be statistically beyond doubt and still
be dropped by this rule alone.

**Lowering it** admits divergent full-length relaxases that share only a
catalytic core — and also genuine fragments. **Raising it** does the reverse.

**Measured** (40 genomes, both searches re-run at each value, caller re-run,
scored with the unmodified `score.py`):

| coverage | curated ICEs | curated IMEs | negative-set calls (32.6 Mb) |
|---|---|---|---|
| **0.5** | 15/18 | 5/12 | 5 union / 2 CONJscan-only |
| 0.4 | 15/18 * | 6/12 | 7 / 4 |
| 0.3 | 15/18 * | 6/12 | 9 / 5 |

> \* **`score.py` will print 16/18 at 0.4 and 0.3, and the table above says 15
> deliberately.** The extra "detection" is a 3,187 bp call clipping the tail of a
> 161 kb element — 1.9% of it — and classed `ime` rather than `ice`. It is
> counted because `MIN_OVERLAP_BP = 1`. See
> [§8.1](#81-lowering-a-threshold-changes-what-detected-means). Honest ICE recall
> is 15/18 at every threshold tested.

Across 30 curated elements, lowering recovers exactly **one**: a 23 kb IME in
*Faecalibacterium duncaniae*, found in full (right edge 53 bp off the published
one) and correctly classed. It appears at 0.4 and gains nothing further at 0.3.

The cost at 0.4 is two more calls, in *E. coli* K-12 MG1655 and *P. aeruginosa*
PAO1 — the two most thoroughly described genomes in bacteriology. Read the class
before panicking: every call added at 0.4 is `cime_or_island` / `passive` /
`low` / `profile_hits_only`. The tiering absorbs them, which is what it is for.
At 0.3 that stops holding: a 21.6 kb IME appears in *S. aureus* N315 at
**medium** confidence with a real mobility claim, unverified either way. (That
one is in the union arm; the CONJscan-only arm's extra call at 0.3 is a
`conjugative_region` in *Aquifex aeolicus*, `low` / `profile_hits_only`.)

**Why the gain is so small when the raw evidence moves a lot.** 0.4 adds 5% more
MacSyFinder hits and changes the hit table in 13 of 40 genomes; 0.3 adds 12% and
changes 19 of 40. Almost none of it reaches the element table, because the caller
needs an integrase within 50 kb of conjugation machinery before it seeds
anything. **Coverage sits upstream of a stronger constraint; co-localisation is
what binds.** See [§7.3](#73-the-one-that-lowering-coverage-does-not-fix) for the
worked case.

**When to lower it anyway.** Across all 395 curated IMEs, **73** carry a relaxase
hit that is clean on E-value (i-evalue ≤ 0.001) and fails only this rule — 60 of
them the `T4SS_MOBT` family at a median 0.32 coverage against a 366-position
profile, largely the *Streptococcus salivarius* clade. If your isolates sit in
that space, 0.4 is worth trying. Nothing measured supports 0.3 for single
isolates.

> **Corrected 2026-07-31: this figure was quoted as 68.** The 68 was a stale
> number in the docstring of the script that builds the table, not in the table.
> Re-running `ime_ceiling/refutation/r07_final_tables.py` unchanged gives
> `coverage_limited=73  no_clean_signal=43`, and 279 + 73 + 43 = 395. `68`
> appears in `config/config.yaml` and in earlier notes; prefer 73. The
> `T4SS_MOBT` breakdown is unaffected — it is 60 of 73.

### 3.3 Element geometry — how anchors become an element

All Layer B, all on `rule conjscan_ice`. These are conventions, not biology; the
measured span is always reported next to the decision.

| Flag | Default | Raising it | Lowering it |
|---|---|---|---|
| `--window-bp` | 15,000 | chains more machinery genes into one element; risks fusing two neighbouring systems | splits one operon into several candidates, each then too small to survive the size floor |
| `--min-element-bp` | 8,000 | fewer, larger ICE candidates | admits short machinery spans; audit reason `cluster_shorter_than_min` stops firing |
| `--min-ime-element-bp` | **2,000** | see below — this is the sharpest knob in the module | admits sub-2 kb clusters, which physically cannot hold both a relaxase and an integrase |
| `--max-element-bp` | 500,000 | admits runaway anchor chains; audit reason `cluster_longer_than_max` | drops genuinely large ICEs |
| `--boundary-bp` | 1,000 | more candidates flagged as probably truncated and capped at **low** | fewer caps, and you stop being told that a call ran into the end of its contig |
| `--integrase-window-bp` | 50,000 | an integrase further from the machinery may still anchor the element | more `integrase_without_conjugation_machinery` orphans, so fewer ICE/IME calls |

> **`--min-ime-element-bp` is an EMPIRICAL CUT, not a published bound.** There is
> no published minimum size for an IME's machinery. 2,000 bp was chosen against a
> handful of Phase 7 observations in one benchmark. Say so in any methods
> write-up rather than implying it is a property of IMEs.

Why the separate floor exists at all: the floor is applied to the **anchor-cluster
span**, and that span scales with the *number of machinery genes*, not with the
element's length. An ICE carries a ~20-gene T4SS operon, so its span is naturally
tens of kb. An IME carries a relaxase and an integrase — two genes — so its span
is 1–6 kb. An 8,000 bp floor therefore selects for ICEs by construction. AICEs
take the same lower floor for the same reason (three genes, not twenty).

**Measured**, over the 28 small-machinery calls (IME + AICE) in the benchmark:

| floor | calls dropped (of 28) | curated IMEs still detected (of 12) |
|---|---|---|
| **2,000 (default)** | 0 | 5 |
| 3,000 | 5 | 5 |
| 5,000 | 7 | 3 |
| 8,000 (= the ICE floor) | 15 | **1** |

The five detected curated IMEs have machinery spans of 3,840 / 4,959 / 5,942 /
6,252 / 8,164 bp. Raising the floor to the ICE value keeps only the last of them.

Note the shape of that table: the recall cliff is between 3,000 and 5,000, not at
the default. **3,000 bp drops five calls and costs no curated IME on this
evidence** — the smallest detected one has a 3,840 bp span. That does not make
3,000 the better default (five calls is a small sample, and none of the five it
drops is a known false positive), but if you want a tighter floor with a measured
cost of zero, that is where it is. It does *not* remove the one IME call on the
negative set, which is 5,633 bp and survives every floor below 8,000.

`--integrase-window-bp` is deliberately much larger than `--window-bp`:
conjugation machinery is an operon and clusters tightly, whereas the integrase
sits at the element boundary, tens of kb away on a large ICE. It was widened to
50 kb after a 15 kb window produced `conjugative_region` instead of ICE on two
clinical *K. pneumoniae* isolates — on both, the ICE*Kp* integrase sits 33,361 bp
from the machinery cluster — and put the ICE label on a different element 1.1 Mb
away. icefinder-opt's README records the same fix upstream.

### 3.4 Boundaries — the *att*-site search

| Knob | Layer | Default | Meaning |
|---|---|---|---|
| `mobilome.require_trna_boundary_for_high` | A | `false` | Whether an element must have tRNA-anchored ends before it may be called **high** confidence. |
| `--att-flank-window-bp` | B | 50,000 | How far beyond the machinery span to look for the element's real ends. |
| `MIN_ATT_REPEAT_BP` | C (`att_search.py`) | 15 | Shortest exact repeat that may be an *att* site. Matches ICEfinder2 and icefinder-opt; DBSCAN-SWA uses 12. |
| `DENOVO_MAX_EXPECTED_CHANCE_MATCHES` | C | 0.05 | How many chance repeats we tolerate when picking the shortest believable length for *this* pair of flanks. Two 50 kb flanks share a 15 bp repeat about 2.3 times by luck, so a fixed 15 bp floor would report noise; the effective floor resolves to 18 bp at both 30 kb and 50 kb. |
| `MAX_ATT_COPIES_OUTSIDE_TRNA` | C | 2 | Copies of the repeat allowed outside any annotated tRNA before it is judged a repeat *family* (rRNA, REP/BIME, paralogues) rather than a single integration scar. |
| `TRNA_ANCHOR_BONUS_BP` | C | 10 | Score bonus for a repeat with a copy inside an annotated tRNA. Modest on purpose: it breaks ties, it does not let a marginal 15 bp repeat outrank a convincing 40 bp one. |
| `MAX_SEED_MATCHES` | C | 200,000 | Compute guard in repeat-dense regions. Reaching it is reported, never silent. |

`require_trna_boundary_for_high` folds two different questions into one number:
*is this an ICE?* (anchors, machinery, one contig) and *where does it end?*
(`boundary_method`). By default they are reported side by side and only the first
sets confidence.

**Measured cost of turning it on**: across the 85 element rows produced by the 52
genomes and element deposits in the benchmark tree — all complete GenBank records,
none of them drafts — 30 reach high confidence, and only **7** of those 30 have
`boundary_method=tRNA`. Turning the flag on demotes 23 of 30 high calls to medium
**even on closed genomes**. On a fragmented
draft, where the flanks are usually simply absent, it caps essentially
everything: confidence stops discriminating and becomes a constant.

Use it on closed long-read assemblies, where an unresolved boundary genuinely is
a warning sign rather than the norm. Either way, cargo is never assigned from an
unresolved boundary: a de novo repeat is reported but never widens an element, so
the interval stays the machinery span — a floor, never an invention.

`--att-flank-window-bp` was 30,000 until commit `4a93d89`. It was widened because
SPI-7's *attR* sits 5,800 bp outside a 30 kb window, so the element was reported
50 kb short for want of somewhere to look. 80 kb, 120 kb and 200 kb change no
element's answer on the benchmark, so 50 kb is where the curve flattens rather
than an arbitrary larger number. The cost of a wider window is more candidate
repeats to rank.

### 3.5 The IS layer

| Knob | Layer | Default | Meaning |
|---|---|---|---|
| `mobilome.contig_boundary_bp` | A | 100 | How close to a contig end counts as "at the boundary" when **flagging** an IS call. Feeds `at_contig_boundary` and `fraction_at_contig_boundary` in the IS summary. |
| `--min-length-bp` | B (`rule isescan_table`) | **0** | Drop IS shorter than this, with a reason. 0 means partials are kept and tiered rather than discarded. |

ISEScan itself is run with **no** `--removeShortIS`, on purpose: we keep partials
and tier them ourselves (spec §6). Note that spec §9's "discard MGEs <500 bp"
convention is **not** applied to IS calls — set `--min-length-bp 500` if you want
it.

`contig_boundary_bp` is deliberately **not** the same number as the window that
caps confidence in the mobility table. That one is a fixed 1,000 bp
(`CONTIG_END_WINDOW_BP` in `colocalise.py`) and answers a different question: was
there enough flanking sequence to have *seen* a neighbouring element at all? An
IS 500 bp from a contig end is not "at the edge", but you still could not have
seen its partner 2 kb away, so absence of evidence there is not evidence of
absence.

A high `fraction_at_contig_boundary` means the assembly broke exactly where the
IS elements are, and the located IS count should be read as a floor. For
reference, the two hybrid clinical assemblies we have run score 0.0089 (1 of 112
IS calls) and 0.0103 (1 of 97) — near-closed genomes. A short-read draft will be
far higher.

### 3.6 AMR × MGE co-localisation

| Knob | Layer | Default | Meaning |
|---|---|---|---|
| `mobilome.max_composite_span_bp` | A | 20,000 | Longest span accepted for a composite transposon (tier 3). Real composites run from ~2.5 kb (IS26 translocatable units) to ~25 kb. |
| `HYBRID_PROMOTER_MAX_BP` | C | 500 | How close an upstream, correctly oriented IS must be before we say it may be driving the gene from an outward-reading promoter (tier 2). ISEcp1 sits 42–266 bp upstream of blaCTX-M. |
| `IS_ADJACENT_MAX_BP` | C | 5,000 | Beyond this the nearest IS is still reported in `distance_bp` but no longer sets `mge_context`. |
| `CONTIG_END_WINDOW_BP` | C | 1,000 | Triggers the boundary flag and the low-confidence cap on an AMR row. |
| `MIN_AMR_FRACTION_INSIDE_ELEMENT` | C | 0.9 | Fraction of the AMR gene inside a named element before we say the element carries it. Borrowed from the EBI pipeline's published thresholds. |
| `MIN_AMR_FRACTION_OVERLAPPED_FOR_DISRUPTION` | C | 0.5 | Fraction overlapped by an IS before we call the CDS disrupted (likely inactivated) rather than merely overlapping. |
| `ABUTTING_IS_MAX_GAP_BP` | C | 25 | How flush an IS must be against a **partial** AMR hit to read the pair as "the IS split this gene". Tight on purpose: it tests "butted up against", not "nearby". |
| `MIN_COVERAGE_FOR_INTACT_GENE` | C | 90.0 | Below this % coverage an AMRFinderPlus hit is treated as a fragment. Only used with the abutting-IS test. |
| `ORIENTATION_EXEMPT_IS_FAMILIES` | C | `{IS6, IS26, …}` | Families exempt from the same-orientation rule for composites. **Do not remove IS26 from this set** — it forms translocatable units with its copies in *direct* orientation, and without the exemption the single most clinically important AMR architecture is the one you miss. |
| `UNINFORMATIVE_FAMILY_VALUES` | C | `{"", NA, NEW, UNKNOWN, ISNCY, …}` | ISEScan family labels that carry no information about *which* element this is. Two elements both labelled "new" must not be called a composite on family grounds. |

None of the Layer C values here have been swept. They are documented so you can
see what a given audit reason is testing, not because we have evidence that a
different value is better.

### 3.7 The naming layers

| Knob | Layer | Default | Meaning |
|---|---|---|---|
| `mobilome.tncentral.min_identity` | A | 90.0 | Percent identity before a BLAST hit may confer a curated **name**. Below it the sequence may be a relative of the element, but the name would claim more than the data shows. |
| `mobilome.tncentral.min_reference_coverage` | A | 0.8 | Fraction of the **reference element** that must be present. Measured against the reference, not your contig: "is the whole of this known transposon here?" |
| `mobilome.iceberg.min_identity` | A | 80.0 | Identity floor for an ICE name. |
| `mobilome.iceberg.min_overlap_fraction` | A | 0.5 | How much of *our* candidate the curated element must cover to be treated as the same thing. Lenient next to the transposon cascade because ICEs are mosaic and their cargo varies between strains. Below 80% of the reference the name is suffixed `-like`. |

`min_reference_coverage` is the one with a measured story, and it is the reason
tier 4 is so hard to reach on fragmented assemblies. On the clinical hybrid run
`TUM24772`, every Tn*Ecp1.1* candidate was refused with:

```
TnEcp1.1: only 49% of the 3417 bp reference element is present (needs 80%).
A fragment of a transposon is not that transposon - on a short-read assembly
this usually means the element is split across contigs.
```

and a second at 12%. That refusal is correct behaviour — a fragment of a
transposon is not that transposon — but it means the gene stays at tier 2 and you
have to look in `{sample}_named_elements_discarded.tsv` before concluding that a
genome has no named transposon in it. Lowering the threshold to admit a 49%
fragment has **not been measured**, and it would mean asserting a named
architecture from half the evidence; if you do it, read the audit and say in the
methods what coverage you accepted.

### 3.8 IS copy number (Illumina/hybrid only)

| Knob | Layer | Default | Meaning |
|---|---|---|---|
| `mobilome.isosdb.min_covered_percent` | A | 90.0 | A database entry must be covered end to end before its depth is believed. A partially covered entry is usually a conserved domain shared with another family, and averaging it in inflates every estimate. |
| `mobilome.isosdb.min_copy_number` | A | 0.5 | Below this multiple of the genome baseline the element is treated as absent. Under 1.0 on purpose: a real single-copy IS sits near 1×, and sampling noise plus mapping loss routinely pushes it to 0.6–0.8×. |

This leg changes no AMR gene's tier and must not: it says nothing about *where*
the extra copies are, only that they exist. It is the honest companion to the
"the located IS count is a floor" warning — it puts a number on how much the
assembler collapsed.

### 3.9 Replicon calling

| Knob | Layer | Default | Meaning |
|---|---|---|---|
| `--min-chromosome-bp` | B (`rule mobilome_replicons`) | 2,000,000 | Above this length an unclassified contig may be called the chromosome. Not Platon's 500 kb, because megaplasmids run to 1–2 Mb; between the two we say `unknown` and let the confidence cap reflect that. |

This matters more than it looks. An ICE integrates into a **chromosome** by
definition, so machinery on a contig called `plasmid` is a conjugative plasmid
region, not an ICE — audit reason `ice_demoted_on_plasmid_replicon`. On an
`unknown` contig you get `ice_on_unclassified_replicon` instead.

---

## 4. The two directions

### 4.1 More conservative

**For:** a regulatory dossier.

> ### Read this before you turn anything: no knob tightens tier 6
>
> The obvious framing — "a false *predicted self-transmissible* is the expensive
> error, so tighten the knobs" — does not survive contact with the code. Tier 6
> is set by exactly two things: an AMR gene inside an element classed `ice`, or a
> gene on a plasmid Platon typed `conjugative`. **Nothing in Layer A or Layer B
> touches either.**
>
> Measured, on the 40 benchmark genomes at coverage 0.5: turning ICEscan off
> leaves the ICE calls **identical** — 37 rows in both arms, the same
> 16 high / 19 medium / 2 low, and not one row's confidence changes. Seven calls
> shift their start coordinate by a few kb because ICEscan contributed an extra
> integrase anchor; the elements, classes and confidences are the same.
> `--min-ime-element-bp` acts only on IME/AICE-architecture clusters, and the
> naming thresholds change names, not tiers.
>
> So every lever below tightens the **tier-5** branch (IME, "mobilisable, needs a
> helper") or the naming layer. That is a real thing to want — it is just not the
> thing the heading promises. **Reporting less is not the same as claiming less.**
>
> There is a second reason this matters. Across **every** run kept on disk —
> the two clinical hybrids and the BAA-2146 positive control — every tier-6 row
> arrived by the *plasmid* route (`mge_context=plasmid`, a replicon Platon typed
> conjugative). **Not one came from `mge_context=ice`.** So the branch these
> knobs could conceivably narrow has not fired yet in practice, while the branch
> that does fire is decided by Platon and is outside the mobilome block entirely.
>
> What actually guards tier 6 is not tunable: three of four anchor classes for
> `high`, the absolute `spans_contigs` cap, the truncation check, and the
> chromosome-vs-plasmid replicon rule. If you want a stricter dossier, **filter
> the deliverable** rather than turn a knob:
>
> ```bash
> # Tier-6 rows you are willing to defend. Columns by name, not by number
> # (see section 1.1). Note this table writes booleans as yes / no / NA,
> # not TRUE / FALSE.
> awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) col[$i] = i; print; next}
>             $col["mobility_tier"] == 6 &&
>             $col["confidence"]    == "high" &&
>             $col["spans_contigs"] != "yes"' \
>     08.mobilome/S1/S1_amr_mobility.tsv
> ```
>
> `machinery_intact` is the column that says which kind of tier 6 you have.
> `yes` means CONJscan typed a complete conjugative system on that replicon; `NA`
> means the tier rests on Platon's counts alone and the row is capped at medium
> for exactly that reason, so the confidence test above already excludes it.
> Report the rows you dropped, with their reasons, alongside the ones you kept.
> That is a defensible position; a narrower config is not.

```yaml
mobilome:
  run: true

  # Keep MacSyFinder's own default. Nothing measured supports going lower for a
  # single isolate, and the two extra negative-control calls at 0.4 are the
  # reason.
  coverage_profile: 0.5

  # Off. Costs three IMEs out of twelve and no ICEs at all; removes three
  # profile-only calls in an archaeon.
  icescan:
    run: false

  # ON ONLY IF YOUR ASSEMBLIES ARE CLOSED. See the note below.
  require_trna_boundary_for_high: false

  # Naming layers: a name is a claim, so demand more of it. Neither of these
  # changes a tier - they change whether an element gets called by name.
  tncentral:
    url: "<your fetch URL>"
    min_identity: 95.0            # not measured; 90.0 is the shipped value
    min_reference_coverage: 0.9   # not measured; 0.8 is the shipped value
  iceberg:
    urls: ["<your fetch URLs>"]
    min_identity: 90.0            # not measured; 80.0 is the shipped value
    min_overlap_fraction: 0.5
```

These snippets are **edits to the shipped `mobilome:` block**, not replacements
for it — leave the keys you are not changing where they are.

Plus one Layer B edit, in `rule conjscan_ice`:

```
              --min-ime-element-bp 8000 \
```

**What it costs, measured.** Raising the IME floor to the ICE floor takes curated
IME recall from **5/12 to 1/12** and removes 15 of the 28 small-machinery calls
in the benchmark. What it buys is small and real: the one IME call on the
negative set — a 5,633 bp element in *S.* Typhimurium LT2 — disappears, taking
CONJscan-only negative calls from 2 to 1. Turning ICEscan off halves the IME rows
(22 → 10 over the 28 genomes run both ways) and removes both AICE rows; ICE
recall, and every individual ICE row, stay exactly as they were.

If your reason for raising the floor is that 2,000 bp feels unevidenced rather
than that you want fewer IME rows, **3,000 is the value with a measured cost of
zero** — see [§3.3](#33-element-geometry--how-anchors-become-an-element). 8,000
is a deliberate choice to stop reporting the IME architecture, not a tightening
of it.

> **Do not turn on `require_trna_boundary_for_high` on draft assemblies.** On the
> 52 *complete* records in the benchmark tree it already demotes 23 of 30
> high-confidence calls to medium, because only 7 of them have tRNA-anchored
> ends. On a short-read draft it caps essentially everything, and a confidence
> column that is always "medium" tells a reviewer nothing. Conservative does not
> mean "cap everything"; it means "make the strong claim rarer and better
> evidenced".

**What the audit will look like.** More `dropped` rows in
`{sample}_ice_discarded.tsv`, dominated by `cluster_shorter_than_min` with the
span and the threshold in `detail`; more `reference_coverage_below_threshold` and
`identity_below_naming_threshold` in `{sample}_named_elements_discarded.tsv`
(where the action word is `discarded`, not `dropped` — see
[§1](#1-start-at-the-audit-file-not-at-the-config)). Fewer `kept_flagged` rows
overall, because fewer candidates survive to be flagged.

**What the output then means.** Small IME-architecture clusters are not reported
at all, so **absence of an IME row is not evidence of absence** — the audit is
where they went. Say that in the dossier, alongside `fraction_at_contig_boundary`
from the IS summary. And say plainly that the ICE/tier-6 calls are the *same*
calls the default configuration makes: the narrowing applies to tier 5, not to
the self-transmissible claim.

### 4.2 Wider net

**For:** an exploratory survey, or an organism nobody has curated, where missing
a real element is the expensive error.

```yaml
mobilome:
  run: true

  # 0.4, not 0.3. 0.4 recovers one curated IME (Faecalibacterium duncaniae,
  # right edge 53 bp off the published coordinates) and every call it adds to
  # the negative set is passive / low / profile_hits_only. 0.3 adds nothing
  # further on the curated sets and produces the first negative-control call
  # that actually asserts mobility.
  coverage_profile: 0.4

  # On. Roughly doubles the IME rows and is the only way to get an AICE call at
  # all. `url` or `dir` is required alongside `run: true` - without one the
  # workflow exits at parse time rather than quietly running without the layer.
  icescan:
    run: true
    url: "https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/icefinder2lite/icf2_dbs.tar.gz"
    sha256: ""      # set it: the URL is unversioned

  require_trna_boundary_for_high: false

  # Wider composite window. A convention, not biology, and the measured span is
  # always in flanking_span_bp next to the call.
  max_composite_span_bp: 25000

  tncentral:
    url: "<your fetch URL>"
    min_identity: 90.0
    min_reference_coverage: 0.8   # leave it: below this you are naming fragments
```

Optional Layer B edits, in `rule conjscan_ice` — each of these is an
**experiment**, not a recommendation; none has been swept:

```
              --min-element-bp 5000 \
              --window-bp 20000 \
```

**What it buys, measured**, against the shipped default (ICEscan off, coverage
0.5) rather than against a half-changed configuration: IMEs **2/12 → 6/12**, ICEs
unchanged at 15/18, negative-set calls **2 → 7**. Both ends of each pair are
like-for-like — the "before" is CONJscan-only at 0.5, the "after" is the union at
0.4.

**What `--min-element-bp` would let in, for calibration.** At the default 8,000,
*P. aeruginosa* PAO1 drops two clusters at **267 bp** and **1,266 bp**, and
*S. aureus* N315 one at **1,263 bp** — all `cluster_shorter_than_min`, all on
genomes with no curated element. A floor of 5,000 does not reach any of them, so
that change is close to free on this evidence; a floor low enough to admit them
is a floor low enough to admit a lone machinery gene. There is no measured value
between 5,000 and 8,000, which is why the two Layer B lines above are labelled an
experiment.

**What the audit will look like.** Far fewer `dropped` rows and many more
`kept_flagged` ones — `no_assembled_system`, `few_anchor_classes`,
`mpf_marker_without_typed_system`, `machinery_not_intact`. That is the shape of a
wider net working correctly: the extra material arrives already labelled weak.

**What the output then means.** The table is now a **candidate list to
triage**, not a set of findings. Filter it before reading:

```bash
# Everything that actually claims mobility, with a system behind it.
# Columns by name, not by number (see section 1.1).
awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) col[$i] = i; print; next}
            $col["evidence_level"] != "profile_hits_only" &&
            $col["mge_class"]      != "cime_or_island"' \
    08.mobilome/S1/S1_ice_candidates.tsv
```

Nothing at `evidence_level=profile_hits_only` should be quoted as a finding
without looking at the genes yourself.

---

## 5. Symptom → knob

Every `reason` below is a literal string from the code — grep for it.

### "An element I know is there was not reported"

| Grep `{sample}_ice_discarded.tsv` for | What it means | Knob |
|---|---|---|
| `cluster_shorter_than_min` | The machinery span was below the size floor. `detail` gives the measured span and which floor applied (`--min-element-bp` or `--min-ime-element-bp`). | `--min-ime-element-bp` for IME/AICE architecture; `--min-element-bp` otherwise |
| `cluster_longer_than_max` | Anchors chained across the contig into one implausible element. | `--max-element-bp`, or **lower** `--window-bp` so the chain breaks |
| `no_conjugation_anchor` | Integrases with no relaxase and no T4SS nearby. A genome's ordinary complement of recombinases. | Not a threshold problem — see [§6](#6-what-is-not-tunable-and-why) |
| `integrase_without_conjugation_machinery` | The integrase marking your element was not attached to any machinery cluster. | `--integrase-window-bp` if the machinery is genuinely far away; but read [§7.3](#73-the-one-that-lowering-coverage-does-not-fix) first, because a *closer* integrase winning the cluster is the commoner cause and no knob fixes it |
| `untrusted_integrase_profile` | An ICEscan integrase profile we deliberately ignore (integron and housekeeping recombinases that would anchor elements onto the wrong gene). | Layer C: `ICESCAN_TRUSTED_INTEGRASE_PROFILES` in `conjscan_to_ice.py` |
| `integrase_match_is_a_transposase` | Bakta's product text matched the integrase regex but the protein is a transposase. | Layer C: `TRANSPOSASE_PRODUCT_PATTERN` |
| `conjscan_found_no_systems` | MacSyFinder assembled nothing at all. | `coverage_profile`, then `icescan.run` |
| `hit_id_not_in_annotation` | A MacSyFinder hit could not be joined to the Bakta GFF3. Not a threshold — a mismatched input pair. | Check that the Bakta run and the `.faa` MacSyFinder read are the same run |

And in `{sample}_named_elements_discarded.tsv`:
`reference_coverage_below_threshold` → `tncentral.min_reference_coverage`;
`identity_below_naming_threshold` → `tncentral.min_identity`. A third,
`tncentral_hit_is_a_plain_is`, is not a threshold at all: an IS carries only what
it needs to move, so it cannot hold a passenger AMR gene, and ISEScan already
inventories it.

If the element is there but no **AMR gene** picked it up, look in
`{sample}_amr_mobility_audit.tsv` for `near_but_outside_element_machinery_span`
(the element's reported interval is the machinery span, and your gene sits in the
cargo beyond it — a boundary problem, not a threshold one) or
`nearest_mobile_element_too_far` (`IS_ADJACENT_MAX_BP`, Layer C).

### "The interval looks far too big / far too small"

| You see | Read | Then |
|---|---|---|
| interval much **smaller** than the published element | `boundary_method` in `{sample}_ice_candidates.tsv` | If `none`, the interval **is** the machinery span — a floor, by design. Compare `start`/`end` with `machinery_start`/`machinery_end`. Try `--att-flank-window-bp`, but note 80/120/200 kb changed no answer on the benchmark |
| interval much **larger** | `n_anchors`, `anchor_classes`, and `evidence_recorded / system_merge_refused_too_long` | Lower `--window-bp` (fuses less) or `--max-element-bp` (drops the chain) |
| a de novo repeat was found but the interval did not move | `kept_flagged / denovo_att_reported_not_applied` | Working as intended: ~16% of arbitrary spans on a real chromosome throw up a de novo repeat, so it is reported and never applied. No knob — this is a design decision, not a threshold |
| ends resolved but at the wrong place | `att_length_bp`, `att_trna`, `attL`, `attR` | Layer C: `MIN_ATT_REPEAT_BP`, `MAX_ATT_COPIES_OUTSIDE_TRNA` |

### "Everything comes out at medium confidence"

First check which cap fired — the audit names one per rule:

```bash
awk -F'\t' '$5=="kept_flagged"' 08.mobilome/S1/S1_ice_discarded.tsv | cut -f6 | sort | uniq -c
```

| Cap reason | Meaning | Knob |
|---|---|---|
| `few_anchor_classes` | fewer than 3 of {integrase, relaxase, coupling protein, T4SS} | `MIN_ANCHOR_CLASSES_FOR_HIGH` (Layer C) — but read the next box before you touch it |
| `machinery_not_intact` | low system wholeness, a decayed model, or a truncated relaxase/VirB4 | `WHOLENESS_INTACT_MIN` / `PROFILE_COVERAGE_INTACT_MIN` (Layer C, both 0.7) |
| `at_contig_boundary` | machinery runs into a contig end → capped **low** | `--boundary-bp` |
| `spans_contigs` | hits on more than one contig → capped **low**, absolutely | none, deliberately |
| `no_assembled_system` | profile hits only, no system → capped **low** | `coverage_profile`, or accept the tier |
| `boundary_not_resolved` / `boundary_denovo_only` | only fires when `require_trna_boundary_for_high` is on | turn that flag off |
| `contig_length_unknown` | no length for the contig, so we cannot say whether the element runs off the end | fix the `--contig-lengths` input |

> ### An IME can never reach high confidence. That is the design, not a knob.
>
> High confidence requires at least **3 of the 4 anchor classes**
> (`MIN_ANCHOR_CLASSES_FOR_HIGH = 3`). An IME is *defined* by having an integrase
> and a relaxase and **no mating-pair apparatus** — that is what makes it
> mobilisable-with-a-helper rather than self-transmissible. So it starts with
> two, and a coupling protein is the only way it ever reaches three.
>
> **Measured:** of the 26 IME rows across the whole benchmark, **none** is high.
> 25 have exactly 2 anchor classes; 1 has 3 and is still medium. 18 of the 26 are
> additionally flagged `machinery_degraded` for `low_system_wholeness` — an
> artefact of CONJscan's MOB model describing more genes than an IME carries, so
> finding the relaxase alone scores 0.333 wholeness. By contrast 30 of the 53
> ICE rows are high.
>
> Lowering `MIN_ANCHOR_CLASSES_FOR_HIGH` to 2 would not even lift the ceiling
> cleanly, and it is worth knowing why before reaching for it. Of the 26 IME rows,
> **only 8 would become eligible**: 17 of the rest are independently capped at
> medium by `machinery_not_intact`, and 1 by `at_contig_boundary`. So the change
> buys 8 relabelled rows — while simultaneously making every two-anchor ICE-ish
> cluster eligible for `high`. It does not make the evidence better; it makes the
> label stop meaning what it says. **The honest read is: for an IME, medium is the
> ceiling, and `mobility` plus `machinery_intact` are the columns that carry the
> information.**

### "I get calls on a genome that should have none"

Read the **class** before the count. On our 12-genome negative set (32.6 Mb, no
curated element), coverage 0.5, ICEscan on: 5 calls. Three are in *Halobacterium*
sp. NRC-1, an archaeon, and every one is `cime_or_island` / `passive` /
`profile_hits_only` / `low`, with `missing_components = relaxase, coupling
protein`. They assert nothing. The tiering absorbed them exactly as intended.

| You see | Do |
|---|---|
| `evidence_level=profile_hits_only` | This is the module's weakest evidence: MacSyFinder saw the profiles but never assembled a system. Filter these out for reporting, or raise `coverage_profile` |
| `mge_class=cime_or_island`, `mobility=passive` | No mobility claim and no tier. Not a false positive in the sense that matters |
| a call at **medium** or better with a real mobility claim | Investigate before dismissing — see [§8.3](#83-a-negative-control-must-be-inspected-not-counted) |
| calls in an archaeon or a very distant lineage | CONJscan and ICEscan models are bacterial. Treat any call as suspect on those grounds alone |

---

## 6. What is not tunable, and why

These are limits of the design and of the data. No parameter moves them, and a
methods section should say so.

**1. The 70.6% IME detection ceiling.** Measured over all 395 curated IMEs at
production settings: 279 carry a relaxase we can detect at all. By size, 277/353
(78.5%) at ≥5 kb and 231/262 (88.2%) at ≥8 kb. So the pilot's 5/12 is roughly 5
of 8–9 *achievable*, not 5 of 12. The remaining 30% are not waiting behind a
threshold — validated by a shuffled-protein control: at production settings 279
of 395 real versus **0 of 395** shuffled.

Size-matched, IMEs are detected as well as or better than ICEs in every bin
(5–10 kb: 62.9% vs 0%; 10–20 kb: 87.6% vs 47.3%; ≥80 kb: 100% vs 95.8%). The IME
branch is not structurally disadvantaged; the earlier appearance of that was a
size confound.

**2. Elements that carry no relaxase at all.** 16 curated IMEs at ≥5 kb have no
clean signal, and they are biologically coherent: NBU1, NBU2, Tn*4555* (*trans*-
mobilised *Bacteroides* units) and the MGI series. They do not encode a relaxase
because **they do not use one** — they are mobilised in *trans* by a co-resident
element. This module is relaxase-anchored, so they are **out of scope by design,
not misses**. No threshold reaches them; detecting them needs a different anchor.

**3. Elements below ~2 kb.** 22 curated IMEs (5.6%) are under 2 kb — 17 of them a
single ORF, the rest two or three — too small to hold both a relaxase and an
integrase. That is what ICEberg deposited, not a detection failure: 18 of the 22
produce no relaxase signal at any threshold. Lowering `--min-ime-element-bp`
below 2,000 does not recover them; it admits noise.

**4. `spans_contigs` caps at low, absolutely.** MacSyFinder is run over the whole
draft as one pseudo-replicon, so genes it joined across a contig break may simply
be unrelated genes that happen to be adjacent in the file. This is not a
threshold and is not exposed.

**5. Boundaries on fragmented assemblies.** If the flanking sequence is not in
the contig, no parameter finds an *att* site. `boundary_method=none` is the
correct answer, and the interval is honestly a floor.

**6. We have never measured this module on a fragmented draft assembly — which
is BacFlux's main input.** Every validation figure on this page comes from
**complete reference genomes**. The only non-reference genomes ever run through
the module are two *hybrid* clinical assemblies of **5 and 6 contigs**. Nothing
with hundreds of contigs has been scored against a truth set.

Two consequences you should assume until measured: `at_contig_boundary` and
`spans_contigs` will fire far more often, dragging confidence down for real
reasons; and tier 4 will be even harder to reach, because a curated element split
across contigs fails `min_reference_coverage` (that is already visible at 12–49%
on 5- and 6-contig assemblies). Read `fraction_at_contig_boundary` in
`{sample}_is_summary.tsv` first, every time.

---

## 7. Worked examples

### 7.1 What a good call looks like — SPI-7 in *Salmonella* Typhi CT18

Config: everything at defaults, `mobilome.run: true`, ICEscan on. No tuning.

From `{sample}_ice_candidates.tsv` (`work/AL513382/ice_elements.tsv`):

```
mge_id              AL513382.1|ice-4423868:4507270
mge_class           ice
mobility            predicted self-transmissible
start / end         4409517 / 4543098        (133,582 bp)
machinery_start/end 4423868 / 4507270
n_anchor_classes    4      integrase,relaxase,t4cp,t4ss
relaxase_type       MOBM        mpf_type  G        mpf_typed_system  TRUE
sys_wholeness_min   0.955       machinery_intact  TRUE
evidence_level      system      missing_components  none
boundary_method     tRNA
attL / attR         4409517..4409543 / 4543072..4543098
att_sequence        GTGGTGCCCGGACTCGGAATCGAACCA   (27 bp, 0 mismatches)
att_trna            tRNA-Phe(gaa)
confidence          high
```

Curated coordinates are 4,409,574–4,543,073 (133,500 bp). The call is 57 bp early
at the left end and 25 bp late at the right.

Audit line that matters:

```
boundaries_resolved   att_pair_found_tRNA
```

Note what the *att* pair did to the interval: the machinery span is
4,423,868–4,507,270 (83,403 bp) and the reported element is
4,409,517–4,543,098 (133,582 bp) — **50,179 bp wider**, 14,351 to the left and
35,828 to the right. That is the whole reason
boundaries are worth resolving: `colocalise.py` intersects `start`/`end` against
the AMR genes to decide what is **cargo**, and an ICE is usually much larger than
its *tra* cluster. Leaving the machinery span in place would systematically
understate what travels with the element — the dangerous direction for an AMR
report.

This element also needed `--att-flank-window-bp` at 50,000: its *attR* sits
5,800 bp outside a 30 kb window, and at the old default SPI-7 was reported 50 kb
short.

### 7.2 A knob that changes the answer — Tn*4451* and the IME size floor

*Clostridium perfringens* U15027, the curated IME Tn*4451* at 277–6,614.

Default config (`--min-ime-element-bp 2000`), from `ice_elements.tsv`:

```
mge_id              U15027.1|ime-352:6293      (5,942 bp; curated 6,338 bp)
mge_class           ime
mobility            mobilisable (needs a helper) - machinery incomplete
n_anchor_classes    2      integrase,relaxase
relaxase_type       MOBV   integrase_products  Recombinase
conjscan_models     CONJScan/Chromosome/MOB, ICEscan/Chromosome/IME
sys_wholeness_min   0.333  degraded_reason  low_system_wholeness
missing_components  coupling protein, mating-pair apparatus
boundary_method     none
contig_length       6730   dist_to_contig_start  351
at_contig_boundary  TRUE
confidence          low
```

Now set `--min-ime-element-bp 8000` (the conservative recipe). The element's
machinery span is 5,942 bp, so it fails the floor, disappears from the table, and
the audit gains:

```
dropped   cluster_shorter_than_min
          machinery span is 5942 bp, below the --min-ime-element-bp threshold
          of 8000 bp. Reported as too small to be a credible integrative
          element; the hits themselves are still in CONJscan's own output.
```

Three things to read out of this one example:

- The default is doing real work. Four of the five detected curated IMEs have
  machinery spans between 3,840 and 6,252 bp; the 8,000 bp ICE floor keeps only
  one of the five.
- The call is already honest about its weakness without any tuning:
  `confidence=low`, because `at_contig_boundary` is TRUE — this deposit is a
  6,730 bp fragment and the element starts 351 bp from its end.
- `sys_wholeness_min = 0.333` is **not** evidence that Tn*4451* is decayed. It is
  CONJscan's MOB model describing three genes where an IME carries one relaxase.
  The `machinery incomplete` wording is conservative-by-construction here, and
  that is the same artefact that keeps 18 of 26 IME rows out of high confidence.

### 7.3 The one that lowering coverage does *not* fix

*Streptococcus salivarius* CP002888.1, curated IME_Ssal57I_tRNAlys at
70,384–75,506. This clade is precisely the `T4SS_MOBT` biology that the IME
ceiling study says is held back by profile coverage — 60 of the 73 borderline
elements. So lowering `coverage_profile` ought to recover it.

**It does not.** Measured at 0.5, 0.4 and 0.3: still missed at every value. The
audit says why in two lines.

```
evidence_recorded  integrase_attached_beyond_cluster_window
    CP002888.1:47906-48934: an integrase at 49258-50742 (Integrase, from
    bakta_product) sits 324 bp away - beyond the machinery clustering window
    but within the 50000 bp integrase window.

not_applicable     integrase_without_conjugation_machinery
    5 integrase(s) were found with no conjugation machinery within 50000 bp,
    so they anchor no element. Examples: CP002888.1:70481-71620; ...
```

The relaxase hits were being found at 0.5 all along. The machinery cluster at
47,906–48,934 took the integrase 324 bp away, and the integrase at 70,481–71,620
— the one that actually marks the curated element — was left with nothing to
attach to. The reported call, `CP002888.1|ime-47906:50742`, does not overlap the
curated interval at all.

**Which knob?** Not `coverage_profile`. Not `--integrase-window-bp` either: the
element's integrase was inside the 50 kb window; it lost to a nearer competitor,
because each integrase is attached to at most one cluster. This is a
co-localisation limitation, not a threshold, and it is the concrete form of the
general finding that **coverage sits upstream of a stronger constraint**.

The general lesson: *read the audit before you turn the knob the literature
points at.* The ceiling study's 73-element figure is a real upper bound on
**hits**; it does not translate 1:1 into **elements**.

### 7.4 What a tier looks like when a naming layer is missing

*K. pneumoniae* ATCC BAA-2146, `blaCTX-M-15`, from
`docs/validation/BAA-2146_amr_mobility.tsv` — a run with **no** TnCentral layer:

```
contig NZ_CP006659.2  replicon chromosome
mge_context is_adjacent   distance_bp 48   orientation same
is_family IS1380         n_flanking_is 2
mobility_tier 2          expression_modulation_not_mobilisation   confidence medium
```

IS*Ecp1* (family IS1380) sits 48 bp upstream, correctly oriented — inside
`HYBRID_PROMOTER_MAX_BP` (500), so tier 2. The same gene's plasmid copy on
`NZ_CP006662.2` scores tier 5 (`plasmid`, `plasmid_mobility=mobilisable`,
confidence high).

Tier 2 is *true* but it sells the situation short: IS*Ecp1* does not merely
modulate `blaCTX-M-15`, it **mobilises** it, and Tn*Ecp1.1* is that unit. Reaching
tier 4 needs `mobilome.tncentral.url` set **and** ≥80% of the reference element
present. On the clinical hybrid runs that did have the layer on, every
Tn*Ecp1.1* candidate was refused at 12% and 49% coverage (see §3.7). See
`mobilome_worked_example.md` for the full treatment and its correction notice —
**no run kept on disk contains a tier 4 row**, so treat tier 4 as a real code
path that is currently unexercised rather than as something a knob will produce.

---

## 8. Heads-up: the traps we hit

### 8.1 Lowering a threshold changes what "detected" means

Recall scores that count **any overlap** flatter a looser setting, because one
large call can swallow several curated elements and each counts as a hit.

Concretely: at `coverage_profile: 0.4` the ICE pilot appears to go 15/18 → 16/18.
It does not. The extra "detection" is ICE*Psy10*, matched by a **3,187 bp** call
clipping the tail of a **161 kb** element — 1.9% recovered — and classed `ime`
rather than `ice`. Our own `score.py` counts it because `MIN_OVERLAP_BP = 1`.
Honest ICE recall is 15/18 at every coverage value tested.

Always report `recovered_fraction` and the class alongside the count. When you
sweep a knob, check that a new "detection" is the element and not a fragment of
it.

(The IME gain at 0.4 **is** real by the same test: IME_FprA2-165_tRNAlys_2 is
called at 666,518–723,176, fully spanning the curated 699,930–723,123 with the
right edge 53 bp off, and correctly classed. It still over-extends 33 kb
leftward — the honest caveat that belongs next to it.)

### 8.2 The permissive HMM pass is about 70% noise

There is a tempting number in the ceiling study: at a permissive setting, 373 of
395 curated IMEs show a relaxase hit, against 279 at production settings. It
looks like a 94% ceiling.

The shuffled-protein control kills it. Run the same search over shuffled
proteins: at production settings **0 of 395** shuffled sequences produce a hit;
at the permissive setting **277 of 395** do. The permissive pass is roughly 70%
noise and **must never be quoted as a ceiling**. 279/395 = 70.6% is the number
with a control behind it.

The general habit: whenever you loosen a threshold to see what you were missing,
run the loosened setting over shuffled or decoy input as well. Otherwise you
cannot tell recovery from noise.

### 8.3 A negative control must be **inspected**, not counted

Our negative set includes *Bacillus subtilis* 168 with the justification "no
curated ICEberg entry". The module called an ICE there:

```
NC_000964.3|ice-529505:545011   ice   529,362-549,932 (20,571 bp)
anchor_classes integrase,relaxase,t4cp,t4ss     relaxase MOBT   mpf FA
boundary_method tRNA    machinery_intact FALSE   confidence medium
```

That is **ICEBs1** — discovered in *B. subtilis* 168 (Auchtung *et al.* 2005).
The genome's inclusion in the negative set was factually wrong; the call is a
**true positive that ICEberg has not curated**. We came close to recording a
textbook-correct detection as a false alarm. The EBI pipeline independently finds
the same element at bp-identical coordinates.

This is why the headline is worded "**2 calls, 0 confirmed errors** — an upper
bound on false positives, **not** a false-positive rate". If you build your own
negative set, budget time to look at every call, and expect the curated databases
to be incomplete rather than the caller to be wrong.

### 8.4 Two numbers that are easy to quote together and should not be

"15/18 ICEs, 5/12 IMEs, 2 calls in 32.6 Mb" mixes arms. The pilot figures come
from runs with **ICEscan on**; the "2 calls" figure comes from the negative
genomes, which had no `icescan/` directory, so the flag was silently ignored. The
like-for-like negative figure for the pilots' configuration is **5**, not 2. Use
5 whenever you quote them in the same sentence.

The same mixing is what produced the stale "3/12 IMEs without ICEscan" in
[§3.1](#31-master-switches): a before/after where the two halves came from
different builds. **Whenever you quote a before and an after, say which arm and
which build each came from**, or re-run both from the same inputs — which is what
`coverage_sweep/run_coverage_sweep.sh` exists to do.

### 8.5 Turning a knob makes a run non-comparable with the last one

Every database this module *downloads by URL* — TnCentral, ICEberg, ISOSDB and
the ICEscan models — comes from an **unversioned** address, so its four download
rules each write a `PROVENANCE.txt` with source, fetch date, observed checksum
and sequence count. (`rule conjscan_models` is the exception both ways: it
installs a versioned package with `macsydata install` and writes no
`PROVENANCE.txt`.)

**Only two of the four can be pinned:** `mobilome.icescan.sha256` and
`mobilome.tncentral.sha256` are checked against the download and fail the rule on
a mismatch. Set both. ICEberg and ISOSDB have **no `sha256` config key** — their
checksum is recorded after the fact in `PROVENANCE.txt`, so you can say later
what you used but cannot demand it up front. Keep those two `PROVENANCE.txt`
files with the results. If you change
a threshold mid-project, re-run everything rather than mixing — the audit file
records the value that was applied, but only per run, and a half-and-half result
set is very hard to unpick later.

### 8.6 A profile-coverage value outside (0, 1] is silently catastrophic

MacSyFinder accepts `--coverage-profile 30` and then finds nothing, which is
indistinguishable from a genome with no conjugative system. `80_mobilome.smk`
now exits with a message rather than letting that happen — but the general point
stands for any knob you add: **validate it at parse time**, because "found
nothing" is a plausible-looking answer.

---

## Related pages

- `mobilome_worked_example.md` — reading a real report end to end
- `methods_icescan_union.md` — why both model sets are run, and what each contributes
- `methods_att_and_small_plasmids.md` — the *att*-site algorithm, its parameters and its licensing position
- `methods_ebi_comparison.md` — head-to-head against the EBI mobilome-annotation-pipeline
- `mobilome_module_SPEC.md` — the design this implements, including §11 on licensing
- `config/config.yaml` — the `mobilome:` block, where every Layer A default is documented next to its measurement
