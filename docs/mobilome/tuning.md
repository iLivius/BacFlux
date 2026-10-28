# Tuning

The mobilome module uses thresholds to select candidates and describe their context.
These are operational choices, rather than biological constants. Start with the audit
files to identify which rule affected a result before changing a parameter.

This page explains where settings are defined, which observations they affect and what
was measured on the benchmark genomes. Retain the configuration and software version
with each analysis so that results remain interpretable.

| Task | Section |
|---|---|
| explain a filtering decision | [Audit files](#start-at-the-audit-file) |
| locate a setting | [Where settings are defined](#three-places-a-knob-can-live) |
| investigate an unexpected result | [Symptoms and relevant settings](#symptom-to-knob) |
| check defaults and effects | [Parameter reference](#the-knobs-one-by-one) |
| choose reporting criteria | [The two directions](#the-two-directions) |
| account for fragmentation | [Draft assemblies](#on-a-draft-assembly) |
| understand design limits | [Limits of tuning](#what-no-knob-fixes) |

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **AMR** | antimicrobial resistance |
    | **IS** | insertion sequence — the smallest mobile element, carrying only the genes it needs to move itself |
    | **ICE** | integrative and conjugative element — carries its own conjugation machinery, so it can move itself into another cell |
    | **IME** | integrative mobilisable element — carries a relaxase but no apparatus of its own, so it needs a helper element |
    | **AICE** | actinomycete integrative and conjugative element |
    | ***att* site** | the short direct repeat left at each end of an element where it integrated |
    | **HMM** | hidden Markov model — a statistical profile of a gene or protein family |

## Start at the audit file

Six files under `08.mobilome/{sample}/`, one per filtering step. Each carries a reason
for every decision it made.

| File | Written by | What it explains |
|---|---|---|
| `{sample}_ice_discarded.tsv` | `conjscan_to_ice.py` | every ICE/IME candidate dropped, downgraded or flagged |
| `{sample}_amr_mobility_audit.tsv` | `colocalise.py` | every AMR gene that did not get the context or tier you expected |
| `{sample}_is_discarded.tsv` | `isescan_to_table.py` | IS calls dropped or flagged |
| `{sample}_named_elements_discarded.tsv` | `name_transposons.py` | TnCentral hits refused a curated name |
| `{sample}_ice_naming.tsv` | `name_ice_elements.py` | why an ICE candidate did or did not get a curated ICEberg name |
| `{sample}_is_copy_number_audit.tsv` | `isosdb_copy_number.py` | database entries excluded from the copy-number estimate |

Four of the six share one layout — `sample, contig, start, end, action, reason, detail`.
Two do not: `{sample}_amr_mobility_audit.tsv` keys on the gene
(`sample, contig, amr_gene, amr_start, amr_end, mobility_tier, mge_context, decision,
reason, detail`), and `{sample}_is_copy_number_audit.tsv` keys on a database entry with
no coordinates at all (`sample, element, action, reason, detail`).

In all six the last three columns work the same way. `action` (or `decision`) says what
happened, `reason` is the token you filter on, and `detail` carries the numbers in plain
language. Only the column *positions* differ, which is why the recipes below say which
column to cut.

```bash
# What happened to the ICE/IME candidates in this sample?
cut -f5,6 08.mobilome/S1/S1_ice_discarded.tsv | sort | uniq -c | sort -rn

# Why did this AMR gene not get a mobile-element context?
# This file has the wider key block: decision is column 8, reason 9, detail 10.
awk -F'\t' '$8=="no_mge_context" || $8=="tier_not_raised"' \
    08.mobilome/S1/S1_amr_mobility_audit.tsv | cut -f3,9,10
```

!!! warning "The action words are not the same in every file"

    In `{sample}_ice_discarded.tsv` and `{sample}_is_discarded.tsv` a thing that is gone
    from the table says **`dropped`**. In `{sample}_named_elements_discarded.tsv` and
    `{sample}_is_copy_number_audit.tsv` the same event says **`discarded`**. A grep
    across all six for one of those words misses half the answer.

    `kept_flagged` is neither: the thing is in the table, but something about it was
    hedged — usually a confidence cap. `evidence_recorded`, `boundaries_resolved`,
    `not_applicable` and `summary` are notes, written so that a decision is never
    invisible. `row_skipped` and `input_missing` are about the inputs rather than the
    biology: a hit that could not be joined to the annotation, or a file that was not
    there.

!!! note "Look columns up by name in the result tables"

    The audit layouts are stable. The two **result** tables —
    `{sample}_amr_mobility.tsv` and `{sample}_ice_candidates.tsv` — have gained columns
    as the module grew, so tables kept from different builds do not line up.

    ```bash
    awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) col[$i] = i; print; next}
                $col["confidence"] == "high"' any_result_table.tsv
    ```

    One catch with that idiom: if the column is **absent**, `col[...]` is 0 and `$0` —
    the whole line — gets tested instead. An `==` test then matches nothing and hands
    you an empty result; a `!=` test matches everything and the filter silently stops
    filtering. Check the header —
    `head -1 file | tr '\t' '\n' | nl` — whenever the row count looks wrong in either
    direction.

    The two tables also spell their booleans differently: `{sample}_amr_mobility.tsv`
    writes `yes` / `no` / `NA`, `{sample}_ice_candidates.tsv` writes `TRUE` / `FALSE` /
    `NA` so that R reads them straight into a logical column
    (`_tsv_bool` in `conjscan_to_ice.py`).

## Where settings are defined { #three-places-a-knob-can-live }

Some settings are exposed in the user configuration; others require changes to workflow
rules or scripts.

| Layer | Where | How to change it |
|---|---|---|
| **A** | your user configuration, for example `config/config_custom.yaml`, under `mobilome:` | Edit the user copy. Keep the distributed template unchanged. |
| **B** | the `shell:` block of a rule in `workflow/rules/shared/80_mobilome.smk` | Add a supported script argument to the rule; record the change and check it locally. |
| **C** | a constant in `workflow/scripts/80_mobilome/` | Edit the source and review the relevant tests. |

The Layer B flags, all on `rule conjscan_ice` unless noted: `--window-bp`,
`--min-element-bp`, `--min-ime-element-bp`, `--max-element-bp`, `--boundary-bp`,
`--integrase-window-bp`, `--att-flank-window-bp`; plus `--min-length-bp` on
`rule isescan_table` and `--min-chromosome-bp` on `rule mobilome_replicons`.

To change one, add the flag to the rule's shell block:

```python
        shell:
            """
            python {params.script} \
              --sample {wildcards.sample} \
              --conjscan-tsv {input.conjscan_dir}/best_solution.tsv \
              ...
              --replicons {input.replicons} \
              --min-ime-element-bp 3000 \
              {params.strict_boundary} \
              --out-table {output.table} \
              --out-audit {output.audit} > {log} 2>&1
            """
```

Every script the mobilome rules run lives in `workflow/scripts/80_mobilome/`, a
directory named after the rule file that runs it, with each script's tests beside it.
The rules run these scripts in `workflow/envs/python.yaml`, a pinned Python with
nothing else in it, so the scripts use the standard library only.
A script sitting loose in `workflow/scripts/` is not part of a run — it is a maintenance
tool you invoke by hand. `workflow/scripts/README.md` has the full map. After a Layer C
edit, `pytest workflow/scripts` is the check that you did not break something else on
the way; the tests encode today's thresholds in places,
so expect to update a test alongside a deliberate change.

!!! warning "Layer B values are not sanity-checked"

    The scripts check only that the value is a whole number. `--min-ime-element-bp
    80000` is accepted, silently drops every IME-architecture call, and the result is
    indistinguishable from a genome that has none.

    `mobilome.coverage_profile` is the exception, and its guard lives in the workflow
    rather than in the script (the check right after `MOBILOME_COVERAGE_PROFILE`, at
    the top of `80_mobilome.smk`): MacSyFinder accepts
    `--coverage-profile 30` and then finds nothing, which looks exactly like a genome
    with no conjugative system, so BacFlux exits at parse time instead — when Snakemake
    reads the config, before any rule runs. The general lesson for any parameter you add:
    check it there too, because "found nothing" is a plausible-looking answer.

!!! note "Most thresholds reach the audit file only when they fire"

    The value appears inside the `detail` of the row it caused —
    `cluster_shorter_than_min` names the floor and the measured span,
    `integrase_without_conjugation_machinery` names the integrase window. The ICE and IS
    audits have no unconditional "settings used" row, so on a sample where nothing was
    dropped neither the audit nor the log says which floors were in force. The one
    exception is the naming layer: `{sample}_named_elements_discarded.tsv` always closes
    with a `summary` / `tncentral_naming_complete` row naming both thresholds it applied.
    The parse-time banner prints only `max_composite_span_bp` and `contig_boundary_bp`
    (the `Mobilome module: ON` line, printed from `00_common.smk`). If you change a
    Layer B parameter, write the value down yourself.

## Symptoms and relevant settings { #symptom-to-knob }

Every `reason` below is a literal string in the code — grep for it.

### An element you know is there was not reported

Grep `{sample}_ice_discarded.tsv`:

| Reason | What it means | Parameter |
|---|---|---|
| `cluster_shorter_than_min` | the machinery span was below the size floor; `detail` gives the measured span and which floor applied | `--min-ime-element-bp` for IME/AICE architecture, `--min-element-bp` otherwise |
| `cluster_longer_than_max` | anchors chained across the contig into one implausible element | `--max-element-bp`, or **lower** `--window-bp` so the chain breaks |
| `no_conjugation_anchor` | integrases with no relaxase and no mating apparatus nearby — a genome's ordinary complement of recombinases | not a threshold problem |
| `integrase_without_conjugation_machinery` | the integrase marking your element was attached to no machinery cluster | `--integrase-window-bp`, but read the next box first |
| `untrusted_integrase_profile` | an ICEscan integrase profile deliberately ignored (integron and housekeeping recombinases that would anchor elements onto the wrong gene) | Layer C: `ICESCAN_TRUSTED_INTEGRASE_PROFILES` in `conjscan_to_ice.py` |
| `integrase_match_is_a_transposase` | the Bakta product text matched the integrase pattern but the protein is a transposase | Layer C: `TRANSPOSASE_PRODUCT_PATTERN` in `conjscan_to_ice.py` |
| `conjscan_found_no_systems` | MacSyFinder assembled nothing at all | `mobilome.coverage_profile`, then the ICEscan layer |
| `hit_id_not_in_annotation` | a MacSyFinder hit could not be joined to the Bakta GFF3 | not a threshold — a mismatched input pair. Check that the Bakta run and the `.faa` MacSyFinder read are the same run |

In `{sample}_named_elements_discarded.tsv`: `reference_coverage_below_threshold` →
`mobilome.tncentral.min_reference_coverage`, and `identity_below_naming_threshold` →
`mobilome.tncentral.min_identity`. A third, `tncentral_hit_is_a_plain_is`, is not a
threshold at all: an insertion sequence carries only what it needs to move, so it cannot
hold a passenger AMR gene, and ISEScan has already inventoried it.

If the element is there but no **AMR gene** picked it up, look in
`{sample}_amr_mobility_audit.tsv` for `near_but_outside_element_machinery_span` (the
element's reported interval is the machinery span and your gene is in the cargo beyond
it — a boundary problem, not a threshold one) or `nearest_mobile_element_too_far`
(`IS_ADJACENT_MAX_BP`, Layer C).

!!! warning "Read the audit before you turn the parameter the literature points at"

    A worked case. *Streptococcus salivarius* CP002888.1 carries the
    curated IME_Ssal57I_tRNAlys at 70,384–75,506. This clade is exactly the `T4SS_MOBT`
    biology that profile coverage is supposed to hold back, so lowering
    `coverage_profile` ought to recover it. Measured at 0.5, 0.4 and 0.3 it is missed at
    every value, and the audit says why in two lines: the relaxase hits were being found
    at 0.5 all along, a machinery cluster took an integrase sitting 324 bp away, and the
    integrase that actually marks the curated element — 20 kb further on — was left with
    nothing to attach to. Each integrase joins at most one cluster.

    Not `coverage_profile`, and not `--integrase-window-bp` either: the element's
    integrase was well inside the 50 kb window. It lost to a nearer competitor. That is
    a co-localisation limit, not a threshold.

### The interval looks far too big, or far too small

| You see | Read | Then |
|---|---|---|
| interval much **smaller** than the published element | `boundary_method` | if `none`, the interval **is** the machinery span — a floor, by design. Compare `start`/`end` against `machinery_start`/`machinery_end`. `--att-flank-window-bp` is the parameter, but see below |
| interval much **larger** | `n_anchors`, `anchor_classes`, and `evidence_recorded / system_merge_refused_too_long` | lower `--window-bp` (fuses less) or `--max-element-bp` (drops the chain) |
| a de novo repeat found but the interval did not move | `kept_flagged / denovo_att_reported_not_applied` | working as intended — a de novo repeat is reported and never applied. No parameter; this is a design decision, not a threshold |
| ends resolved but at the wrong place | `att_length_bp`, `att_trna`, `attL`, `attR` | Layer C: `MIN_ATT_REPEAT_BP`, `MAX_ATT_COPIES_OUTSIDE_TRNA` |

### Everything comes out at medium confidence

Check which cap fired — the audit names one per rule:

```bash
awk -F'\t' '$5=="kept_flagged"' 08.mobilome/S1/S1_ice_discarded.tsv | cut -f6 | sort | uniq -c
```

| Cap reason | Meaning | Parameter |
|---|---|---|
| `few_anchor_classes` | fewer than 3 of the five anchor classes were present | `MIN_ANCHOR_CLASSES_FOR_HIGH` (Layer C, in `conjscan_to_ice.py`) — read the box first |
| `machinery_not_intact` | low system wholeness, a decayed model, or a truncated relaxase/VirB4 | `WHOLENESS_INTACT_MIN` / `PROFILE_COVERAGE_INTACT_MIN` (Layer C, both 0.7, in `conjscan_to_ice.py`) |
| `at_contig_boundary` | machinery runs into a contig end → capped **low** | `--boundary-bp` |
| `spans_contigs` | hits on more than one contig → capped **low**, absolutely | none, deliberately |
| `no_assembled_system` | profile hits only, no assembled system → capped **low** | `mobilome.coverage_profile`, or accept the tier |
| `boundary_not_resolved` / `boundary_denovo_only` | only fires when `require_trna_boundary_for_high` is on | turn that key off |
| `contig_length_unknown` | no length for the contig, so it cannot be said whether the element runs off the end | fix the `--contig-lengths` input |

!!! note "IME confidence on the benchmark"

    The confidence rule is the same for every class. High confidence requires at
    least three distinct anchor classes (`MIN_ANCHOR_CLASSES_FOR_HIGH = 3`), intact
    machinery and no other confidence cap. The five anchor classes are integrase,
    relaxase, coupling protein, mating apparatus and AICE translocase. An IME has an
    integrase and a relaxase; a third class can come only from a CONJscan coupling
    protein, or from a VirB4 that CONJscan's MOB model lists as accessory. ICEscan's
    coupling-protein and mating-pair hits are never used as anchors.

    None of the 26 IME rows in the benchmark reached high confidence: 25 had two
    classes and one had three but remained medium. Eighteen also carried
    `machinery_degraded` because of `low_system_wholeness`: a relaxase-only match to
    a larger CONJscan model can score 0.333 wholeness even when the IME architecture
    is expected. In comparison, 30 of 53 ICE rows were high confidence.

    Lowering the anchor-class minimum to two would make only eight of the 26 IMEs
    eligible: 17 others are also capped by `machinery_not_intact` and one by
    `at_contig_boundary`. It would also change the criterion for two-anchor ICE-like
    clusters. Read `mobility`, `machinery_intact` and the audit reasons together.

### Calls on a genome that should have none

Read the **class** before the count.

| You see | Do |
|---|---|
| `evidence_level=profile_hits_only` | the module's weakest evidence: the profiles were seen, no system was assembled. Filter these out for reporting, or raise `coverage_profile` |
| `mge_class=cime_or_island`, `mobility=passive` | no transfer capability is asserted and no tier is raised: a gene inside stays at tier 1 with the audit reason `inside_context_only_element`. The candidate interval still needs evaluation |
| a call at medium or better with a real mobility claim | investigate before dismissing |
| calls in an archaeon or a very distant lineage | the CONJscan and ICEscan models are bacterial. Treat any call as suspect on those grounds alone |

## Parameter reference { #the-knobs-one-by-one }

Defaults are as shipped. "Measured" means there is a run behind the claim; where nothing
has been swept, this page says so rather than offering a plausible-sounding
recommendation.

### Profile coverage: `coverage_profile` { #coverage_profile-the-one-layer-a-knob-that-changes-what-is-detected }

The fraction of an HMM profile a protein must align to before the hit is kept
(`mobilome.coverage_profile`, default **0.5**, read into `MOBILOME_COVERAGE_PROFILE` at
the top of `80_mobilome.smk`). Not identity, not
an E-value: a hit can be statistically beyond doubt and still be dropped by this rule
alone. It is applied to **both** MacSyFinder searches deliberately — their hit tables
are merged, so different stringencies would make an element's class depend on which
model set was looser.

Lowering it admits divergent full-length relaxases that share only a catalytic core, and
also genuine fragments. Measured end to end over 18 curated ICEs, 12 curated IMEs and 12
genomes with no curated element (32.6 Mb), with the ICEscan layer on (the same table
sits in the comment above `coverage_profile` in `config/config.yaml`):

| `coverage_profile` | curated ICEs | curated IMEs | calls on the negative set |
|---|---|---|---|
| **0.5** (shipped) | 15 of 18 | 5 of 12 | 5 (2 without ICEscan) |
| 0.4 | 15 of 18 | 6 of 12 | 7 (4 without ICEscan) |
| 0.3 | 15 of 18 | 6 of 12 | 9 (5 without ICEscan) |

Across 30 curated elements, lowering recovers exactly **one**: a 23 kb IME in
*Faecalibacterium duncaniae*, found in full (right edge 53 bp off the published
coordinates) and correctly classed. It appears at 0.4 and gains nothing further at 0.3.

The cost at 0.4 is two more calls, in *E. coli* K-12 MG1655 and *P. aeruginosa* PAO1.
Read the class first: every call added at 0.4 is `cime_or_island` /
`passive` / `low` / `profile_hits_only`, claiming no mobility and carrying no tier. At
0.3 that stops holding — a 21.6 kb IME appears in *S. aureus* N315 at **medium**
confidence with a real mobility claim, unverified either way.

**When to lower it anyway.** Across all 395 curated IMEs, 73 carry a relaxase hit that
is clean on E-value and fails only this rule — 60 of them the `T4SS_MOBT` family at
about 0.32 coverage, largely the *Streptococcus salivarius* clade (the same comment
in `config/config.yaml`). If your isolates sit in that space, 0.4 is worth trying.
Nothing measured supports 0.3 for single isolates. The EBI mobilome-annotation-pipeline
runs 0.3, which reflects a metagenome pipeline's priorities — recall first, because MAG
proteins are fragmentary anyway — not a verdict on what a single-isolate workflow should
do.

!!! note "Why the gain is so small when the raw evidence moves a lot"

    0.4 adds 5% more MacSyFinder hits and changes the hit table in 13 of 40 genomes; 0.3
    adds 12% and changes 19 of 40. Almost none of it reaches the element table, because
    the caller needs an integrase within 50 kb of conjugation machinery before it seeds
    anything. **Coverage sits upstream of a stronger constraint, and co-localisation is
    what binds.**

### Element geometry: script arguments { #element-geometry-real-knobs-but-not-config-keys }

`rule conjscan_ice` passes none of these, so `conjscan_to_ice.py`'s own defaults are
what every measured result in this project was produced at. They are Layer B.

| Flag | Default | Raising it | Lowering it |
|---|--:|---|---|
| `--window-bp` | 15,000 | chains more machinery genes into one element; risks fusing two neighbouring systems | splits one operon into several candidates, each then too small to survive the size floor |
| `--min-element-bp` | 8,000 | fewer, larger ICE candidates | admits short machinery spans; `cluster_shorter_than_min` stops firing |
| `--min-ime-element-bp` | 2,000 | the parameter with the largest measured effect — see below | admits machinery spans below the empirical 2 kb floor |
| `--max-element-bp` | 500,000 | admits runaway anchor chains | drops genuinely large ICEs |
| `--boundary-bp` | 1,000 | more candidates flagged as probably truncated and capped at low | fewer caps, and you stop being told that a call ran into the end of its contig |
| `--integrase-window-bp` | 50,000 | an integrase further from the machinery may still anchor the element | more `integrase_without_conjugation_machinery` orphans, so fewer ICE/IME calls |
| `--att-flank-window-bp` | 50,000 | more candidate repeats to rank | the element's real ends may sit outside the window |

The defaults are the `DEFAULT_…` constants of `conjscan_to_ice.py`
(`DEFAULT_WINDOW_BP`, `DEFAULT_MIN_ELEMENT_BP`, `DEFAULT_MIN_IME_ELEMENT_BP`,
`DEFAULT_MAX_ELEMENT_BP`, `DEFAULT_BOUNDARY_BP`, `DEFAULT_INTEGRASE_WINDOW_BP`) and
`DEFAULT_FLANK_WINDOW_BP` in `att_search.py`.

**Why there is a separate IME floor.** The floor is applied to the *anchor-cluster span*,
and that span scales with the number of machinery genes, not with the element's length.
An ICE carries a roughly 20-gene mating-apparatus operon, so its span is naturally tens
of kb. An IME carries a relaxase and an integrase — two genes — so its span is 1–6 kb.
An 8,000 bp floor selects for ICEs by construction. AICEs take the lower floor for the
same reason.

Measured over the 28 small-machinery calls (IME + AICE) in the benchmark:

| floor | calls dropped (of 28) | curated IMEs still detected (of 12) |
|--:|--:|--:|
| **2,000** (default) | 0 | 5 |
| 3,000 | 5 | 5 |
| 5,000 | 7 | 3 |
| 8,000 (= the ICE floor) | 15 | 1 |

The recall cliff is between 3,000 and 5,000, not at the default: the five detected
curated IMEs have machinery spans of 3,840 / 4,959 / 5,942 / 6,252 / 8,164 bp (listed
in the note on the ICE/IME geometry thresholds in `config/config.yaml`). **3,000 bp
drops five calls at a measured cost of zero curated IMEs** — that is where to go if
2,000 feels unevidenced. It does not make 3,000 the better default: five calls is a small sample and none of the five is a known false
positive.

!!! warning "2,000 bp is an empirical cut, not a published bound"

    There is no published minimum size for an IME's machinery. The value was chosen
    against a handful of observations in one benchmark. Say so in a methods write-up
    rather than implying it is a property of IMEs.

!!! note "The att search is always given the 8,000 bp floor"

    It uses the 8,000 bp ICE floor rather than the 2,000 bp IME machinery floor,
    and rejects any repeat pair implying an element below 8 kb, for every candidate.
    A small machinery span does not imply an element of the same length, but a
    genuinely smaller element cannot have its own boundaries accepted. Three of the
    five detected curated IMEs are shorter than 8 kb (IME_CdiR20291_ND 5,551 bp,
    Tn*4451* 6,338 bp, MTn*Pi*10 7,591 bp), so a pair at their true ends would be
    rejected; the search then ranks the remaining pairs, and a pair implying 8 kb or
    more can still win. Neither floor is a configuration key, and a per-class *att*
    floor has not been measured. Read `boundary_method` and the reported coordinates
    rather than assuming that admission as an IME also means its complete extent was
    resolved.

`--integrase-window-bp` is deliberately much larger than `--window-bp`: conjugation
machinery is an operon and clusters tightly, whereas the integrase sits at the element
boundary, tens of kb away on a large ICE. It was widened to 50 kb after a 15 kb window
produced `conjugative_region` instead of ICE on two clinical *K. pneumoniae* isolates —
on both, the ICE*Kp* integrase sits 33,361 bp from the machinery cluster — and put the
ICE label on a different element 1.1 Mb away.

`--att-flank-window-bp` was 30,000 until SPI-7 showed why that was too narrow: its
*attR* sits 5,800 bp outside a 30 kb window, so the element was reported 50 kb short.
80 kb, 120 kb and 200 kb change no element's answer on the benchmark, so 50 kb is
where the curve flattens rather than an arbitrary larger number.

### The composite span, and the co-localisation constants

`mobilome.max_composite_span_bp` (Layer A, default **20,000**, `MOBILOME_MAX_COMPOSITE_SPAN`
in `00_common.smk`) is
the longest span accepted for a composite transposon, which is ladder tier 3. Real
composites run from about 2.5 kb (IS*26* translocatable units) to about 25 kb. The
measured span is always reported next to the call, in `flanking_span_bp`.

Everything else in the co-localisation step is Layer C, in `colocalise.py`:

| Constant | Default | Meaning |
|---|--:|---|
| `HYBRID_PROMOTER_MAX_BP` | 500 | how close an upstream, correctly oriented IS must be before the module says it may be driving the gene from an outward-reading promoter (tier 2). IS*Ecp1* sits 42–266 bp upstream of *bla*<sub>CTX-M</sub> |
| `IS_ADJACENT_MAX_BP` | 5,000 | beyond this the nearest IS is still reported in `distance_bp` but no longer sets `mge_context` |
| `CONTIG_END_WINDOW_BP` | 1,000 | triggers the boundary flag and the low-confidence cap on an AMR row |
| `MIN_AMR_FRACTION_INSIDE_ELEMENT` | 0.9 | fraction of the AMR gene inside a named element before the element is said to carry it |
| `MIN_AMR_FRACTION_OVERLAPPED_FOR_DISRUPTION` | 0.5 | fraction overlapped by an IS before the CDS is called disrupted (likely inactivated) rather than merely overlapping |
| `ABUTTING_IS_MAX_GAP_BP` | 25 | how flush an IS must be against a **partial** AMR hit to read the pair as "the IS split this gene". Tight on purpose: it tests "butted up against", not "nearby" |
| `MIN_COVERAGE_FOR_INTACT_GENE` | 90.0 | below this percent coverage an AMRFinderPlus hit is treated as a fragment. Used only with the abutting-IS test |
| `ORIENTATION_EXEMPT_IS_FAMILIES` | `{IS6, IS26, …}` | families exempt from the same-orientation rule for composites |
| `UNINFORMATIVE_FAMILY_VALUES` | `{"", NA, NEW, UNKNOWN, ISNCY, …}` | ISEScan family labels carrying no information about *which* element this is |

!!! warning "IS26: keep IS6 in `ORIENTATION_EXEMPT_IS_FAMILIES`"

    For every other family, a flanking pair qualifies for tier 3 only when both
    copies carry a strand and the strands agree (`same_orientation = yes`). Opposite
    strands are rejected with `flanking_is_pair_inverted_orientation`, a missing strand
    with `flanking_is_strand_unknown`, both under `tier_not_raised`; composites with
    inverted IS copies, such as Tn*5* and Tn*10*, are therefore not called. Pairs of
    family IS6, or copies named IS26, skip this test. Direct pairs pass it anyway; the
    exemption adds inverted pairs (`same_orientation = no`) and pairs with an unknown
    strand (`NA`). The family and span tests still apply, a partial flanking copy still
    caps confidence at `medium` (`flanking_is_partial`), and every such composite
    records `is26_orientation_exemption_applied`. Removing the family from
    the set would drop IS26-bounded genes with inverted or strand-unknown copies below
    tier 3. A matching IS pair is evidence of a candidate architecture; it does not
    demonstrate active transposition.

    `UNINFORMATIVE_FAMILY_VALUES` prevents two unclassified IS calls from supporting
    a composite solely because both have the same uninformative label.

None of the Layer C values here have been swept. They are documented so that you can see
what a given audit reason is testing, not because there is evidence that a different
value is better.

### The naming thresholds

Both naming layers are *optional* and stay off until you give each a source — see
[Enabling the module](enabling.md). Once on, four Layer A keys decide what gets a name.

| Key | Default | Source | Meaning |
|---|--:|---|---|
| `mobilome.tncentral.min_identity` | 90.0 | `TNCENTRAL_MIN_IDENTITY`, `00_common.smk` | percent identity before a BLAST hit may confer a curated **name**. Below it the sequence may be a relative of the element, but the name would claim more than the data shows |
| `mobilome.tncentral.min_reference_coverage` | 0.80 | `TNCENTRAL_MIN_COVERAGE`, `00_common.smk` | fraction of the **reference element** that must be present — "is the whole of this known transposon here?" |
| `mobilome.iceberg.min_identity` | 80.0 | `ICEBERG_MIN_IDENTITY`, `00_common.smk` | identity floor for an ICE name |
| `mobilome.iceberg.min_overlap_fraction` | 0.50 | `ICEBERG_MIN_OVERLAP`, `00_common.smk` | how much of *our* candidate the curated element must cover to be treated as the same thing. Lenient next to the transposon cascade, because ICEs are mosaic and their cargo varies between strains |

ICEberg naming labels an existing candidate without changing any AMR gene's tier.
TnCentral naming can affect tier 4: accepting a named transposon can raise the tier of
a gene it contains when no higher-priority context applies. Changing these thresholds
therefore affects more than names in the TnCentral branch.

`min_reference_coverage` is the one with a measured story, and it is why tier 4 is hard
to reach on a fragmented assembly. On a clinical hybrid run, every Tn*Ecp1.1* candidate
was refused:

```text
TnEcp1.1: only 49% of the 3417 bp reference element is present (needs 80%).
A fragment of a transposon is not that transposon - on a short-read assembly
this usually means the element is split across contigs.
```

and a second at 12%. That refusal is correct behaviour, but it means the gene stays at
tier 2 and you have to open `{sample}_named_elements_discarded.tsv` before concluding
that a genome holds no named transposon. Lowering the threshold to admit a 49% fragment
has **not been measured**, and it means asserting a named architecture from half the
evidence; if you do it, read the audit and state in the methods what coverage you
accepted.

### The three contig-end windows

They answer different questions about different objects, which is why they are separate
numbers and only one of them is a config key.

| Window | Value | Where | Question it answers |
|---|--:|---|---|
| `mobilome.contig_boundary_bp` | 100 | Layer A, `MOBILOME_BOUNDARY_BP` in `00_common.smk`; passed to `isescan_to_table.py` as `--boundary-bp` | is an **IS call** sitting on the edge of its contig? Feeds `at_contig_boundary` and `fraction_at_contig_boundary` in the IS summary |
| `CONTIG_END_WINDOW_BP` | 1,000 | Layer C, `colocalise.py` | was there enough flanking sequence to have **seen** a neighbouring element at all? Caps confidence in the mobility table |
| `--boundary-bp` on `conjscan_ice` | 1,000 | Layer B, `DEFAULT_BOUNDARY_BP` in `conjscan_to_ice.py` | the same question again, for an **ICE/IME candidate** rather than an AMR gene |

An IS 500 bp from a contig end is not "at the edge", but you still could not have seen
its partner 2 kb away — so absence of evidence there is not evidence of absence. That is
why the second window is ten times the first. Note that the first and third share the
flag name `--boundary-bp` while meaning different things, and only the first is wired to
a config key.

A high `fraction_at_contig_boundary` indicates that many IS calls lie at contig ends,
where fragmentation may obscure further copies or context. The located count can
underestimate the true count. For reference, the hybrid assemblies of the clinical *K.
pneumoniae* isolates TUM24772 and K3, near-closed genomes, score 0: none of their 95 and
110 IS calls sits at a contig end. A fragmented draft may have a higher value; inspect
the sample-specific summary.

ISEScan itself is run with **no** `--removeShortIS`, on purpose: partials are kept and
tiered rather than discarded. `--min-length-bp` on `rule isescan_table` defaults to **0**
(the flag's default in `isescan_to_table.py`) — set it to 500 if you want the
discard-below-500-bp convention applied to IS calls too.

Two kinds of call are dropped for a biological reason rather than a threshold. The
first is an IS call lying over the *dnaA* gene. *dnaA* is essential, so no living isolate
carries an insertion inside it. ISEScan reports a partial IS21 there, most likely because
IS21 elements carry IstB, a helper protein related to DnaA, so the IS21 profile matches
DnaA itself. `rule isescan_table` reads Bakta's annotation table to find *dnaA*, and a
call of any family goes to `{sample}_is_discarded.tsv` with reason `overlaps_dnaA` when
at least half of it lies inside a CDS Bakta calls *dnaA*, or at least half of that CDS
lies inside the call. The second case is the call that runs on into *dnaN*, the next
gene, because ISEScan's gene finder joined the two. That half is
`DNAA_OVERLAP_FRACTION`, a constant in `isescan_to_table.py`, not a configuration key.

Two of ISEScan's IS21 models (IS21_259 and IS21_288) also match the chromosome's
replicative helicase loader, DnaC in *E. coli*, *Salmonella*, *Klebsiella* and their close
relatives and DnaI in the Bacilli, and that gene is just as essential. Its name
alone cannot single it out: phages, integrative elements and plasmids carry loaders of
their own, and Bakta can name the IstB of a real IS21 element *dnaC*. What marks the
chromosome's own loader is its neighbour. A call goes to `{sample}_is_discarded.tsv`
with reason `overlaps_helicase_loader` only when all four of these hold:

1. ISEScan calls it a **partial** IS21. A complete IS21, or a call of any other family,
   is never dropped this way.
2. It lies over a CDS that Bakta names *dnaC* or *dnaI*, by the same half-of-the-shorter
   test as *dnaA*.
3. That CDS has its chromosomal partner on the same contig within 5 kb: *dnaT* with the
   product "primosomal protein DnaT", the other half of the *dnaC*–*dnaT* operon of
   *E. coli* and its relatives, which sits 0–2 bp away; or *nrdR*, at the head of the
   *nrdR*–*dnaB*–*dnaI* cluster of the Bacilli, about 1.2–3.2 kb away. The DnaT-like
   proteins of mobile elements, which Bakta also names *dnaT* but describes as "DnaT
   DNA-binding domain-containing protein", do not count.
4. No IS21 transposase lies on the same contig within 5 kb of the loader: no CDS named
   *istA*, and none whose product mentions IS21, the word IstA, Mu transposase,
   integrase catalytic or integrase core. When one is there, the "loader" may be a
   real IS21 element's IstB, and the call is kept. A transposase of another family does
   not count, because unrelated transposases do sit a few kb from true chromosomal
   loaders.

The audit row names the loader, its locus tag and coordinates, the partner and its
distance, and how much of the call and of the gene they share. Every other partial
IS21 over *dnaC* or *dnaI* is kept and counts in the IS totals: the loaders of phages,
ICEs and plasmids, which sit nowhere near a *dnaT* or *nrdR*; a real IS21 element whose
IstB Bakta names *dnaC*; and the chromosome's own loader where the partner is not on
the same contig, either because the taxon lacks it (Paenibacillaceae, Clostridia,
*Aquifex*) or because a fragmented draft put it on another contig. The run log of `rule
isescan_table` says, for every *dnaC* and *dnaI*, which way it went and why. The two
5 kb windows and the gene names are constants in `isescan_to_table.py`
(`LOADER_PARTNER_*`, `IS21_TRANSPOSASE_*`), not settings meant to be tuned.

### Boundaries: the att search

| Parameter | Layer | Default | Source | Meaning |
|---|:-:|--:|---|---|
| `mobilome.require_trna_boundary_for_high` | A | `false` | `MOBILOME_REQUIRE_TRNA_BOUNDARY`, `00_common.smk` | whether an element must have tRNA-anchored ends before it may be called **high** confidence |
| `MIN_ATT_REPEAT_BP` | C | 15 | `att_search.py` | shortest exact repeat that may be an *att* site. Matches ICEfinder2; DBSCAN-SWA uses 12 |
| `DENOVO_MAX_EXPECTED_CHANCE_MATCHES` | C | 0.05 | `att_search.py` | how many chance repeats to tolerate when picking the shortest believable length for *this* pair of flanks |
| `MAX_ATT_COPIES_OUTSIDE_TRNA` | C | 2 | `att_search.py` | copies of the repeat allowed outside any annotated tRNA before it is judged a repeat *family* rather than a single integration scar |
| `MAX_SEED_MATCHES` | C | 200,000 | `att_search.py` | compute guard in repeat-dense regions |

`require_trna_boundary_for_high` folds two different questions into one number: *is this
an ICE?* (anchors, machinery, one contig) and *where does it end?* (`boundary_method`).
By default they are reported side by side and only the first sets confidence.

Measured cost of turning it on: across the 85 element rows produced by the 52 complete
records in the benchmark, 30 reach high confidence and only **7** of those have
`boundary_method=tRNA`. The key demotes 23 of 30 high calls to medium **even on closed
genomes**. Use it on closed long-read assemblies, where an unresolved boundary is a
warning sign rather than the norm.

Either way, cargo is never assigned from an unresolved boundary: a de novo repeat is
reported but never widens an element, so the interval stays the machinery span, a floor.

!!! note "Limits and ranking of the att search"

    **The length floor is not the guard that does the work.** The chance-match
    calculation resolves to 18 bp at the shipped 50 kb window — and returned 18 at the
    old 30 kb window too, so widening the window does not move it
    (`minimum_informative_repeat_length` in `att_search.py`). It also does not filter
    much: measured on the *K. pneumoniae* positive-control chromosome over 300
    randomly placed non-ICE spans,
    the length threshold alone let 22% of them return a confident boundary against the
    1.3% the formula predicts. The guard that actually works is the copy count,
    `MAX_ATT_COPIES_OUTSIDE_TRNA` — an rRNA operon repeat clears every length threshold
    and fails the copy count cleanly.

    **Ranking is by anchoring first, then length** — there is no tRNA score bonus. A
    fixed bonus was tried and removed, because the two properties are not comparable
    quantities: repeat length says how unlikely a match is by chance, sitting at a tRNA
    3′ end says the match is where integration actually happens. Measured on SPI-7, a
    51 bp repeat in ordinary sequence beat the real 24 bp *att* pair at tRNA-Phe and the
    element came out 50 kb short (`search_maximal_repeat` in `att_search.py`).

    **`MAX_SEED_MATCHES` truncates silently.** Hitting it breaks out of the scan, so a
    real repeat sitting among the unexamined seeds is missed and the element reports
    `boundary_method=none` — indistinguishable from "the flanks carry no repeat". No
    audit row says the search was cut short. Across 246 att searches over 52 benchmark
    genomes no window reached the cap, so this is a latent gap rather than a live
    problem (the note above `MAX_SEED_MATCHES` in `att_search.py`).

### IS copy number and replicon calling

| Parameter | Layer | Default | Source | Meaning |
|---|:-:|--:|---|---|
| `mobilome.isosdb.min_covered_percent` | A | 90.0 | `ISOSDB_MIN_COVERED`, `00_common.smk` | a database entry must be covered end to end before its depth is believed. A partially covered entry is usually a conserved domain shared with another family, and averaging it in inflates every estimate |
| `mobilome.isosdb.min_copy_number` | A | 0.5 | `ISOSDB_MIN_COPIES`, `00_common.smk` | below this multiple of the genome baseline the element is treated as absent. Under 1.0 on purpose: a real single-copy IS sits near 1×, and sampling noise plus mapping loss routinely push it to 0.6–0.8× |
| `--min-chromosome-bp` | B | 2,000,000 | `DEFAULT_MIN_CHROMOSOME_BP`, `platon_replicons.py` | above this length an unclassified contig may be called the chromosome |

The copy-number analysis is *optional* (`mobilome.isosdb.fasta_url`) and runs in
Illumina and hybrid modes only, because it needs reads. It changes no AMR gene's tier
and must not: it says nothing about *where* the extra copies are, only that their
estimated abundance exceeds the located count.

`--min-chromosome-bp` matters more than it looks. An ICE integrates into a chromosome by
definition, so machinery on a contig called `plasmid` is a conjugative plasmid region,
not an ICE — audit reason `ice_demoted_on_plasmid_replicon`. On an `unknown` contig you
get `ice_on_unclassified_replicon` instead. The floor is not Platon's 500 kb because
megaplasmids run to 1–2 Mb; between the two the module says `unknown` and lets the
confidence cap reflect that.

## The two directions

### Tightening for a dossier

!!! note "Tier 6 and stricter reporting"

    Tier 6 means the gene lies inside an element classed `ice` (at least 90% of the
    gene inside the interval) or on a plasmid Platon typed conjugative. Settings that
    change which ICE intervals exist, or how far they reach, change it: in Layer A
    `coverage_profile` and `icescan.run`; in Layer B, on `rule conjscan_ice`,
    `--window-bp`, `--min-element-bp`, `--max-element-bp`, `--integrase-window-bp`
    and `--att-flank-window-bp`; and `--min-chromosome-bp` on
    `rule mobilome_replicons`, because an ICE on a contig called a plasmid becomes a
    conjugative region. `--min-ime-element-bp` acts only on IME and AICE clusters,
    and the naming thresholds act on tier 4 (TnCentral) or on names alone (ICEberg).
    A high-confidence tier-6 call rests on rules that are not settings: three anchor
    classes, intact machinery, and the `spans_contigs` and `at_contig_boundary`
    caps. No threshold provides a calibrated false-positive rate for this tier.

    In the 40-genome comparison at profile coverage 0.5, enabling ICEscan left the 37
    ICE rows and their confidence labels unchanged: 16 high, 19 medium and two low.

    For a report restricted to stronger tier-6 predictions, use explicit criteria
    and disclose the excluded rows. For example, after confirming these columns
    are present:

    ```bash
    # The mobility table uses yes / no / NA for boolean fields.
    awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) col[$i] = i; print; next}
                $col["mobility_tier"] == 6 &&
                $col["confidence"]    == "high" &&
                $col["spans_contigs"] != "yes"' \
        08.mobilome/S1/S1_amr_mobility.tsv
    ```

    This selects predictions meeting the stated criteria. It does not validate
    transfer experimentally or establish a false-positive rate.

The levers that do help, and what each costs:

| Change | Effect | Measured cost |
|---|---|---|
| leave `coverage_profile: 0.5` | keeps MacSyFinder's own default | nothing measured supports going lower for a single isolate |
| raise the naming thresholds (`tncentral.min_identity` to 95.0, `min_reference_coverage` to 0.9, `iceberg.min_identity` to 90.0) | requires stronger similarity or coverage for a curated name | not measured; TnCentral changes can also affect tier 4 |
| `--min-ime-element-bp 3000` (Layer B) | drops 5 of 28 small-machinery calls | **zero** curated IMEs on this evidence |
| `--min-ime-element-bp 8000` (Layer B) | stops reporting the IME architecture | curated IME recall 5/12 → 1/12, 15 of 28 calls gone. Buys the removal of one IME call on the negative set |

The audit will then hold more `dropped` rows, dominated by `cluster_shorter_than_min`
with the span and the threshold in `detail`, and more
`reference_coverage_below_threshold` and `identity_below_naming_threshold` in
`{sample}_named_elements_discarded.tsv` — where the action word is `discarded`, not
`dropped`. Small IME-architecture clusters stop being reported at all, so **absence of an
IME row is not evidence of absence**; the audit is where they went, and a dossier should
say so.

### Widening for a survey

Use this where missing a real element is the expensive error — an exploratory survey, or
an organism nobody has curated.

```yaml
mobilome:
  run: true

  # 0.4, not 0.3. 0.4 recovers one curated IME and every call it adds to the
  # negative set is passive / low / profile_hits_only. 0.3 adds nothing further
  # and produces the first negative-control call that actually asserts mobility.
  coverage_profile: 0.4

  # A wider composite window. A convention, not biology, and the measured span
  # is always in flanking_span_bp next to the call.
  max_composite_span_bp: 25000
```

Turning the ICEscan layer on roughly doubles the IME rows and is the only way to get an
AICE call at all; it needs its own source, which [Enabling the module](enabling.md) covers.
Measured against an otherwise unchanged configuration rather than a half-changed one — the
"before" being CONJscan alone at 0.5 and the "after" the union at 0.4 — IMEs go 2/12 to
6/12, ICEs stay at 15/18, and negative-set calls go 2 to 7.

The table is now a **candidate list to triage**, not a set of findings. Filter it before
reading, and treat nothing at `evidence_level=profile_hits_only` as a finding without
looking at the genes yourself:

```bash
# Everything that actually claims mobility, with a system behind it.
awk -F'\t' 'NR==1 {for (i = 1; i <= NF; i++) col[$i] = i; print; next}
            $col["evidence_level"] != "profile_hits_only" &&
            $col["mge_class"]      != "cime_or_island"' \
    08.mobilome/S1/S1_ice_candidates.tsv
```

The audit will hold far fewer `dropped` rows and many more `kept_flagged` ones —
`no_assembled_system`, `few_anchor_classes`, `mpf_marker_without_typed_system`,
`machinery_not_intact`. That is the shape of a wider net working correctly: the extra
material arrives already labelled weak.

## On a draft assembly

Fragmentation is not a threshold problem, and almost nothing on this page fixes it. The
benchmark cut 40 genomes to three contiguity levels and reran the analysis. Element
classes were generally more stable than their reported boundaries in this experiment;
individual predictions can still be lost or change class.
[Draft assemblies](draft-assemblies.md) has the full treatment; what follows is only
what changes about tuning.

**Read the first row of the audit before anything else.** `{sample}_ice_discarded.tsv`
opens with an `assembly_qc` row carrying the sample's own contig count and N50, plus one
sentence saying which of three measured bands it falls in. The band edges are 150 kb and
30 kb (`CONTIGUITY_N50_TRUSTED_BP` and `CONTIGUITY_N50_CLASS_ONLY_BP` in
`conjscan_to_ice.py`); 30 kb rather than 50 kb because the divider has
to sit *between* the two arms it separates, and detection held at the 50 kb arm.

**Do not turn on `require_trna_boundary_for_high`.** On the 52 complete records it
already demotes 23 of 30 high-confidence calls, because only 7 of them have
tRNA-anchored ends. On a short-read draft, where the flanks are usually simply absent, it
caps essentially everything, and a confidence column that always reads "medium" tells an
assessor nothing. Conservative does not mean "cap everything"; it means make the strong
claim rarer and better evidenced.

**Expect `min_reference_coverage` to block tier 4, and do not lower it to unblock it.**
A curated element split across contigs fails the coverage test by construction — already
visible at 12% and 49% on 5- and 6-contig assemblies. Admitting the fragment names a
transposon from half of it.

**Check `fraction_at_contig_boundary` in `{sample}_is_summary.tsv`.** A high value
indicates that fragmentation may affect the located IS count. With Illumina reads,
the optional copy-number analysis (`mobilome.isosdb.fasta_url`) provides an independent
abundance estimate. Its difference from the located count is not a direct measurement
of assembly collapse.

## Limits of parameter tuning { #what-no-knob-fixes }

These are limits of the design and of the data, and a methods section should say so.

**Relaxase signal in 70.6% of the tested IMEs.** Measured over all 395 curated IMEs at
production settings: 279 carry a relaxase that can be detected at all. The remaining
entries had no detectable relaxase signal at these settings. A shuffled-protein control
produced **0 of 395** hits at production settings. This measures the tested profiles on
this dataset; it is not a universal biological ceiling or a measure of complete-element
recovery. Size-matched, IMEs are detected as well as or better than ICEs in every size
bin; the earlier appearance otherwise was a size confound.

**Elements that carry no relaxase.** NBU1, NBU2, Tn*4555* and the MGI series are
mobilised in *trans* by a co-resident element, so they encode no relaxase because they do
not use one. This module is relaxase-anchored, which makes them out of scope by design
rather than misses. No threshold reaches them.

**Elements below about 2 kb.** 22 curated IMEs (5.6%) are under 2 kb, 17 of them a single
ORF. 18 of the 22 produce no relaxase signal at any threshold. Lowering
`--min-ime-element-bp` below 2,000 does not recover them; it admits noise.

**`spans_contigs` caps at low, absolutely.** MacSyFinder is run over the whole draft as
one pseudo-replicon, so genes it joined across a contig break may simply be unrelated
genes that happen to be adjacent in the file. Not a threshold, and not exposed.

**Boundaries on fragmented assemblies.** If the flanking sequence is not on the contig,
no parameter finds an *att* site. `boundary_method=none` is the correct answer, and the
interval is explicitly a floor.

## Traps

**Lowering a threshold changes what "detected" means.** A recall score that counts any
overlap flatters a looser setting, because one large call can cover several curated
elements and each counts as a hit. At `coverage_profile: 0.4` the internal scorer prints
16/18 ICEs rather than 15/18 — the extra "detection" is a 3,187 bp call clipping the
tail of a 161 kb element, 1.9% of it, and classed `ime` rather than `ice`. ICE recovery
with the expected class is 15/18 at every coverage value tested. Always report the
recovered fraction and the class next to the count.

**Whenever you loosen a threshold to see what you were missing, run the loosened setting
over shuffled or decoy input as well.** A permissive HMM pass finds a relaxase hit in
373 of 395 curated IMEs, which looks like a 94% ceiling until the shuffled control is
run at the same setting: **277 of 395** shuffled sequences also produce a hit. About 70%
of the shuffled controls produced a hit. This is a control hit rate, not an estimate
that 70% of the real calls are false. The comparison shows why additional permissive
hits cannot be treated as recovered elements without evaluation.

**A negative control must be inspected, not counted.** The negative set here included
*Bacillus subtilis* 168 on the grounds of having no curated ICEberg entry, and the
module called an ICE at 529,362–549,932. That is ICEBs1, discovered in *B. subtilis* 168
(Auchtung *et al.* 2005) — a true positive ICEberg has not curated, and it came close to
being recorded as a false alarm. An uncurated call may be a real element or a false
prediction; inspect the evidence before assigning either label.

**Say which arm and which build a before/after came from.** "15/18 ICEs, 5/12 IMEs, 2
calls in 32.6 Mb" mixes arms: the pilot figures come from runs with the ICEscan layer on,
while the "2 calls" figure comes from negative genomes that had no ICEscan directory. The
like-for-like negative figure for that configuration is **5**, not 2.

**Turning a parameter makes a run non-comparable with the last one.** The audit records the
value that was applied, but only per run, and a half-and-half result set is very hard to
unpick later. If you change a threshold mid-project, re-run everything — and keep the
`PROVENANCE.txt` each download rule writes, since the databases come from unversioned
addresses ([Enabling the module](enabling.md)).

## Related pages

| Page | Read it for |
|---|---|
| [Reading the output](output.md) | what every column of the mobility table means |
| [Draft assemblies](draft-assemblies.md) | the full fragmentation measurement and which columns survive it |
| [Enabling the module](enabling.md) | the master switch and the four optional layers |
| [The mobility ladder](mobility-ladder.md) | what each tier claims, and what it does not |
| [Validation](validation.md) | the benchmarks these numbers come from, and what they do not show |
| [Worked example](worked-example.md) | one genome end to end, with the "why not" trail |
| [Configuration](../reference/configuration.md) | every config key and its default, in the file's own order |
