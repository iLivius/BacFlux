# Reading the output

Everything lands in `08.mobilome/{sample}/`. Two files carry the result:

| File | What it is |
|---|---|
| `{sample}_amr_mobility.tsv` | **the deliverable** — one row per AMR gene, 46 columns |
| `{sample}_amr_mobility_audit.tsv` | **why** — one row per decision: every tier not raised, every structure rejected, every confidence cap |

Read them in this order: `{sample}_is_summary.tsv` first (how fragmented is this
assembly?), then the mobility table, then the audit for any row you intend to quote.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **AMR** | antimicrobial resistance |
    | **IS** | insertion sequence — the smallest mobile element, carrying only the genes it needs to move itself |
    | **IME** | integrative mobilisable element — integrated in the chromosome, mobilisable only with a helper element |

## The 46 columns

Coordinates are 1-based and inclusive on both ends, which is what AMRFinderPlus and
ISEScan both emit, so nothing is converted anywhere. A distance is the number of bases
strictly **between** two features — adjacent features are 0 bp apart. `NA` means the
question did not apply to this row.

### Where the gene is

| Column | Contents |
|---|---|
| `sample`, `contig` | which genome, which sequence |
| `replicon` | `chromosome` / `plasmid` / `unknown` |
| `replicon_id` | the plasmid's identifier when there is one |

### What the gene is

Straight from AMRFinderPlus, which is why these columns exist rather than ABRicate's.

| Column | Contents |
|---|---|
| `amr_gene`, `amr_name` | the element symbol and its full name |
| `amr_element_type` | `AMR` / `STRESS` / `VIRULENCE` — the module reasons about all three |
| `amr_class`, `amr_subclass` | drug class and subclass |
| `amr_start`, `amr_end`, `amr_strand` | position on the contig |
| `amrfinder_method` | `EXACTX`, `BLASTP`, `HMM`, `POINTX`, `PARTIAL_CONTIG_ENDX` … — a free confidence tier from the curator |
| `amr_pct_identity`, `amr_pct_coverage` | against the reference gene |
| `amr_partial_at_contig_end` | `yes` when the method says the hit runs off the end of the contig. Caps confidence |

### What is around it

| Column | Contents |
|---|---|
| `mge_context` | `none`, `is_adjacent`, `composite`, `unit_transposon`, `integron`, `ice`, `ime`, `plasmid` — plus the context-only types (`conjugative_region`, `genomic_island`, `aice`) and their `near_*` variants, which report a nearby element without raising the tier |
| `mge_id` | the element that set the context |
| `mge_name` | its curated name, when a naming layer supplied one |
| `named_element`, `named_element_type` | a curated transposon or integron **containing** the gene, kept even when a higher rung won the tier |
| `distance_bp` | the **measured** distance to that element. 0 when the gene is inside it |
| `orientation` | the element's strand relative to the gene's: `same` / `opposite` |

### The insertion-sequence evidence

| Column | Contents |
|---|---|
| `is_family`, `is_cluster` | ISEScan's family and cluster for the IS in question |
| `n_flanking_is` | IS copies on either side within the composite window |
| `same_orientation` | whether the two flanking copies point the same way. `NA` when a strand is unknown |
| `flanking_span_bp` | the measured span of a called composite — compare it against `mobilome.max_composite_span_bp` |
| `is_inside_amr_cds` | `yes` = an IS has landed **inside** the coding sequence. Read as likely **inactivation**, never as mobilisation |
| `is_amr_overlap_bp` | how much of the gene it overlaps |

### The conjugation machinery behind a tier 5 or 6 call

| Column | Contents |
|---|---|
| `plasmid_mobility` | `conjugative` / `mobilisable` / `non_mobilisable` / `unknown`, from Platon |
| `mobility_evidence` | Platon's own gene counts, passed through verbatim |
| `relaxase_type` | e.g. `MOBH`, `MOBF` — CONJscan's typing |
| `mpf_type` | the mating-pair apparatus type, e.g. `F`, `T` |
| `machinery_intact` | `no` = truncated or degraded |
| `machinery_source` | **read this one.** `containing_element:<id>` = the gene sits inside that element. `same_replicon:<id>` = the machinery is elsewhere on the same plasmid — still the right evidence, because a mating apparatus makes the whole replicon transferable, but it is not wrapped around this gene, and the column says so rather than implying containment |

!!! note "Platon counts genes; CONJscan types the machine"

    Platon counts hits to conjugation-related genes — cheap, always available, and not
    proof of a working apparatus. CONJscan asks whether a **complete, typed** mating-pair
    system is present. So Platon's replicon-level call sets the **tier**, and CONJscan
    corroborates or contradicts it, which moves the **confidence** and is named in the
    audit:

    | Platon | CONJscan | Result |
    |---|---|---|
    | conjugative | complete typed system | tier 6, **high** — `plasmid_conjugation_machinery_verified` |
    | conjugative | nothing, or relaxase only | tier 6, **medium** — `plasmid_conjugation_from_hit_counts_only` |
    | *not* conjugative | complete typed system | tier **unchanged**, confidence dropped — `conjscan_typed_system_but_platon_did_not` |

    That last row is deliberate. A disagreement is never silently promoted to "predicted
    self-transmissible": the tier follows the replicon-level call, the conflict is stated
    plainly, and a person decides.

### Element boundaries

| Column | Contents |
|---|---|
| `boundary_method` | `tRNA` / `denovo` / `none` / `NA` |
| `attL`, `attR` | the direct-repeat sequences, when a pair was found |

**Read `boundary_method` before you read any coordinate.**

| Value | Meaning | Does it move the element's edges? |
|---|---|---|
| `tRNA` | attL/attR at an annotated tRNA 3′ end — the arrangement integration actually produces | **yes** |
| `denovo` | a bracketing direct repeat with no tRNA behind it | **no** — reported as a lead only |
| `none` | no pair found; the interval is the machinery span and nothing more | no |
| `NA` | the gene has no ICE/IME context, so there were no boundaries to look for | — |

A de novo repeat never widens an element, and the restraint is measured: on a clinical
chromosome, 300 randomly placed 15 kb non-ICE spans produced a confident de novo
"boundary" 22% of the time — real chromosomes are full
of rRNA operons, REP elements and paralogues. Requiring the repeat to occur exactly twice
cut that to 16%, and refusing to *act* on one cuts what can reach the AMR table to 1%.
Widening an element makes every gene inside it cargo, so a fabricated 50 kb boundary
would turn a *gyrA* point mutation — the textbook non-transferable determinant — into
"predicted self-transmissible".

### How much the assembly can support

| Column | Contents |
|---|---|
| `contig_length` | the yardstick the next three are measured against |
| `dist_to_contig_end` | bases from the gene to the nearer end of its contig |
| `is_at_contig_boundary` | the context element runs into a contig end |
| `spans_contigs` | the element's evidence is on more than one contig. **An absolute cap at low confidence**, whatever else the row says |

### The verdict

| Column | Contents |
|---|---|
| `mobility_tier` | 1–6, or `NA` when the gene could not be placed |
| `mobility_tier_label` | the plain-language label ([The mobility ladder](mobility-ladder.md)) |
| `confidence` | `high` / `medium` / `low` |

Tier and confidence answer different questions and must be read together. `NA` /
`not_assessable` is not tier 1: "we could not assess this" and "we assessed it and found
nothing" are different results.

## The audit file

One row per **decision**, not per gene — a single gene can generate several: no context
found, plus two rejected structures, plus a confidence cap. Ten decision words, and they
are what you filter on:

```bash
cut -f8,9 {sample}_amr_mobility_audit.tsv | sort | uniq -c | sort -rn
```

| `decision` | Means |
|---|---|
| `input_missing` | a whole input table was absent or empty. Written once per missing table, **before any gene is looked at**, so "no evidence" can be told from "no data" |
| `input_assumed` | a table was present but missing an optional column, so a documented default was applied to all its rows |
| `row_skipped` | an AMRFinderPlus row had no contig coordinates, so it gets no tier |
| `no_mge_context` | the gene was assessed and nothing mobile was near it. A **result**, not a failure — this is what a tier 1 looks like |
| `tier_not_raised` | evidence was found but did not meet the bar for the next rung. The commonest word in the file |
| `evidence_rejected` | a candidate structure was thrown out (IS inside the coding sequence, unusable element type …) |
| `evidence_recorded` | kept as a note without changing the tier — e.g. a curated name that a higher rung outranked |
| `evidence_corroborated` | two independent tools agreed. Agreement is evidence and should be visible |
| `tools_disagree` | two tools disagreed. **Never resolved silently** |
| `confidence_capped` | the call was downgraded, and the reason says which guard fired |

`reason` is a fixed keyword you can filter on, and `detail` carries the numbers. A real example, for
*bla*NDM-1:

```text
tier_not_raised   flanking_is_pair_span_too_long
  1 candidate flanking pair(s) rejected. Examples:
  ...insertion_sequence-108774:110090 (IS110) / ...insertion_sequence-137590:140825 (new):
  span 32052 bp > 20000 bp limit
```

The composite detector found two flanking IS copies and rejected them, because 32 kb
is beyond `max_composite_span_bp`. That is a conservative miss rather than a silent one:
the measured span is in the file, so a reader who thinks 32 kb is defensible for their
organism can raise the threshold and re-run.

!!! note "The five ways a flanking pair fails"

    `flanking_is_pair_span_too_long`, `flanking_is_different_family`,
    `flanking_is_family_unknown`, `flanking_is_strand_unknown`,
    `flanking_is_pair_inverted_orientation`. Every AMR gene is tested against every
    nearby IS pair, so on an IS-rich clinical genome these are among the most common
    lines in the file — `flanking_is_different_family` was the second
    most common reason of all across a 12-genome clinical set.

## Everything else stage 08 writes

| File | What it is for |
|---|---|
| `{sample}_is_summary.tsv` | how many IS calls sit within `mobilome.contig_boundary_bp` of a contig end, and what fraction that is. Read it first |
| `{sample}_is_elements.tsv`, `_is_discarded.tsv` | one tidy row per insertion sequence, and the dropped rows with a reason |
| `{sample}_ice_candidates.tsv`, `_ice_discarded.tsv` | ICE/IME candidates with anchors, boundaries, class and confidence — and the discard trail. The first row of the discard file records this assembly's own contig count and N50, with a sentence saying which contiguity band it falls in |
| `{sample}_amrfinderplus.tsv`, `_amrfinderplus_mutations.tsv` | the raw AMR calls, and the point mutations separately |
| `{sample}_amrfinder_organism.txt`, `_amrfinder_organism_audit.tsv` | which `--organism` the GTDB-Tk call mapped to, or none, and why |
| `{sample}_replicon_calls.tsv` | chromosome or plasmid per contig, with the evidence |
| `{sample}_contig_lengths.tsv` | the yardstick behind every contig-edge flag |

!!! warning "`conjscan_output_missing` means the search did not run"

    A genome with no ICE or IME calls has usually been searched and found to carry no
    conjugation machinery, which is the normal result for an environmental isolate. But
    if the machinery search itself fails, the module records that and reports no
    calls — so the two look identical in the results table.

    They are distinguished in `{sample}_ice_discarded.tsv`. If it carries the reason
    `conjscan_output_missing`, nothing searched that genome: read the absence of ICE and
    IME calls as **unknown**, not as absence. This is the same distinction the mobility
    ladder makes between `not_assessable` and tier 1.

| `isescan/`, `conjscan/` (and `icescan/`) | the tools' own output trees, kept as raw evidence |

The optional layers add their own tables — curated names, ICE names, IS copy number — and
each is listed on [Optional layers](optional-layers.md).

## An ICE cannot sit on a plasmid

`{sample}_ice_candidates.tsv` reports elements on **every** replicon, because CONJscan
runs over the whole proteome. A conjugative element found on a plasmid is typed
`conjugative_region`, never `ice`, and that is a definition rather than a threshold: an
**I**ntegrative and **C**onjugative **E**lement integrates into a chromosome, and
something that is already its own replicon has nothing to integrate into. The classifier
enforces it through the integrase requirement:

| Integrase | Relaxase | Mating apparatus | Class |
|:-:|:-:|:-:|---|
| ✓ | ✓ | ✓ | **ICE** — predicted self-transmissible |
| ✓ | ✓ | ✗ | **IME** — mobilisable, needs a helper |
| ✓ | ✗ | ✗ | **CIME / genomic island** — integrated, no transfer machinery left. Passive: it moves only if something else carries it |
| ✗ | ✓ | ✓ | **conjugative region** — report it, do not call it an ICE |

A fifth class, **AICE** (the actinomycete class), is called only when the optional
ICEscan model set is on: an AICE conjugates, but it moves double-stranded DNA through
a translocase rather than nicking it with a relaxase, so the default models cannot see
it. It carries no tier, and [Validation](validation.md) reports it as unvalidated —
treat an AICE call as a hypothesis ([Optional layers](optional-layers.md)).

!!! warning "The class turns on one annotation, and an annotation is a string match"

    A `conjugative_region` row once appeared on the module's own positive control
    because Bakta labelled the integrase *DNA integration/recombination/inversion
    protein* and the product regex did not match that wording. The gene was 5.7 kb away
    and comfortably inside the clustering window. The element was under-called on the
    genome used to validate the module. Worth remembering before concluding
    that a region has no integrase.

## Next

- [Worked example](worked-example.md) — every tier, on one genome, with the real rows.
- [Draft assemblies](draft-assemblies.md) — which of these columns still hold on a draft.
- [Tuning](tuning.md) — symptom to knob.
- [Output files](../reference/output.md) — where stage 08 sits in the run.
