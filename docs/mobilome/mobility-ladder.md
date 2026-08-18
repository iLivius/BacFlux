# The mobility ladder

Every AMR gene the module assesses gets exactly **one** tier, from 1 to 6. The tier
answers a single question — *if this gene can move, what would move it?* — and the
label beside the number is written into the output table, so nobody reading the TSV has
to look a number up.

!!! note "The five terms the table uses"

    An **IS** (insertion sequence) is the smallest mobile element: it encodes only what
    it needs to move itself and carries no passenger genes, so an IS can never *be* the
    resistance gene — only its neighbour.

    A **composite transposon** is what you get when two copies of the same IS land
    either side of something: the whole block, cargo included, can then move as a unit.

    An **integron** is a capture system with a cassette array; a resistance gene sitting
    in a cassette has an architecture somebody has characterised and named.

    An **ICE** (integrative and conjugative element) sits in the chromosome and carries
    its own conjugation machinery, so it can move itself into another cell.

    An **IME** (integrative mobilisable element) also sits in the chromosome and carries
    a **relaxase** — the enzyme that nicks the DNA to start a transfer — but no
    apparatus of its own, so it needs a helper element to move.

    That last distinction is the whole difference between tier 5 and tier 6.

| Tier | `mobility_tier_label` | Context | What it means |
|:--:|---|---|---|
| **1** | `intrinsic_candidate` | chromosomal, nothing mobile nearby | nothing was found that could move it |
| **2** | `expression_modulation_not_mobilisation` | an IS sits upstream, on the gene's strand | the IS can supply an outward-reading ("hybrid") promoter and raise expression. **Not mobilisation** — the gene still cannot move |
| **3** | `composite_mobilisable_within_cell` | the gene sits between two copies of the same IS family | a composite transposon: the whole block can hop within the cell |
| **4** | `named_element_mobilisable` | the gene sits inside a curated unit transposon or integron cassette | mobilisable, with an architecture somebody has characterised and named |
| **5** | `mobilisable_needs_helper` | on a plasmid, or inside an IME | it can reach another cell, but only with a helper element |
| **6** | `predicted_self_transmissible` | inside an ICE, or on a conjugative plasmid | **predicted** self-transmissible |

## How one tier is chosen

The rungs are tested **from the top down, and the first match wins**. Order matters
more than it looks, so it is written out here exactly as the code applies it:

1. inside an ICE → **6**
2. on a plasmid Platon typed conjugative → **6**
3. inside an IME → **5**
4. on any other plasmid → **5**
5. inside a curated named element → **4**
6. between two IS copies of one family → **3**
7. an IS within 500 bp upstream, same strand → **2**
8. none of the above → **1**

Two consequences follow. **A more specific call beats a broader one**: an ICE (an
interval) is preferred over a conjugative plasmid (a whole replicon) when both would
give tier 6, and the winning element is named in `mge_id`. And **a curated transposon on
a plasmid scores 5 or 6, not 4**, because the replicon evidence is tested first — which
is why tier 4 is reached only when a curated name is the *strongest* thing available,
in practice on the chromosome.

Where a lower rung's evidence existed but lost, it is not thrown away: it goes to the
audit file (`named_element_overrides_composite_pattern`, and so on), and a curated name
is still carried in `mge_name` even when a higher rung set the tier.

### The three windows that define the rungs

None of them is a property of a cell. Every row therefore reports the **measured**
distance next to the tier, so a reader who disagrees with a threshold can re-judge the
row without re-running anything.

| Window | Value | Used for |
|---|--:|---|
| upstream promoter distance | 500 bp | tier 2. IS*Ecp1* sits 42–266 bp upstream of *bla*CTX-M; beyond 500 bp an intervening gene is likely and the mechanism becomes speculative |
| context window | 5,000 bp | how close any IS must be before it is reported as `mge_context` at all |
| composite span | 20,000 bp (`mobilome.max_composite_span_bp`) | tier 3. Real composites run from ~2.5 kb to ~25 kb |

The only one that is a config key is the composite span. See [Tuning](tuning.md).

## Two things the ladder deliberately does not do

**A gene with no usable coordinates is `NA` / `not_assessable`, never tier 1.** "We
could not assess this" and "we assessed it and found nothing" are different results, and
collapsing them would quietly turn missing data into an intrinsic-resistance claim.

**An IS that has landed *inside* the coding sequence is reported separately.** That
usually inactivates the gene, so it gets its own `is_inside_amr_cds` flag rather than
being folded into a mobility call. An insertion is not a mobilisation.

## Intrinsic versus acquired — the framing this exists to serve

An AMR gene list tells you what an isolate can resist. It does not tell you whether that
resistance can move. A chromosomal efflux pump that every member of the species carries,
and a *bla*CTX-M on a conjugative plasmid, produce the same kind of line in an ABRicate
table and mean very different things.

That distinction is what the regulatory framing rests on, and the definitions are
explicit:

> **Intrinsic AMR gene** — "Gene inherent to strains of a bacterial species … An AMR
> gene is considered 'intrinsic' when it is shared by the vast majority of wild type
> strains of the same species (or subspecies) and is restricted to those located on the
> chromosome."
>
> **Acquired AMR gene** — "A resistance gene novel for the strain under assessment,
> acquired through horizontal transfer … Acquired AMR genes could be integrated in the
> bacterial chromosome **or** harboured on a separate genetic element."
>
> — EFSA Scientific Committee (2025), *Guidance on the characterisation of
> microorganisms in support of the risk assessment of products used in the food chain*,
> Glossary. <https://doi.org/10.2903/j.efsa.2025.9705>

The method for deciding the intrinsic side is published separately by
[EFSA BIOHAZ (2023)](https://doi.org/10.2903/j.efsa.2023.8323), together with a tool for
the automated analysis of gene distribution across a species.

!!! warning "Why tier 1 says *candidate*"

    Both halves of the definition are population-level, and **BacFlux implements no
    population comparison of any kind** — every analysis is per genome. Two things
    follow:

    - *Intrinsic* cannot be established from one genome. It is a claim about the
      species.
    - *Chromosomal* does not mean *intrinsic*. An acquired gene may perfectly well sit
      on the chromosome.

    So tier 1 is `intrinsic_candidate`, and the qualifier is doing real work: the module
    observed a chromosomal gene with no mobile-element context, which is *consistent
    with* intrinsic and does not establish it. The confirmatory work is a species-wide
    distribution analysis plus phenotypic testing.

The module produces **supporting evidence** for the intrinsic-versus-acquired judgement.
It does not produce the judgement, it is not "EFSA-compliant" and does not claim to be.
The full argument, with sources, is in
[`methods_amr_intrinsic_acquired.md`](https://github.com/iLivius/BacFlux/blob/main/docs/methods_amr_intrinsic_acquired.md).

### The one place a single genome speaks to the intrinsic side

Point mutations. A substitution in *gyrA*, *rpoB* or *rpsL* is a change to a gene the
strain already had — chromosomal and not transferable by definition. A BLAST screen
cannot see it: query a susceptible genome for *gyrA* and you get ~99.9% identity, query a
resistant one and you get ~99.9% identity, because the one residue that decides the
phenotype is averaged away inside the percentage.

AMRFinderPlus with `--organism` asks the other question — not "is *gyrA* here" but "is
residue 83 of GyrA a serine or a leucine" — and reports the hit with `Method: POINTX` or
`POINTP`. That is why the module maps the GTDB-Tk call onto an AMRFinderPlus organism
where it can, and why the mutations file is a separate deliverable
([Turning it on](enabling.md)).

An empty `{sample}_amrfinderplus_mutations.tsv` means **not assessed**, not "no mutations
found": AMRFinderPlus curates mutations for roughly three dozen mostly clinical
organisms, and the refusal, with its reason, is written to
`{sample}_amrfinder_organism_audit.tsv`.

## Language discipline

Tier 6 is `predicted_self_transmissible`, and the qualifier is never dropped. The
prediction rests on finding conjugation machinery that looks complete, not on watching a
transfer happen. **The confirmatory experiment is a filter or broth mating assay**, and
no output of this module substitutes for one.

## How far to trust a tier

Tiers 1, 2, 3 and 6 are the ones you will meet in practice, and they are the ones the
benchmark exercises. Two rungs are rarer than their position in the table suggests, and
it is worth knowing why before one turns up in a report.

**Tier 4 is structurally uncommon, not broken.** It needs three things at once: the
opt-in [TnCentral naming layer](optional-layers.md) configured, a curated transposon or
integron actually containing the gene, and the gene **not** on a plasmid and **not**
inside an ICE or IME — because those are tested first and the first match wins. A
resistance gene in a curated transposon on a plasmid scores tier 5, not tier 4.

Nothing is lost when that happens: the curated name is still written to
`named_element` and `named_element_type`, and the audit file records
`inside_named_element_but_higher_tier_applies`. So you can always see that a gene sat
in a named element even when its tier came from somewhere else.

It is also strict about evidence: a curated element must be present at ≥80% of its
reference length before it may name anything. That is why tier 4 turns up on closed
genomes much more readily than on drafts, where the same element often survives only as a
fragment. Worked through on a real genome, both halves of that — the override and the
80% refusal — are on [Worked example](worked-example.md).

**Tier 5 has two routes, and only one is well exercised.** The plasmid route — the gene
is on a plasmid — is common and benchmarked. The IME route, where a chromosomal gene
sits inside an integrative mobilisable element, is much rarer, and the benchmark has not
put an AMR gene inside one.

**What to do with that.** Tier 4 rests on a single measured instance, and the tier 5 IME
route on none. If either reaches something you are signing, check the evidence columns
yourself rather than the tier alone — `named_element`, `mge_id`, `boundary_method` and
`confidence` are all there for exactly this. The measured performance of the tiers that
*are* benchmarked is on [Validation](validation.md), including the boundary caveat that
applies to tier 1.

## Next

- [Reading the output](output.md) — the columns that carry the evidence behind a tier.
- [Turning it on](enabling.md) — the switch, and what runs.
- [Worked example](worked-example.md) — every tier on one real genome.
- [Draft assemblies](draft-assemblies.md) — what fragmentation does to a tier.
