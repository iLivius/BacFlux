# The mobility ladder

Every AMR gene the module assesses gets exactly **one** tier, from 1 to 6. The tier
answers a single question — *if this gene can move, what would move it?* — and the
label beside the number is written into the output table, so nobody reading the TSV has
to look a number up.

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

## Which rungs are actually measured

| Tier | Status |
|:--:|---|
| 1 | benchmarked, and see the boundary caveat on [Validation](validation.md) |
| 2, 3 | exercised on real genomes; both are *inferences* from IS position and orientation |
| 4 | **a working code path that has never fired.** It needs the opt-in TnCentral layer *and* a curated element present at ≥80% of its reference length. No run kept on disk contains a tier 4 row |
| 5 | exercised by its plasmid route only. **No AMR gene in any run has been given an IME context**, so half of this rung is untested end to end |
| 6 | benchmarked, both routes |

Treat a tier 4 or tier 5 IME call as an unvalidated code path rather than a measured one,
and say so if it reaches a dossier. The numbers behind this table are on
[Validation](validation.md).

## Next

- [Reading the output](output.md) — the columns that carry the evidence behind a tier.
- [Turning it on](enabling.md) — the switch, and what runs.
- [Worked example](worked-example.md) — every tier on one real genome.
- [Draft assemblies](draft-assemblies.md) — what fragmentation does to a tier.
